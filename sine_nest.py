#!/usr/bin/env python3
"""Nested, split and clustered SINE copies; TinT-style chronology of SINE (sub)families.

  scan  --genome G.fa --copies COPIES.bed --library LIB.fa -o PREFIX [--window 400] [--close 300] [--cpu N]
        Searches every annotated copy and its surroundings (+-window) with nhmmer (HMMER) for SINE pieces
        of >= 20 nt, including fragments below the annotation thresholds, and classifies each locus:
          nested     a SINE (insert) sits between two pieces of an older SINE (host) whose consensus
                     coordinates join up; the host bases at the break are duplicated (TSD) on the far side
          split      two pieces of one copy whose consensus coordinates join up, nothing in between
          dimer      two near-full copies head to tail on one strand within 30 bp
          satellite  >= 3 pieces on one strand with a regular period and similar spacers (tandem array)
          close      other loci with several elements
          single     one element
        Writes PREFIX.elements.tsv (every piece, its role), PREFIX.loci.tsv, PREFIX.nest.tsv (one row per
        nesting event: TinT input), PREFIX.reassembled.bed (BED12 of hosts rebuilt from their pieces) and
        PREFIX.caution.bed (annotated copies whose flanks are SINE or that belong to a satellite: their
        flank-based orthology calls need caution).

  orth  --a SP=G.fa --b SP=G.fa --scan PREFIX -o ORTH.tsv [--sine-length L] [--threads N]
        Orthology for the elements of compound loci of genome a (nested, close, dimer: from scan PREFIX).
        The flanks *outside* the whole locus (unique sequence, not the host SINE) are placed in genome b as
        a concordant pair (bwa mem -a); in the b window between them each element is tested by its left
        and right junctions (present) against its empty-site junction with one TSD copy removed (absent).
        Writes rows in the orth_<a>-<b>.tsv format (status SINE/PM, cluster N<k>R) for sine_registry.py,
        plus ORTH.tsv.report with the outcome per element.

  supersede ORTH.tsv --scan SP=PREFIX [SP=PREFIX ...] -o OUT.tsv
        Main-pipeline rows (flank-based, clusters C...) with an insertion anchor inside a compound locus
        (nested, close, dimer, satellite) are moved to OUT.tsv.superseded: in such loci the consensus aligns
        to whichever element is most similar, so a flank-based row can pair an old host in one genome with a
        young insert in the other. The compound-locus rows (N...) of `orth` replace them.

  tint  NEST.tsv [NEST.tsv ...] -o PREFIX [--min-events 3] [--boot 200]
        TinT (Kriegs et al. 2007; Churakov et al. 2010) chronology from nesting events. Each type T (the
        library sequence that matched, e.g. a subfamily) is active around time mu_T; an element can only
        land in a host that already exists. With normally distributed activity (one period per type, no
        target preference -- the TinT assumptions) the probability that, in a nesting between A and B,
        A is the insert is Phi((mu_A - mu_B) / sqrt(2)) (Thurstone case V, sigma = 1 for every type).
        mu is fitted by maximum likelihood (first type fixed at 0); 95% intervals by parametric bootstrap.
        This is a reconstruction from the published model assumptions, not the original TinT code.
        Writes PREFIX.matrix.tsv (insert x host counts), PREFIX.chronology.tsv and PREFIX.chronology.svg.

Requires nhmmer (HMMER 3) for scan; Python standard library only otherwise.
"""
import argparse
import collections
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rescue_multi import Fasta, read_hits, local_align, revcomp  # noqa: E402

MIN_PIECE = 20          # TinT: elements shorter than 20 nt are ignored
JOIN_TOL = 30           # consensus coordinates of two host pieces may overlap (TSD) or gap by this much
ADJ_TOL = 30            # max genomic gap between an insert and the host pieces around it
FULL = 0.8              # near-full copy: >= 80% of the consensus


def die(msg):
    sys.exit(f'sine_nest.py: {msg}')


# ---------------------------------------------------------------------------------------------- scan
def read_library(path):
    lib, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]; lib[name] = []
            elif name:
                lib[name].append(line.upper())
    return {k: ''.join(v) for k, v in lib.items()}


def read_copies(path):
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 3 or f[0].startswith(('#', 'track')):
                continue
            name = f[3] if len(f) > 3 else f'{f[0]}:{f[1]}-{f[2]}'
            strand = f[5] if len(f) > 5 else '+'
            out.append((f[0], int(f[1]), int(f[2]), name, strand))
    out.sort()
    return out


def run_nhmmer(library, windows_fa, tbl, cpu):
    if not shutil.which('nhmmer'):
        die('nhmmer (HMMER 3) not found in PATH')
    cmd = ['nhmmer', '--dna', '--qformat', 'fasta', '--tformat', 'fasta', '--notextw', '--noali', '--cpu', str(cpu), '-E', '0.01', '--tblout', tbl, library, windows_fa]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def parse_tbl(tbl):
    hits = collections.defaultdict(list)
    with open(tbl) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split()
            win, query = f[0], f[2]
            hf, ht, af, at = map(int, f[4:8])
            strand = f[11]
            ev, score = float(f[12]), float(f[13])
            lo, hi = (af, at) if strand == '+' else (at, af)
            if hi - lo + 1 < MIN_PIECE:
                continue
            hits[win].append(dict(q=query, cf=hf, ct=ht, s=lo - 1, e=hi, strand=strand, ev=ev, score=score))
    return hits


def best_nonoverlapping(hits):
    """keep the best-scoring hit where hits overlap by > 50% of the shorter one"""
    keep = []
    for h in sorted(hits, key=lambda x: -x['score']):
        ok = True
        for k in keep:
            ov = min(h['e'], k['e']) - max(h['s'], k['s'])
            if ov > 0.5 * min(h['e'] - h['s'], k['e'] - k['s']):
                ok = False; break
        if ok:
            keep.append(h)
    keep.sort(key=lambda x: x['s'])
    # partial overlaps: an old host piece aligned by nhmmer often runs into the younger SINE inserted
    # next to it (both are SINE sequence); trim the weaker piece back to the boundary
    for a, b in zip(keep, keep[1:]):
        o = a['e'] - b['s']
        if o <= 0:
            continue
        w, side = (a, 'right') if a['score'] < b['score'] else (b, 'left')
        if side == 'right':
            w['e'] -= o
            if w['strand'] == '+': w['ct'] -= o
            else: w['cf'] += o
        else:
            w['s'] += o
            if w['strand'] == '+': w['cf'] += o
            else: w['ct'] -= o
    return [h for h in keep if h['e'] - h['s'] >= MIN_PIECE and h['ct'] >= h['cf']]


def host_order(a, b):
    """a, b pieces on one strand: return (5' piece, 3' piece) in host (consensus) orientation"""
    return (a, b) if a['strand'] == '+' else (b, a)


def joins(p5, p3):
    """consensus coordinates of the 5' and 3' piece continue each other (overlap = TSD, small gap)"""
    ov = p5['ct'] - p3['cf'] + 1
    return -JOIN_TOL <= ov <= JOIN_TOL, ov


def tsd_check(seq, ins_start, ins_end, ov=0):
    """direct repeat flanking the insert: the target bases before it recur right after it.
    nhmmer boundaries are fuzzy, so both ends may shift by up to 8 bp; the longest repeat
    (6-25 bp, <= 1 mismatch per 8 bp) wins"""
    best = ''
    for d1 in range(-8, 9):
        for d2 in range(-8, 9):
            a, b = ins_start + d1, ins_end + d2
            for t in range(25, len(best), -1):
                if t < 6:
                    break
                left = seq[a - t:a]; right = seq[b:b + t]
                if len(left) == t == len(right) and sum(x != y for x, y in zip(left, right)) <= t // 8:
                    best = left
                    break
    return best


def periodic(pieces):
    starts = [p['s'] for p in pieces]
    d = [b - a for a, b in zip(starts, starts[1:])]
    if len(d) < 2:
        return False
    m = sum(d) / len(d)
    sd = (sum((x - m) ** 2 for x in d) / len(d)) ** 0.5
    return sd < 0.1 * m


def spacer_similarity(seq, pieces, k=5):
    """mean 5-mer Jaccard similarity of consecutive spacers (padded by 10 bp): units of a tandem array
    share their spacer, independent insertions do not; insensitive to fuzzy piece boundaries"""
    sp = [seq[max(0, a['e'] - 10):b['s'] + 10] for a, b in zip(pieces, pieces[1:])]
    sp = [x for x in sp if len(x) >= 20]
    if len(sp) < 2:
        return 1.0
    km = [{x[i:i + k] for i in range(len(x) - k + 1)} for x in sp]
    sims = [len(x & y) / len(x | y) for x, y in zip(km, km[1:]) if x | y]
    return sum(sims) / len(sims) if sims else 0.0


def classify(pieces, seq, lib):
    """roles per piece index and the locus class; nested events"""
    roles = {i: 'single' for i in range(len(pieces))}
    events = []
    clen = {q: len(s) for q, s in lib.items()}
    full = [p['ct'] - p['cf'] + 1 >= FULL * clen[p['q']] for p in pieces]
    used = set()
    # nested: host piece - insert(s) - host piece
    for i in range(len(pieces)):
        for j in range(i + 2, min(len(pieces), i + 4)):
            a, b = pieces[i], pieces[j]
            if a['strand'] != b['strand'] or i in used or j in used:
                continue
            p5, p3 = host_order(a, b)
            ok, ov = joins(p5, p3)
            if not ok:
                continue
            mids = list(range(i + 1, j))
            gap_lo, gap_hi = a['e'], b['s']
            ins_lo, ins_hi = pieces[mids[0]]['s'], pieces[mids[-1]]['e']
            if ins_lo - gap_lo > ADJ_TOL or gap_hi - ins_hi > ADJ_TOL + 30:
                continue
            tsd = tsd_check(seq, ins_lo, ins_hi, ov)
            for m in mids:
                roles[m] = 'insert'; used.add(m)
            roles[i] = 'host5' if a is p5 else 'host3'; roles[j] = 'host3' if b is p3 else 'host5'
            used.update((i, j))
            events.append(dict(host=(i, j), inserts=mids, ov=ov, tsd=tsd))
    # split annotation: two joining pieces with (almost) nothing between them
    for i in range(len(pieces) - 1):
        a, b = pieces[i], pieces[i + 1]
        if i in used or i + 1 in used or a['strand'] != b['strand'] or b['s'] - a['e'] > 10:
            continue
        p5, p3 = host_order(a, b)
        ok, _ = joins(p5, p3)
        if ok:
            roles[i] = roles[i + 1] = 'split'; used.update((i, i + 1))
    rest = [i for i in range(len(pieces)) if i not in used]
    if len(rest) >= 3 and len({pieces[i]['strand'] for i in rest}) == 1:
        ps = [pieces[i] for i in rest]
        if periodic(ps) and spacer_similarity(seq, ps) > 0.3:
            for i in rest:
                roles[i] = 'unit'
            return roles, 'satellite', events
    if events:
        cls = 'nested'
    elif any(r == 'split' for r in roles.values()):
        cls = 'split'
    elif len(rest) == 2 and all(full[i] for i in rest) and pieces[rest[0]]['strand'] == pieces[rest[1]]['strand'] \
            and pieces[rest[1]]['s'] - pieces[rest[0]]['e'] <= ADJ_TOL:
        cls = 'dimer'
        for i in rest:
            roles[i] = 'dimer'
    elif len(pieces) > 1:
        cls = 'close'
        for i in rest:
            roles[i] = 'close'
    else:
        cls = 'single'
    return roles, cls, events


def identity(gseq, cons):
    """banded global identity of a piece to its consensus segment (p-distance proxy for age)"""
    n, m = len(gseq), len(cons)
    if not n or not m:
        return 0.0
    band = abs(n - m) + 30
    INF = -10 ** 9
    prev = {j: -j for j in range(0, min(m, band) + 1)}
    prev_m = {j: 0 for j in prev}
    for i in range(1, n + 1):
        cur, cur_m = {}, {}
        lo, hi = max(0, i - band), min(m, i + band)
        for j in range(lo, hi + 1):
            best, bm = INF, 0
            if j == 0:
                best, bm = -i, 0
            else:
                d = prev.get(j - 1)
                if d is not None:
                    s = 1 if gseq[i - 1] == cons[j - 1] else -1
                    if d + s > best:
                        best, bm = d + s, prev_m[j - 1] + (s == 1)
                l = cur.get(j - 1)
                if l is not None and l - 1 > best:
                    best, bm = l - 1, cur_m[j - 1]
            u = prev.get(j)
            if u is not None and u - 1 > best:
                best, bm = u - 1, prev_m[j]
            cur[j], cur_m[j] = best, bm
        prev, prev_m = cur, cur_m
    return 100.0 * prev_m.get(m, 0) / max(n, m)


def cmd_scan(args):
    lib = read_library(args.library)
    if not lib:
        die(f'no sequences in {args.library}')
    g = Fasta(args.genome)
    copies = read_copies(args.copies)
    # loci: copies closer than --close, as in SINE_orth_loc.bash
    loci = []
    for c in copies:
        if loci and loci[-1]['chrom'] == c[0] and c[1] - loci[-1]['end'] <= args.close:
            loci[-1]['copies'].append(c); loci[-1]['end'] = max(loci[-1]['end'], c[2])
        else:
            loci.append(dict(chrom=c[0], start=c[1], end=c[2], copies=[c]))
    tmp = tempfile.mkdtemp(prefix='sine_nest.')
    try:
        wfa = os.path.join(tmp, 'windows.fa'); tbl = os.path.join(tmp, 'hits.tbl')
        with open(wfa, 'w') as fh:
            for k, L in enumerate(loci, 1):
                L['id'] = f'L{k}'
                L['ws'] = max(0, L['start'] - args.window)
                L['we'] = min(g.length(L['chrom']), L['end'] + args.window)
                L['seq'] = g.fetch(L['chrom'], L['ws'], L['we'])
                fh.write(f">{L['id']}\n{L['seq']}\n")
        run_nhmmer(args.library, wfa, tbl, args.cpu)
        hits = parse_tbl(tbl)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    pre = args.output
    n_cls = collections.Counter(); n_ev = 0
    with open(f'{pre}.elements.tsv', 'w') as el, open(f'{pre}.loci.tsv', 'w') as lo, \
            open(f'{pre}.nest.tsv', 'w') as ne, open(f'{pre}.reassembled.bed', 'w') as rb, \
            open(f'{pre}.caution.bed', 'w') as cb:
        el.write('locus\telement\tchrom\tstart\tend\tstrand\ttype\tcons_from\tcons_to\tevalue\tscore\trole\tcopy\n')
        lo.write('locus\tchrom\tstart\tend\tclass\telements\tannotated_copies\n')
        ne.write('locus\tchrom\thost_type\tinsert_type\thost_identity\tinsert_identity\thost_break\t'
                 'tsd_length\ttsd\tinsert_orientation\thost\tinsert\n')
        for L in loci:
            pieces = best_nonoverlapping(hits.get(L['id'], []))
            if not pieces:
                lo.write(f"{L['id']}\t{L['chrom']}\t{L['start']}\t{L['end']}\tno_hit\t0\t{len(L['copies'])}\n")
                n_cls['no_hit'] += 1
                continue
            roles, cls, events = classify(pieces, L['seq'], lib)
            n_cls[cls] += 1
            gpos = lambda p: (L['ws'] + p['s'], L['ws'] + p['e'])
            def copy_of(p):
                s, e = gpos(p)
                for c in L['copies']:
                    if min(e, c[2]) - max(s, c[1]) > 0.5 * (c[2] - c[1]):
                        return c[3]
                return '.'
            for i, p in enumerate(pieces):
                s, e = gpos(p)
                el.write(f"{L['id']}\tE{i + 1}\t{L['chrom']}\t{s}\t{e}\t{p['strand']}\t{p['q']}\t{p['cf']}\t{p['ct']}\t"
                         f"{p['ev']:.2g}\t{p['score']:.1f}\t{roles[i]}\t{copy_of(p)}\n")
            lo.write(f"{L['id']}\t{L['chrom']}\t{L['ws'] + pieces[0]['s']}\t{L['ws'] + pieces[-1]['e']}\t{cls}\t"
                     f"{len(pieces)}\t{len(L['copies'])}\n")
            for evn in events:
                i, j = evn['host']; a, b = pieces[i], pieces[j]
                p5, p3 = host_order(a, b)
                host_seq = L['seq'][a['s']:a['e']] + L['seq'][b['s']:b['e']]
                if a['strand'] == '-':
                    host_seq = revcomp(L['seq'][b['s']:b['e']]) + revcomp(L['seq'][a['s']:a['e']])
                hcons = lib[p5['q']][p5['cf'] - 1:p3['ct']]
                hid = identity(host_seq, hcons)
                for m in evn['inserts']:
                    y = pieces[m]
                    yseq = L['seq'][y['s']:y['e']]
                    if y['strand'] == '-':
                        yseq = revcomp(yseq)
                    yid = identity(yseq, lib[y['q']][y['cf'] - 1:y['ct']])
                    orient = 'same' if y['strand'] == a['strand'] else 'opposite'
                    ne.write(f"{L['id']}\t{L['chrom']}\t{p5['q']}\t{y['q']}\t{hid:.1f}\t{yid:.1f}\t{p5['ct']}\t"
                             f"{len(evn['tsd'])}\t{evn['tsd'] or '-'}\t{orient}\t{L['ws'] + a['s']}-{L['ws'] + b['e']}\t"
                             f"{L['ws'] + y['s']}-{L['ws'] + y['e']}\n")
                    n_ev += 1
                s0, e0 = L['ws'] + a['s'], L['ws'] + b['e']
                rb.write(f"{L['chrom']}\t{s0}\t{e0}\t{L['id']}_host\t0\t{a['strand']}\t{s0}\t{e0}\t0,0,0\t2\t"
                         f"{a['e'] - a['s']},{b['e'] - b['s']},\t0,{b['s'] - a['s']},\n")
            for i, p in enumerate(pieces):
                c = copy_of(p)
                if c != '.' and roles[i] in ('insert', 'unit', 'host5', 'host3', 'dimer', 'close'):
                    s, e = gpos(p)
                    cb.write(f"{L['chrom']}\t{s}\t{e}\t{c}\t{roles[i]}\t{p['strand']}\n")
    print(f'{len(loci)} loci: ' + ', '.join(f'{k} {v}' for k, v in n_cls.most_common()) +
          f'; {n_ev} nesting events -> {pre}.nest.tsv', file=sys.stderr)


# ---------------------------------------------------------------------------------------------- orth
JUNCTION = 25          # bases on each side of an element boundary in a junction probe
SUPPORT_P = 0.75       # identity of both element junctions needed for "present"
SUPPORT_A = 0.75       # identity of the empty-site junction needed for "absent"
MARGIN = 0.1           # the winning hypothesis must beat the other by this much
ORTH_COLS = ['alignment', 'cluster', 'status', 'species1', 'locus1', 'sine1', 'species2', 'locus2', 'sine2',
             'LF', 'FL', 'LFcp', 'RFlength', 'FR', 'RFcp', 'OneTwo', 'OneSINE', 'TwoSINE', 'SL', 'SO', 'ST', 'FI',
             'anchor1', 'anchor2']


def semiglobal_identity(q, t):
    """identity of the whole probe q aligned to any part of t (free end gaps in t); returns (identity, end in t)"""
    n, m = len(q), len(t)
    prev = [(0, 0)] * (m + 1)                       # (score, matches); free leading gaps in t
    for i in range(1, n + 1):
        cur = [(-2 * i, 0)] + [None] * m
        for j in range(1, m + 1):
            d = prev[j - 1]; s_ = 1 if q[i - 1] == t[j - 1] else -1
            best = (d[0] + s_, d[1] + (s_ == 1))
            u = prev[j]
            if u[0] - 2 > best[0]:
                best = (u[0] - 2, u[1])
            l = cur[j - 1]
            if l[0] - 2 > best[0]:
                best = (l[0] - 2, l[1])
            cur[j] = best
        prev = cur
    j = max(range(m + 1), key=lambda x: prev[x][0])
    return prev[j][1] / n, j


def probe_score(probe, target):
    """identity of a junction probe at its best place in target (seeded, then semi-global)"""
    if len(probe) < 20:
        return 0.0, None
    r = local_align(probe, target, k=8, band=20)
    if not r:
        return 0.0, None
    _, ts, te, _ = r
    lo = max(0, ts - len(probe)); seg = target[lo:te + len(probe)]
    ident, end = semiglobal_identity(probe, seg)
    return ident, lo + end - len(probe) // 2


def window_locus(chrom, x, strand, slop, flank=300):
    """the pipeline's window string for an element whose 5' junction is at x (registry anchor = x)"""
    if strand == '+':
        return f'{chrom}:{max(0, x - flank)}-{x + slop}(+)'
    return f'{chrom}:{max(0, x - slop)}-{x + flank}(-)'


def cmd_orth(args):
    (sa, fa_a), (sb, fa_b) = args.a.split('=', 1), args.b.split('=', 1)
    ga, gb = Fasta(fa_a), Fasta(fa_b)
    loci = {r['locus']: r for r in _read_tsv(f'{args.scan}.loci.tsv') if r['class'] in ('nested', 'close', 'dimer')}
    elems = collections.defaultdict(list)
    for r in _read_tsv(f'{args.scan}.elements.tsv'):
        if r['locus'] in loci:
            elems[r['locus']].append(r)
    tsds = collections.defaultdict(list)
    for r in _read_tsv(f'{args.scan}.nest.tsv'):
        tsds[r['locus']].append(int(r['tsd_length']))
    slop = args.sine_length + args.flank
    tmp = tempfile.mkdtemp(prefix='sine_nest_orth.')
    try:
        fq = os.path.join(tmp, 'flanks.fa')
        with open(fq, 'w') as fh:
            for lid, L in loci.items():
                c, s0, e0 = L['chrom'], int(L['start']), int(L['end'])
                left = ga.fetch(c, s0 - args.flank, s0); right = ga.fetch(c, e0, e0 + args.flank)
                if len(left) >= 100 and len(right) >= 100:
                    fh.write(f'>{lid}|L\n{left}\n>{lid}|R\n{right}\n')
        sam = os.path.join(tmp, 'flanks.sam')
        with open(sam, 'w') as out:
            subprocess.run(['bwa', 'mem', '-a', '-t', str(args.threads), fa_b, fq], check=True, stdout=out,
                           stderr=subprocess.DEVNULL)
        hits = read_hits(sam)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    rows, report = [], collections.Counter()
    with open(args.output + '.report', 'w') as rep:
        rep.write('locus\telement\trole\toutcome\tleft_junction\tright_junction\tempty_site\n')
        for k, (lid, L) in enumerate(sorted(loci.items(), key=lambda x: int(x[0][1:])), 1):
            c, s0, e0 = L['chrom'], int(L['start']), int(L['end'])
            span = e0 - s0
            pairs = []
            for lh in hits.get(f'{lid}|L', []):
                for rh in hits.get(f'{lid}|R', []):
                    if lh[0] != rh[0] or lh[3] != rh[3]:
                        continue
                    if lh[3] == '+' and 0 <= rh[1] - lh[2] <= span + 3000:
                        pairs.append((lh[4] + rh[4], lh[0], lh[2], rh[1], '+'))
                    elif lh[3] == '-' and 0 <= lh[1] - rh[2] <= span + 3000:
                        pairs.append((lh[4] + rh[4], lh[0], rh[2], lh[1], '-'))
            pairs.sort(reverse=True)
            if not pairs:
                report['no_concordant_pair'] += 1
                rep.write(f'{lid}\t.\t.\tno_concordant_pair\t.\t.\t.\n'); continue
            if len(pairs) > 1 and pairs[1][0] >= 0.95 * pairs[0][0]:
                report['ambiguous'] += 1
                rep.write(f'{lid}\t.\t.\tambiguous:{len(pairs)}\t.\t.\t.\n'); continue
            _, bc, bs, be, bstr = pairs[0]
            ext = 40
            T = gb.fetch(bc, bs - ext, be + ext)
            if bstr == '-':
                T = revcomp(T)
            tb0 = bs - ext                          # genome-b coordinate of T[0] (+ orientation)
            def to_b(i):                            # T index -> genome b coordinate
                return tb0 + i if bstr == '+' else be + ext - i
            A = ga.fetch(c, s0 - 60, e0 + 60); a0 = s0 - 60
            test = [e for e in elems[lid] if e['role'] not in ('split', 'single', 'host5', 'host3')]
            hosts = [e for e in elems[lid] if e['role'] in ('host5', 'host3')]
            if hosts:                               # a host is tested as one element, from its outer ends
                test.append(dict(element='H' + '+'.join(h['element'] for h in hosts), role='host',
                                 start=str(min(int(h['start']) for h in hosts)), end=str(max(int(h['end']) for h in hosts)),
                                 strand=hosts[0]['strand']))
            for e in test:
                es, ee = int(e['start']) - a0, int(e['end']) - a0
                jl = A[es - JUNCTION:es + JUNCTION]; jr = A[ee - JUNCTION:ee + JUNCTION]
                sl, pl = probe_score(jl, T); sr, pr = probe_score(jr, T)
                best0, p0 = 0.0, None
                for t in sorted(set(tsds.get(lid, []) + list(range(0, 21)))):
                    j0 = A[es - 30:es] + A[ee + t:ee + t + 30]
                    sc, pos = probe_score(j0, T)
                    if sc > best0:
                        best0, p0 = sc, pos
                if min(sl, sr) >= SUPPORT_P and min(sl, sr) > best0 + MARGIN:
                    outcome, pos = 'present', (pl if e['strand'] == '+' else pr)
                elif best0 >= SUPPORT_A and best0 > max(sl, sr) + MARGIN:
                    outcome, pos = 'absent', p0
                else:
                    outcome, pos = 'unresolved', None
                report[f"{e['role']}:{outcome}"] += 1
                rep.write(f"{lid}\t{e['element']}\t{e['role']}\t{outcome}\t{sl:.2f}\t{sr:.2f}\t{best0:.2f}\n")
                if pos is None:
                    continue
                xa = int(e['start']) if e['strand'] == '+' else int(e['end'])
                xb = to_b(pos); bstrand = e['strand'] if bstr == '+' else ('-' if e['strand'] == '+' else '+')
                status = 'SINE' if outcome == 'present' else 'PM'
                name = f"N{k}_{e['element']}R"
                rows.append([f'{name}.nested.{status}', name, status, sa, window_locus(c, xa, e['strand'], slop, args.flank),
                             '1', sb, window_locus(bc, xb, bstrand, slop, args.flank), '1' if outcome == 'present' else '0']
                            + [''] * (len(ORTH_COLS) - 11) + [xa, xb])
    with open(args.output, 'w') as fh:
        fh.write('\t'.join(ORTH_COLS) + '\n')
        for r in rows:
            fh.write('\t'.join(map(str, r)) + '\n')
    print(f'{len(loci)} compound loci of {sa} in {sb}: ' + ', '.join(f'{k} {v}' for k, v in sorted(report.items())) +
          f' -> {args.output}', file=sys.stderr)


def _read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        return [dict(zip(header, l.rstrip('\n').split('\t'))) for l in fh]


# ---------------------------------------------------------------------------------------------- supersede
def cmd_supersede(args):
    comp = collections.defaultdict(list)            # (species, chrom) -> [(start, end)]
    for spec in args.scan:
        sp, pre = spec.split('=', 1)
        for r in _read_tsv(f'{pre}.loci.tsv'):
            if r['class'] in ('nested', 'close', 'dimer', 'satellite'):
                comp[(sp, r['chrom'])].append((int(r['start']) - 30, int(r['end']) + 30))
    rows = _read_tsv(args.orth)
    # tables made before the anchor columns: anchor from the window and the most common window length
    lens = collections.Counter()
    for r in rows:
        for k in ('locus1', 'locus2'):
            m = re.match(r'^(.*):(\d+)-(\d+)\(([+-])\)$', r[k])
            if m:
                lens[int(m.group(3)) - int(m.group(2))] += 1
    slop = lens.most_common(1)[0][0] - 300 if lens else 0
    def inside(sp, locus, anchor):
        m = re.match(r'^(.*):(\d+)-(\d+)\(([+-])\)$', locus)
        if not m:
            return False
        chrom = m.group(1)
        if anchor not in ('', None):
            a = int(anchor)
        else:
            a = int(m.group(3)) - slop if m.group(4) == '+' else int(m.group(2)) + slop
        return any(s0 <= a <= e0 for s0, e0 in comp.get((sp, chrom), ()))
    with open(args.orth) as fh:
        header = fh.readline()
    kept, dropped = [], []
    for r in rows:
        main = not r['cluster'].startswith('N')
        if main and (inside(r['species1'], r['locus1'], r.get('anchor1')) or inside(r['species2'], r['locus2'], r.get('anchor2'))):
            dropped.append(r)
        else:
            kept.append(r)
    cols = header.rstrip('\n').split('\t')
    for path, rs in ((args.output, kept), (args.output + '.superseded', dropped)):
        with open(path, 'w') as fh:
            fh.write(header)
            for r in rs:
                fh.write('\t'.join(r.get(c, '') for c in cols) + '\n')
    print(f'{len(dropped)} flank-based rows inside compound loci moved to {args.output}.superseded; '
          f'{len(kept)} rows kept', file=sys.stderr)


# ---------------------------------------------------------------------------------------------- tint
def Phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def fit_mu(types, n, iters=2000, lr=0.05):
    """maximise sum n[a][b] * log Phi((mu_a - mu_b)/sqrt 2) over mu (mu[0] = 0) by gradient ascent"""
    mu = [0.0] * len(types)
    k = len(types)
    s2 = math.sqrt(2)
    for it in range(iters):
        grad = [0.0] * k
        for a in range(k):
            for b in range(k):
                if a == b or not n[a][b]:
                    continue
                d = (mu[a] - mu[b]) / s2
                p = max(Phi(d), 1e-12)
                dens = math.exp(-d * d / 2) / math.sqrt(2 * math.pi)
                g = n[a][b] * dens / p / s2
                grad[a] += g; grad[b] -= g
        tot = sum(sum(r) for r in n) or 1
        step = lr * 50 / tot
        for a in range(1, k):
            mu[a] += step * grad[a] - 1e-3 * mu[a] * step   # weak ridge keeps all-one-way pairs finite
        if max(abs(x) for x in grad[1:]) * step < 1e-7:
            break
    return mu


def cmd_tint(args):
    counts = collections.Counter(); types = collections.Counter()
    for path in args.nest:
        with open(path) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            for line in fh:
                r = dict(zip(header, line.rstrip('\n').split('\t')))
                counts[(r['insert_type'], r['host_type'])] += 1
                types[r['insert_type']] += 1; types[r['host_type']] += 1
    ts = sorted(t for t, c in types.items() if c >= args.min_events)
    if len(ts) < 2:
        die(f'fewer than two types with >= {args.min_events} nesting events')
    ix = {t: i for i, t in enumerate(ts)}
    n = [[0] * len(ts) for _ in ts]
    for (a, b), c in counts.items():
        if a in ix and b in ix:
            n[ix[a]][ix[b]] += c
    mu = fit_mu(ts, n)
    # parametric bootstrap
    rng = random.Random(1)
    boots = [[] for _ in ts]
    for _ in range(args.boot):
        nb = [[0] * len(ts) for _ in ts]
        for a in range(len(ts)):
            for b in range(a + 1, len(ts)):
                m = n[a][b] + n[b][a]
                if not m:
                    continue
                p = Phi((mu[a] - mu[b]) / math.sqrt(2))
                x = sum(rng.random() < p for _ in range(m))
                nb[a][b], nb[b][a] = x, m - x
        mb = fit_mu(ts, nb)
        for i, v in enumerate(mb):
            boots[i].append(v)
    pre = args.output
    with open(f'{pre}.matrix.tsv', 'w') as fh:
        fh.write('insert\\host\t' + '\t'.join(ts) + '\tas_insert\n')
        for a in range(len(ts)):
            fh.write(ts[a] + '\t' + '\t'.join(map(str, n[a])) + f'\t{sum(n[a])}\n')
        fh.write('as_host\t' + '\t'.join(str(sum(n[a][b] for a in range(len(ts)))) for b in range(len(ts))) + '\n')
    order = sorted(range(len(ts)), key=lambda i: mu[i])
    ref = order[0]                                   # times relative to the oldest type, in every replicate
    rows = []
    for i in order:
        d = sorted(boots[i][r] - boots[ref][r] for r in range(len(boots[i])))
        lo = d[int(0.025 * len(d))] if d else 0.0; hi = d[max(0, int(0.975 * len(d)) - 1)] if d else 0.0
        younger = [ts[b] for b in range(len(ts)) if n[i][b] and not n[b][i]]
        older = [ts[b] for b in range(len(ts)) if n[b][i] and not n[i][b]]
        both = [ts[b] for b in range(len(ts)) if n[i][b] and n[b][i]]
        note = ('two-way' if both else
                'one-way: only inserted into ' + ','.join(younger) + ' (time is a lower bound)' if younger and not older else
                'one-way: only hosts ' + ','.join(older) + ' (time is an upper bound)' if older and not younger else 'one-way')
        rows.append((ts[i], mu[i] - mu[ref], lo, hi, sum(n[i]), sum(r[i] for r in n), note))
    with open(f'{pre}.chronology.tsv', 'w') as fh:
        fh.write('type\tmu\tci95_low\tci95_high\tas_insert\tas_host\tcomparisons\n')
        for r in rows:
            fh.write(f'{r[0]}\t{r[1]:.3f}\t{r[2]:.3f}\t{r[3]:.3f}\t{r[4]}\t{r[5]}\t{r[6]}\n')
    # SVG: one bar per type, oldest on top; bar = 95% CI, dot = mu, ellipse = +-0.674 sd (75% of activity)
    lo_all = min(r[2] for r in rows) - 1.5; hi_all = max(r[3] for r in rows) + 1.5
    W, H, top = 640, 30 * len(rows) + 50, 30
    X = lambda v: 120 + (v - lo_all) / (hi_all - lo_all) * (W - 140)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="sans-serif" font-size="12">',
           f'<text x="120" y="18">older  &#8594;  younger (relative TinT time, activity sd = 1)</text>']
    for k, r in enumerate(rows):
        y = top + 30 * k + 10
        svg.append(f'<text x="8" y="{y + 4}">{r[0]}</text>')
        svg.append(f'<line x1="{X(r[1] - 1.96):.1f}" y1="{y}" x2="{X(r[1] + 1.96):.1f}" y2="{y}" stroke="#999"/>')
        svg.append(f'<ellipse cx="{X(r[1]):.1f}" cy="{y}" rx="{X(r[1] + 0.674) - X(r[1]):.1f}" ry="7" fill="#9ecae1" stroke="#3182bd"/>')
        svg.append(f'<line x1="{X(r[2]):.1f}" y1="{y - 9}" x2="{X(r[3]):.1f}" y2="{y - 9}" stroke="#d62728" stroke-width="2"/>')
    svg.append('</svg>')
    with open(f'{pre}.chronology.svg', 'w') as fh:
        fh.write('\n'.join(svg) + '\n')
    print('chronology (oldest first): ' + ' < '.join(r[0] for r in rows) +
          f' -> {pre}.chronology.tsv .matrix.tsv .chronology.svg', file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('scan', help='find nested, split, dimeric, satellite and close SINE copies')
    s.add_argument('--genome', required=True, help='genome FASTA with .fai')
    s.add_argument('--copies', required=True, help='annotated copies (BED, e.g. sear2k)')
    s.add_argument('--library', required=True, help='consensus sequence(s): one family or its subfamilies')
    s.add_argument('-o', '--output', required=True)
    s.add_argument('--window', type=int, default=400, help='searched sequence around each locus [400]')
    s.add_argument('--close', type=int, default=300, help='copies closer than this form one locus [300]')
    s.add_argument('--cpu', type=int, default=2)
    s.set_defaults(func=cmd_scan)
    o = sub.add_parser('orth', help='orthology of the elements of compound loci (outer flanks + junction tests)')
    o.add_argument('--a', required=True, metavar='SP=FASTA', help='genome with the compound loci (scan PREFIX)')
    o.add_argument('--b', required=True, metavar='SP=FASTA', help='genome searched (bwa index required)')
    o.add_argument('--scan', required=True, metavar='PREFIX', help='scan output prefix of genome a')
    o.add_argument('-o', '--output', required=True)
    o.add_argument('--sine-length', type=int, default=360, help='consensus length, for window strings [360]')
    o.add_argument('--flank', type=int, default=300)
    o.add_argument('--threads', type=int, default=2)
    o.set_defaults(func=cmd_orth)
    u = sub.add_parser('supersede', help='move flank-based rows inside compound loci aside')
    u.add_argument('orth', help='orth_<a>-<b>.tsv (with anchor1/anchor2 columns)')
    u.add_argument('--scan', nargs='+', required=True, metavar='SP=PREFIX')
    u.add_argument('-o', '--output', required=True)
    u.set_defaults(func=cmd_supersede)
    t = sub.add_parser('tint', help='TinT-style chronology from nesting events')
    t.add_argument('nest', nargs='+', help='PREFIX.nest.tsv file(s) from scan')
    t.add_argument('-o', '--output', required=True)
    t.add_argument('--min-events', type=int, default=3, help='types with fewer nesting events are left out [3]')
    t.add_argument('--boot', type=int, default=200, help='bootstrap replicates for the intervals [200]')
    t.set_defaults(func=cmd_tint)
    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
