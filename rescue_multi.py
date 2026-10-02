#!/usr/bin/env python3
"""Resolve multi-copy clusters of SINE_orth_loc with both flanks, and account for every cluster.

Clusters of 3-10 loci that the multi stage could not reduce to one pair, and clusters
of more than 10 loci, are usually caused by a repetitive left flank: the left flank of
a SINE copy matches several places in the other genome. The right flank is often unique.

  prepare  For every annotated copy (the *_uniq.bed copies the run started from) that lies
           in an unresolved cluster, write its left and right 300-bp flanks (SINE
           orientation) as FASTA, one file per genome, to be mapped with bwa mem -a.
  resolve  From the SAM files, take the best hits of each flank as anchors and align the
           other flank in the window next to each anchor where it must lie (mate rescue:
           same strand, right flank 3' of the left one at the distance of an empty site or
           of a SINE). Accept the placement whose two flanks score clearly better than any
           other placement. Writes the accepted pairs as cluster lines
           (M<n>R) in the same format as <sp1>-<sp2>_double, for alignment and ComPair.sh.
  account  One row per cluster of the run: class (double/multi/poly), number of loci and
           what became of it, including the outcome of the rescue.

Only the Python standard library is used.
"""

import argparse
import collections
import os
import re
import sys

FLANK = 300
MIN_ALIGNED = 100  # same as the pipeline's filter on mapped flanks


def die(msg):
    sys.exit(f'rescue_multi.py: {msg}')


class Fasta:
    """Random access to a FASTA file through its samtools .fai index."""

    def __init__(self, path):
        self.fh = open(path, 'rb')
        self.idx = {}
        with open(path + '.fai') as fai:
            for line in fai:
                name, length, offset, linebases, linewidth = line.split('\t')[:5]
                self.idx[name] = (int(length), int(offset), int(linebases), int(linewidth))

    def length(self, chrom):
        return self.idx[chrom][0]

    def fetch(self, chrom, start, end):
        length, offset, lb, lw = self.idx[chrom]
        start, end = max(0, start), min(length, end)
        if end <= start:
            return ''
        a = offset + (start // lb) * lw + start % lb
        b = offset + ((end - 1) // lb) * lw + (end - 1) % lb + 1
        self.fh.seek(a)
        return self.fh.read(b - a).decode().replace('\n', '').replace('\r', '').upper()


def revcomp(s):
    return s[::-1].translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))


def kv(specs, what):
    out = {}
    for spec in specs:
        k, _, v = spec.partition('=')
        if not v:
            die(f'{what} expects NAME=PATH, got "{spec}"')
        out[k] = v
    return out


LOCUS = re.compile(r'>([^>]*?):(\d+)-(\d+)\(([+-])\)')


def read_clusters(path):
    """Cluster lines: 'CnR >chrom:s-e(strand)::info SeqStartSEQ,>...' -> {cluster: [(chrom, s, e, strand)]}."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            name = line.split(' ', 1)[0]
            out[name] = [(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)) for m in LOCUS.finditer(line)]
    return out


def resolved_clusters(statpairs):
    done = set()
    if statpairs and os.path.exists(statpairs):
        with open(statpairs) as fh:
            for line in fh:
                done.add(line.split('\t', 1)[0].split('.', 1)[0])
    return done


def read_bed(path, sp):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 6:
                rows.append((sp, f[0], int(f[1]), int(f[2]), f[5]))
    return rows


def copy_id(c):
    sp, chrom, s, e, strand = c
    return f'{sp}|{chrom}:{s}-{e}({strand})'


def parse_copy_id(cid):
    sp, locus = cid.split('|', 1)
    m = re.match(r'^(.*):(\d+)-(\d+)\(([+-])\)$', locus)
    return sp, m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


# ── prepare ────────────────────────────────────────────────────────────────────

def cmd_prepare(args):
    genomes = {k: Fasta(v) for k, v in kv(args.genome, '--genome').items()}
    copies = kv(args.copies, '--copies')
    if len(genomes) != 2 or set(genomes) != set(copies):
        die('give --genome and --copies for the same two species')
    done = resolved_clusters(args.resolved)
    clusters = {}
    for path in args.clusters:
        clusters.update(read_clusters(path))
    open_clusters = {c: loci for c, loci in clusters.items() if c not in done}

    windows = collections.defaultdict(list)          # chrom -> [(start, end, cluster)]
    for cl, loci in open_clusters.items():
        for chrom, s, e, _ in loci:
            windows[chrom].append((s, e, cl))
    for v in windows.values():
        v.sort()

    selected = collections.defaultdict(list)          # species -> [(copy, cluster)]
    for sp, path in copies.items():
        for c in read_bed(path, sp):
            _, chrom, s, e, _ = c
            for ws, we, cl in windows.get(chrom, ()):
                if ws >= e:
                    break
                if we > s:
                    selected[sp].append((c, cl))
                    break

    with open(f'{args.out}.copies.tsv', 'w') as tsv:
        tsv.write('copy\tcluster\n')
        for sp, g in genomes.items():
            with open(f'{args.out}_{sp}.fa', 'w') as fa:
                for c, cl in selected[sp]:
                    _, chrom, s, e, strand = c
                    if strand == '+':
                        left, right = g.fetch(chrom, s - FLANK, s), g.fetch(chrom, e, e + FLANK)
                    else:
                        left, right = revcomp(g.fetch(chrom, e, e + FLANK)), revcomp(g.fetch(chrom, s - FLANK, s))
                    cid = copy_id(c)
                    if len(left) >= MIN_ALIGNED:
                        fa.write(f'>{cid}|L\n{left}\n')
                    if len(right) >= MIN_ALIGNED:
                        fa.write(f'>{cid}|R\n{right}\n')
                    tsv.write(f'{cid}\t{cl}\n')
    n = sum(len(v) for v in selected.values())
    print(f'{len(open_clusters)} unresolved clusters, {n} copies to rescue', file=sys.stderr)


# ── resolve ────────────────────────────────────────────────────────────────────

CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


def read_hits(sam):
    """{read name: [(chrom, start, end, strand, score)]} for primary and secondary hits."""
    hits = collections.defaultdict(list)
    with open(sam) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            f = line.split('\t')
            flag = int(f[1])
            if flag & 4 or flag & 2048:
                continue
            ref_len = q_len = 0
            for n, op in CIGAR.findall(f[5]):
                if op in 'MDN=X':
                    ref_len += int(n)
                if op in 'MI=X':
                    q_len += int(n)
            if q_len < MIN_ALIGNED:
                continue
            score = next((int(t[5:]) for t in f[11:] if t.startswith('AS:i:')), 0)
            start = int(f[3]) - 1
            hits[f[0]].append((f[2], start, start + ref_len, '-' if flag & 16 else '+', score))
    return hits


def local_align(q, t, k=11, band=25):
    """Best local alignment of q in t (match 1, mismatch -4, gap -6), banded around the
    diagonal with most shared k-mers. Returns (score, t_start, t_end, q_aligned) or None."""
    if len(q) < k or len(t) < k:
        return None
    index = collections.defaultdict(list)
    for i in range(len(t) - k + 1):
        index[t[i:i + k]].append(i)
    diag = collections.Counter()
    for j in range(len(q) - k + 1):
        for i in index.get(q[j:j + k], ()):
            diag[i - j] += 1
    if not diag:
        return None
    d, votes = diag.most_common(1)[0]
    if votes < 3:
        return None
    best = (0, 0, 0, 0)
    prev = {}
    start_of = {}
    for j in range(len(q)):
        cur, cur_start = {}, {}
        lo, hi = max(0, j + d - band), min(len(t), j + d + band + 1)
        for i in range(lo, hi):
            m = 1 if q[j] == t[i] else -4
            diag_s = prev.get(i - 1, 0) + m
            up = prev.get(i, 0) - 6
            left = cur.get(i - 1, 0) - 6
            sc = max(0, diag_s, up, left)
            if sc == 0:
                continue
            cur[i] = sc
            if sc == diag_s:
                cur_start[i] = start_of.get(i - 1, (i, j)) if prev.get(i - 1, 0) > 0 else (i, j)
            elif sc == up:
                cur_start[i] = start_of.get(i, (i, j))
            else:
                cur_start[i] = cur_start.get(i - 1, (i, j))
            if sc > best[0]:
                ts, qs = cur_start[i]
                best = (sc, ts, i + 1, j + 1 - qs)
        prev, start_of = cur, cur_start
    return best if best[0] > 0 else None


def cmd_resolve(args):
    genomes = {k: Fasta(v) for k, v in kv(args.genome, '--genome').items()}
    sams = kv(args.sam, '--sam')
    sp_list = list(genomes)
    other = {sp_list[0]: sp_list[1], sp_list[1]: sp_list[0]}
    slop = args.sine_length + FLANK
    max_gap = 2 * args.sine_length + 100

    clusters = {}
    with open(f'{args.prep}.copies.tsv') as fh:
        fh.readline()
        for line in fh:
            cid, cl = line.rstrip('\n').split('\t')
            clusters[cid] = cl
    flanks = {}
    for sp in genomes:
        path = f'{args.prep}_{sp}.fa'
        if os.path.exists(path):
            with open(path) as fh:
                name = None
                for line in fh:
                    if line.startswith('>'):
                        name = line[1:].strip()
                    else:
                        flanks[name] = line.strip()

    def place(g, flank, chrom, lo, hi, strand):
        """Align a flank (SINE orientation) to target[lo:hi] on `strand`; genome coordinates back."""
        lo, hi = max(0, lo), min(g.length(chrom), hi)
        t = g.fetch(chrom, lo, hi)
        if strand == '-':
            t = revcomp(t)
        r = local_align(flank, t)
        if not r or r[3] < MIN_ALIGNED:
            return None
        sc, ts, te, _ = r
        if strand == '+':
            return sc, lo + ts, lo + te
        return sc, hi - te, hi - ts

    report, pairs = [], []
    for sp, sam in sams.items():
        if not os.path.exists(sam):
            continue
        hits = read_hits(sam)
        tgt = other[sp]
        g = genomes[tgt]
        for cid, cl in clusters.items():
            if not cid.startswith(sp + '|'):
                continue
            fl, fr = flanks.get(cid + '|L', ''), flanks.get(cid + '|R', '')
            left, right = hits.get(cid + '|L', []), hits.get(cid + '|R', [])
            # anchors: the best hits of each flank; the other flank is searched next to each
            anchors = ([('L', h) for h in sorted(left, key=lambda h: -h[4])[:args.anchors]] +
                       [('R', h) for h in sorted(right, key=lambda h: -h[4])[:args.anchors]])
            placements = []
            for side, (chrom, hs, he, strand, _) in anchors:
                if side == 'L':
                    lp = place(g, fl, chrom, hs - 20, he + 20, strand)
                    if not lp:
                        continue
                    if strand == '+':
                        rp = place(g, fr, chrom, lp[2] - 60, lp[2] + max_gap + FLANK, '+')
                    else:
                        rp = place(g, fr, chrom, lp[1] - max_gap - FLANK, lp[1] + 60, '-')
                else:
                    rp = place(g, fr, chrom, hs - 20, he + 20, strand)
                    if not rp:
                        continue
                    if strand == '+':
                        lp = place(g, fl, chrom, rp[1] - max_gap - FLANK, rp[1] + 60, '+')
                    else:
                        lp = place(g, fl, chrom, rp[2] - 60, rp[2] + max_gap + FLANK, '-')
                if not lp or not rp:
                    continue
                gap = rp[1] - lp[2] if strand == '+' else lp[1] - rp[2]
                if -60 <= gap <= max_gap:
                    placements.append((lp[0] + rp[0], chrom, lp[1], lp[2], strand, gap))
            best_at = []
            for c in sorted(placements, reverse=True):   # one per placement of the left flank
                if not any(c[1] == o[1] and c[4] == o[4] and abs(c[2] - o[2]) <= 20 for o in best_at):
                    best_at.append(c)
            ranked = best_at
            if not left and not right:
                outcome = 'no_hit'
            elif not ranked:
                outcome = 'no_concordant_pair'
            elif len(ranked) > 1 and ranked[0][0] - ranked[1][0] < args.margin:
                outcome = f'ambiguous:{len(ranked)}'
            else:
                outcome = 'resolved'
            report.append((cid, cl, outcome, len(left), len(right), len(ranked)))
            if outcome != 'resolved':
                continue
            _, chrom, ls, le, strand, gap = ranked[0]
            tlen = g.length(chrom)
            t_win = (chrom, ls, min(tlen, le + slop), strand) if strand == '+' else (chrom, max(0, ls - slop), le, strand)
            _, qchrom, qs, qe, qstrand = parse_copy_id(cid)
            qlen = genomes[sp].length(qchrom)
            q_win = ((qchrom, max(0, qs - FLANK), min(qlen, qs + slop), '+') if qstrand == '+'
                     else (qchrom, max(0, qe - slop), min(qlen, qe + FLANK), '-'))
            pairs.append((sp, q_win, tgt, t_win, cid, cl, gap))

    def anchor(w):
        chrom, s, e, strand = w
        return chrom, (e - slop if strand == '+' else s + slop)

    def outcome_of(cid, text):
        i = next(i for i, r in enumerate(report) if r[0] == cid)
        report[i] = report[i][:2] + (text,) + report[i][3:]

    seen = []
    written = 0
    with open(args.out, 'w') as out:
        for sp, q_win, tgt, t_win, cid, cl, gap in pairs:
            a_q, a_t = anchor(q_win), anchor(t_win)
            dup = None
            for name, s_q, s_t in seen:   # the same pair found from both genomes
                if ((s_q[0] == a_q[0] and abs(s_q[1] - a_q[1]) <= 60 and s_t[0] == a_t[0] and abs(s_t[1] - a_t[1]) <= 60) or
                        (s_q[0] == a_t[0] and abs(s_q[1] - a_t[1]) <= 60 and s_t[0] == a_q[0] and abs(s_t[1] - a_q[1]) <= 60)):
                    dup = name
                    break
            if dup:
                outcome_of(cid, f'resolved:{dup}')
                continue
            written += 1
            name = f'M{written}R'
            seen.append((name, a_q, a_t))
            loci = []
            for g, (chrom, s, e, strand) in ((genomes[sp], q_win), (genomes[tgt], t_win)):
                seq = g.fetch(chrom, s, e)
                if strand == '-':
                    seq = revcomp(seq)
                loci.append(f'>{chrom}:{s}-{e}({strand})::RESCUE={cl}SeqStart{seq}')
            out.write(f'{name} ' + ','.join(loci) + '\n')
            outcome_of(cid, f'resolved:{name}')

    with open(args.report, 'w') as rep:
        rep.write('copy\tcluster\toutcome\tleft_hits\tright_hits\tconcordant_placements\n')
        for r in report:
            rep.write('\t'.join(map(str, r)) + '\n')
    c = collections.Counter(r[2].split(':')[0] for r in report)
    print(f'rescue: {len(report)} copies -> ' + ', '.join(f'{k} {v}' for k, v in c.most_common())
          + f'; {written} new locus pairs', file=sys.stderr)


# ── account ────────────────────────────────────────────────────────────────────

def read_stat(path):
    out = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 2:
                    bank = f[0][2:] if f[0].startswith('./') else f[0]
                    out[bank] = f[1]
    return out


def cmd_account(args):
    stats = {}
    for path in args.stat:
        stats.update(read_stat(path))
    final = collections.defaultdict(list)       # cluster -> statuses of final alignments
    if os.path.exists(args.statpairs):
        with open(args.statpairs) as fh:
            for line in fh:
                aln = line.split('\t', 1)[0]
                final[aln.split('.', 1)[0]].append(aln.rsplit('.', 1)[1])
    rescue = collections.defaultdict(list)
    if args.rescue_report and os.path.exists(args.rescue_report):
        with open(args.rescue_report) as fh:
            fh.readline()
            for line in fh:
                cid, cl, outcome = line.split('\t')[:3]
                if outcome.startswith('resolved:'):
                    m = outcome.split(':', 1)[1]
                    st = final.get(m)
                    outcome = 'rescued:' + st[0] if st else 'rescued:qc_failed'
                rescue[cl].append(outcome.split(':')[0] if outcome.startswith('ambiguous') else outcome)

    totals = collections.Counter()
    with open(args.out, 'w') as out:
        out.write('cluster\tclass\tloci\tfate\trescue\n')
        for cls, path in (('double', args.double), ('multi', args.multi), ('poly', args.poly)):
            for cl, loci in read_clusters(path).items():
                if final.get(cl):
                    fate = 'resolved:' + ','.join(sorted(final[cl]))
                elif cls == 'double':
                    st = stats.get(f'{cl}.cl.dbl')
                    fate = f'rejected:{st}' if st else 'rejected'
                elif cls == 'multi':
                    fate = 'multi_unresolved'
                else:
                    fate = 'poly_unanalysed'
                rs = collections.Counter(rescue.get(cl, []))
                rtext = ','.join(f'{k}={v}' for k, v in sorted(rs.items())) or '.'
                out.write(f'{cl}\t{cls}\t{len(loci)}\t{fate}\t{rtext}\n')
                totals[(cls, fate.split(':')[0])] += 1
                for k, v in rs.items():
                    totals[('rescue', k.split(':')[0] + (':' + k.split(':')[1] if k.startswith('rescued') else ''))] += v
    for (cls, fate), n in sorted(totals.items()):
        print(f'  {cls:7s} {fate:20s} {n}', file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('prepare')
    p.add_argument('--clusters', nargs='+', required=True, help='<sp1>-<sp2>_multi and _poly cluster files')
    p.add_argument('--resolved', help='statpairs of the run (clusters that already have a final alignment)')
    p.add_argument('--copies', nargs=2, required=True, metavar='SP=BED', help='the *_uniq.bed copies of both species')
    p.add_argument('--genome', nargs=2, required=True, metavar='SP=FASTA', help='both genomes (with .fai)')
    p.add_argument('--out', required=True, help='output prefix')
    p.set_defaults(func=cmd_prepare)

    r = sub.add_parser('resolve')
    r.add_argument('--prep', required=True, help='prefix given to prepare')
    r.add_argument('--sam', nargs=2, required=True, metavar='SP=SAM',
                   help='bwa mem -a output of the flanks of species SP mapped to the other genome')
    r.add_argument('--genome', nargs=2, required=True, metavar='SP=FASTA')
    r.add_argument('--sine-length', type=int, required=True)
    r.add_argument('--anchors', type=int, default=10,
                   help='best hits of each flank used as anchors to search the other flank next to [10]')
    r.add_argument('--margin', type=int, default=20,
                   help='min score lead of the best concordant placement over the next one [20]')
    r.add_argument('--out', required=True, help='cluster lines for alignment (<sp1>-<sp2>_rescue)')
    r.add_argument('--report', required=True, help='per-copy outcome table')
    r.set_defaults(func=cmd_resolve)

    a = sub.add_parser('account')
    a.add_argument('--double'), a.add_argument('--multi'), a.add_argument('--poly')
    a.add_argument('--stat', nargs='+', default=[], help='stat_doubles, stat_multi, stat_rescue files')
    a.add_argument('--statpairs', required=True)
    a.add_argument('--rescue-report')
    a.add_argument('--out', required=True)
    a.set_defaults(func=cmd_account)

    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
