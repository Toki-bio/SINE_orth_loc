#!/usr/bin/env python3
"""Compare SINE_orth_loc runs of one genome pair made with different flank settings, on real data.

Usage: compare_real.py --pair sp1 sp2 --chunks DIR --beds DIR RUN_DIR [RUN_DIR ...]
  --chunks  directory with <a>_on_<b>.junction.bed: the insertion junction of every input copy, lifted to the other
            genome through an alignment of the copy +-10 kb (independent of the flank settings tested)
  --beds    directory with <sp>-*.bed, the input copies (column 6 strand)
Only the main flank-mapping rows (clusters C<n>R, status SINE/PM/MP) are compared; rescue (M) and nested (N) rows
are left aside. A row is checked on the genome that carries the SINE (both genomes for SINE rows): the copy whose
junction is the row's junction (within 40 bp; junction = window start + FLANK, see rows_of) must lift to within 50 bp
of the row's junction in the other genome.
Per run: rows, agree / disagree / untested; then rows shared by runs and rows found by one run only, with the
same check, and rows where two runs pair one locus with different partners.
--out DIR: also write DIR/loci.tsv (one line per locus that differs between runs: the category, the check, the
reference row = subfamily or span the locus was judged against, and the alignment name of that locus in every run) and DIR/<category>/<run>/aln_<sp1>-<sp2>_<type>.aln.gz, small bundles with
only those alignments, to open in sine_loci_browser.html. Categories: disagree (some run's row disagrees with the lift),
partner (runs pair the locus with different loci), status (runs give different statuses), only_<run> (one run only).
For a run without a row at the locus, its alignments whose window holds the locus are taken from all its bundles
(rejected ones too) and named in loci.tsv with a leading "~"."""
import argparse, bisect, collections, csv, glob, gzip, os, re, sys
ap = argparse.ArgumentParser()
ap.add_argument('--pair', nargs=2, required=True)
ap.add_argument('--chunks', required=True)
ap.add_argument('--beds', required=True)
ap.add_argument('--out')
ap.add_argument('runs', nargs='+')
a = ap.parse_args()
A, B = a.pair
LOC = re.compile(r'^(.*):(\d+)-(\d+)\(([+-])\)$')
TOL_COPY, TOL_LIFT, TOL_MATCH = 40, 50, 40

# input copies: per genome, sorted junctions -> copy id
junc = {}
for sp in (A, B):
    d = collections.defaultdict(list)
    for f in glob.glob(os.path.join(a.beds, f'{sp}-*.bed')):
        for l in open(f):
            c, s, e, *_r = l.rstrip('\n').split('\t'); st = _r[2]
            d[c].append((int(s) if st == '+' else int(e), f'{c}:{s}-{e}({st})'))
    for c in d: d[c].sort()
    junc[sp] = d
def copy_at(sp, chrom, x):
    L = junc[sp].get(chrom, [])
    i = bisect.bisect_left(L, (x - TOL_COPY, ''))
    best = None
    while i < len(L) and L[i][0] <= x + TOL_COPY:
        if best is None or abs(L[i][0] - x) < abs(best[0] - x): best = L[i]
        i += 1
    return best[1] if best else None
lift = {}
for q, t in ((A, B), (B, A)):
    m = {}
    for l in open(os.path.join(a.chunks, f'{q}_on_{t}.junction.bed')):
        f = l.rstrip('\n').split('\t')
        cid = re.sub(r'_\d+_\d+$', '', f[3])            # paftools liftover appends _<start>_<end>
        m.setdefault(cid, []).append((f[0], int(f[1])))
    lift[q] = m

def rows_of(run):
    """rows as (status, A side, B side, alignment); a side = (chrom, junction, sine). The junction is the anchor
    column, unless it lies > 50 bp from the window start + FLANK (+) or end - FLANK (-) (FLANK = the run's most
    common anchor offset): about 5 % of anchors lie inside the flank by an amount that changes with FLANK, and then
    the window position is used."""
    f = glob.glob(os.path.join(run, f'orth_{A}-{B}.tsv'))[0]
    rows = [r for r in csv.DictReader(open(f), delimiter='\t')
            if r['cluster'].startswith('C') and r['status'] in ('SINE', 'PM', 'MP') and r['anchor1'] and r['anchor2']]
    # reference row of each cluster (subfamily or span): orth column, or refs_<A>-<B>.tsv of runs made before it
    refs = {}
    rf = os.path.join(run, f'refs_{A}-{B}.tsv')
    if os.path.exists(rf):
        refs = {r['cluster']: r['reference'] for r in csv.DictReader(open(rf), delimiter='\t')}
    off = collections.Counter()
    for r in rows:
        m = LOC.match(r['locus1'])
        off[int(r['anchor1']) - int(m.group(2)) if m.group(4) == '+' else int(m.group(3)) - int(r['anchor1'])] += 1
    F = off.most_common(1)[0][0]
    print(f'{os.path.basename(os.path.normpath(run))}: FLANK {F}', file=sys.stderr)
    out = []
    for r in rows:
        side = {}
        for k in ('1', '2'):
            c, s, e, st = LOC.match(r['locus' + k]).groups()
            w = int(s) + F if st == '+' else int(e) - F; x = int(r['anchor' + k])
            side[r['species' + k]] = (c, x if abs(x - w) <= 50 else w, r['sine' + k])
        if A in side and B in side:
            out.append((r['status'], side[A], side[B], r['alignment'], r.get('reference') or refs.get(r['cluster'], '')))
    return out

def check(sa, sb):
    """agree / disagree / untested for one row (sa, sb: chrom, anchor, sine of A and B)"""
    res = []
    for sp, me, other in ((A, sa, sb), (B, sb, sa)):
        if me[2] != '1': continue
        cid = copy_at(sp, me[0], me[1])
        hits = lift[sp].get(cid) if cid else None
        if not hits: res.append('untested'); continue
        ok = any(h[0] == other[0] and abs(h[1] - other[1]) <= TOL_LIFT for h in hits)
        res.append('agree' if ok else 'disagree')
    if 'disagree' in res: return 'disagree'
    return 'agree' if 'agree' in res else 'untested'

runs = {os.path.basename(os.path.normpath(r)): rows_of(r) for r in a.runs}
print('run\trows\tSINE\tPM\tMP\tagree\tdisagree\tuntested')
for n, rows in runs.items():
    st = collections.Counter(r[0] for r in rows); ch = collections.Counter(check(r[1], r[2]) for r in rows)
    print(n, len(rows), st['SINE'], st['PM'], st['MP'], ch['agree'], ch['disagree'], ch['untested'], sep='\t')

# rows keyed by the A-side anchor; match across runs within TOL_MATCH on A, then compare the B side
def key(sa): return (sa[0], sa[1] // 100)
index = {}
for n, rows in runs.items():
    d = collections.defaultdict(list)
    for r in rows: d[key(r[1])].append(r)
    index[n] = d
def find(n, sa):
    for k in ((sa[0], sa[1] // 100 + o) for o in (-1, 0, 1)):
        for r in index[n].get(k, []):
            if abs(r[1][1] - sa[1]) <= TOL_MATCH: return r
    return None
names = list(runs)
print('\nrun\tonly_this_run\tagree\tdisagree\tuntested\tother_partner_in_some_run\tother_status_in_some_run')
for n in names:
    only = collections.Counter(); partner = status = 0
    for r in runs[n]:
        found = [find(m, r[1]) for m in names if m != n]
        if not any(found): only[check(r[1], r[2])] += 1; only['n'] += 1
        for f in found:
            if not f: continue
            if f[2][0] != r[2][0] or abs(f[2][1] - r[2][1]) > TOL_MATCH: partner += 1; break
        for f in found:
            if f and f[0] != r[0]: status += 1; break
    print(n, only['n'], only['agree'], only['disagree'], only['untested'], partner, status, sep='\t')

if a.out:
    # loci across runs: rows matched on the A-side anchor (within TOL_MATCH), first row of each run kept
    loci = []; lindex = collections.defaultdict(list)
    for n in names:
        for r in runs[n]:
            hit = None
            for k in ((r[1][0], r[1][1] // 100 + o) for o in (-1, 0, 1)):
                for L in lindex.get(k, []):
                    if abs(L['a'][1] - r[1][1]) <= TOL_MATCH: hit = L; break
                if hit: break
            if not hit:
                hit = {'a': r[1], 'rows': {}}; loci.append(hit); lindex[key(r[1])].append(hit)
            hit['rows'].setdefault(n, r)
    def cats(L):
        rs = L['rows']; c = []
        if any(check(r[1], r[2]) == 'disagree' for r in rs.values()): c.append('disagree')
        bs = [r[2] for r in rs.values()]
        if any(b[0] != bs[0][0] or abs(b[1] - bs[0][1]) > TOL_MATCH for b in bs): c.append('partner')
        if len({r[0] for r in rs.values()}) > 1: c.append('status')
        if len(rs) == 1: c.append('only_' + next(iter(rs)))
        return c
    # every alignment header of every run: window on the A genome -> (alignment name, bundle type)
    HDR = re.compile(r'^>(.*):(\d+)-(\d+)\([+-]\)')
    windows = {}
    for run in a.runs:
        n = os.path.basename(os.path.normpath(run)); d = collections.defaultdict(list)
        for b in glob.glob(os.path.join(run, f'aln_{A}-{B}_*.aln.gz')):
            t = b.rsplit('_', 1)[1].split('.')[0]; name = None
            for l in gzip.open(b, 'rt'):
                if l.startswith('##FILE '): name = l[7:].rstrip('\n')
                elif l.startswith('>') and name:
                    m = HDR.match(l)
                    if m and m.group(1) in junc[A]: d[m.group(1)].append((int(m.group(2)), int(m.group(3)), name, t))
        for c in d: d[c].sort()
        windows[n] = d
    def at(n, sa):
        L = windows[n].get(sa[0], []); i = bisect.bisect_left(L, (sa[1] - 5000,)); out = []
        while i < len(L) and L[i][0] <= sa[1]:
            if L[i][1] >= sa[1] and (L[i][2], L[i][3]) not in out: out.append((L[i][2], L[i][3]))
            i += 1
        return out
    want = collections.defaultdict(set)        # (category, run) -> alignment names
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, 'loci.tsv'), 'w') as fo:
        fo.write('\t'.join(['locus', 'category', A + '_anchor', 'check', 'reference'] + names) + '\n')
        k = 0
        for L in loci:
            c = cats(L)
            if not c: continue
            k += 1; cells = []
            for n in names:
                r = L['rows'].get(n)
                if r:
                    al = [r[3]]; cells.append(f"{r[0]} {r[2][0]}:{r[2][1]} {r[3]} ({check(r[1], r[2])})")
                else:
                    al = [x[0] for x in at(n, L['a'])]; cells.append('~' + ','.join(al) if al else '-')
                for x in c:
                    for nm in al: want[(x, n)].add(nm)
            chk = ','.join(sorted({check(r[1], r[2]) for r in L['rows'].values()}))
            refc = ','.join(sorted({r[4] for r in L['rows'].values() if r[4]}))
            fo.write('\t'.join([f'L{k}', ','.join(c), f"{L['a'][0]}:{L['a'][1]}", chk, refc] + cells) + '\n')
    for (x, n), nms in sorted(want.items()):
        run = next(r for r in a.runs if os.path.basename(os.path.normpath(r)) == n)
        d = os.path.join(a.out, x, n); os.makedirs(d, exist_ok=True); got = 0
        for b in glob.glob(os.path.join(run, f'aln_{A}-{B}_*.aln.gz')):
            o = gzip.open(os.path.join(d, os.path.basename(b)), 'wt'); o.write('##SINE_orth_loc bundle v1\n'); p = False
            for l in gzip.open(b, 'rt'):
                if l.startswith('##FILE '): p = l[7:].rstrip('\n') in nms; got += p
                if p: o.write(l)
            o.close()
        print(f'{x}\t{n}\t{got} alignments', file=sys.stderr)
    print(f'\n{k} loci written to {a.out}/loci.tsv')
