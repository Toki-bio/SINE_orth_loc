#!/usr/bin/env python3
"""Compare SINE_orth_loc runs of one genome pair made with different flank settings, on real data.

Usage: compare_real.py --pair sp1 sp2 --chunks DIR --beds DIR RUN_DIR [RUN_DIR ...]
  --chunks  directory with <a>_on_<b>.junction.bed: the insertion junction of every input copy, lifted to the other
            genome through an alignment of the copy +-10 kb (independent of the flank settings tested)
  --beds    directory with <sp>-*.bed, the input copies (column 6 strand)
Only the main flank-mapping rows (clusters C<n>R, status SINE/PM/MP) are compared; rescue (M) and nested (N) rows
are left aside. A row is checked on the genome that carries the SINE (both genomes for SINE rows): the copy whose
junction is the row's anchor (within 40 bp) must lift to within 50 bp of the row's anchor in the other genome.
Per run: rows, agree / disagree / untested; then rows shared by runs and rows found by one run only, with the
same check, and rows where two runs pair one locus with different partners."""
import argparse, bisect, collections, csv, glob, os, re
ap = argparse.ArgumentParser()
ap.add_argument('--pair', nargs=2, required=True)
ap.add_argument('--chunks', required=True)
ap.add_argument('--beds', required=True)
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
    f = glob.glob(os.path.join(run, f'orth_{A}-{B}.tsv'))[0]
    out = []
    for r in csv.DictReader(open(f), delimiter='\t'):
        if not r['cluster'].startswith('C') or r['status'] not in ('SINE', 'PM', 'MP'): continue
        if not r['anchor1'] or not r['anchor2']: continue
        side = {}
        for k in ('1', '2'):
            c = r['locus' + k].rsplit(':', 1)[0]
            side[r['species' + k]] = (c, int(r['anchor' + k]), r['sine' + k])
        if A in side and B in side: out.append((r['status'], side[A], side[B]))
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
