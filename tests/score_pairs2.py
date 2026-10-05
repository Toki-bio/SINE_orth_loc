#!/usr/bin/env python3
"""Score orth_<a>-<b>.tsv against simulation v2 truth, per element class.
Usage: score_pairs2.py SIMDIR RUNDIR sp1 sp2
A row is correct if both loci resolve to the same element (copy for sine=1, empty site for sine=0)
and the presence states agree with truth; 'paralog' = a haplotig/duplicate copy was used."""
import sys, re, collections, csv
S, run, a, b = sys.argv[1:5]
truth = {r['event']: r for r in csv.DictReader(open(f'{S}/truth.tsv'), delimiter='\t')}
def base(n): return re.sub(r'(_h[12]|_hap|_u\d+)$', '', n)
copies = collections.defaultdict(list); sites = collections.defaultdict(list)
for sp in (a, b):
    for l in open(f'{S}/{sp}-SINEX.bed'):
        ch, s, e, n, _, st = l.rstrip().split('\t'); s, e = int(s), int(e)
        copies[(sp, ch)].append((s if st == '+' else e, n))
    for l in open(f'{S}/{sp}.sites.tsv'):
        n, ch, cur, st = l.rstrip().split('\t')
        if n in truth and truth[n][sp] == 'A':
            sites[(sp, ch)].append((int(cur), n))
rows = list(csv.DictReader(open(f'{run}/orth_{a}-{b}.tsv'), delimiter='\t'))
L = collections.Counter(int(m.group(3)) - int(m.group(2)) for r in rows for k in ('locus1', 'locus2')
                        for m in [re.match(r'(.*):(\d+)-(\d+)\((.)\)', r[k])])
slop = L.most_common(1)[0][0] - 300 if L else 0
def ev(sp, locus, sine):
    m = re.match(r'(.*):(\d+)-(\d+)\((.)\)', locus); ch, s, e, st = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
    anc = e - slop if st == '+' else s + slop
    pool = copies[(sp, ch)] if sine == '1' else sites[(sp, ch)]
    best = min(pool, key=lambda x: abs(x[0] - anc), default=None)
    return best[1] if best and abs(best[0] - anc) <= 80 else None
res = collections.defaultdict(collections.Counter); found = collections.defaultdict(set)
ends = set(l.strip() for l in open(f'{S}/contig_ends.tsv')) if __import__('os').path.exists(f'{S}/contig_ends.tsv') else set()
def ev_any(sp, locus):
    return ev(sp, locus, '1') or ev(sp, locus, '0')
missing = collections.Counter()
for r in rows:
    if r['status'] == 'MISSING':
        side = '1' if r['sine1'] == '?' else '2'
        e = ev_any(r['species' + side], r['locus' + side])
        missing['at simulated contig end' if e and base(e) in ends else 'elsewhere'] += 1
        continue
    e1 = ev(r['species1'], r['locus1'], r['sine1']); e2 = ev(r['species2'], r['locus2'], r['sine2'])
    if e1 is None or e2 is None:
        res['?']['unmatched'] += 1; continue
    b1, b2 = base(e1), base(e2); cls = truth.get(b1, {}).get('class', '?')
    if b1 != b2: res[cls]['wrong'] += 1
    elif '_hap' in e1 + e2: res[cls]['paralog'] += 1
    elif cls == 'satellite': res[cls]['satellite_call'] += 1
    else: res[cls]['correct'] += 1; found[cls].add(b1)
det = collections.defaultdict(set)
for n, t in truth.items():
    if t['class'] in ('otherTE', 'satellite'): continue
    if 'P' in (t[a], t[b]): det[t['class']].add(n)
print(f'{a}-{b}: rows={len(rows)} unmatched={res["?"]["unmatched"]}' + (f' MISSING rows: {dict(missing)}' if missing else ''))
for cls in sorted(set(det) | set(res) - {'?'}):
    c = res[cls]
    print(f'  {cls:10s} recall {len(found[cls]):3d}/{len(det[cls]):3d}  correct={c["correct"]} wrong={c["wrong"]} paralog={c["paralog"]}'
          + (f' satellite_calls={c["satellite_call"]}' if c['satellite_call'] else ''))
