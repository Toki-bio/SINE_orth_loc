#!/usr/bin/env python3
"""Score read-based genotypes against the simulation truth.

score_genotypes.py [--min-gq=N] SIMDIR GROUPS.tsv GENOTYPES.tsv HAP_SPECIES...  (e.g. aaa ccc for a hybrid)
A group's event is found through its annotated copies (BED names); the true genotype is the
number of the sample's haplotypes whose species carries the event.
"""
import collections
import csv
import sys

args = sys.argv[1:]
min_gq = 0
if args[0].startswith('--min-gq='):
    min_gq = int(args.pop(0).split('=')[1])
sim, groups_path, gt_path, *haps = args
truth = {r['event']: r for r in csv.DictReader(open(f'{sim}/truth.tsv'), delimiter='\t')}
names = {}
for sp in ('aaa', 'bbb', 'ccc', 'ddd'):
    for line in open(f'{sim}/{sp}-SINEX.bed'):
        f = line.rstrip('\n').split('\t')
        names[f'{f[0]}:{f[1]}-{f[2]}({f[5]})'] = f[3].replace('_dup', '')
event = {}
with open(groups_path) as fh:
    header = fh.readline().rstrip('\n').split('\t')
    species = header[header.index('flags') + 1:]
    for line in fh:
        g = dict(zip(header, line.rstrip('\n').split('\t')))
        evs = {names[c] for sp in species for c in g[sp].split(',') if c in names}
        if len(evs) == 1:
            event[g['group']] = evs.pop()
res = collections.Counter()
conf = collections.Counter()
rep = collections.Counter()
for r in csv.DictReader(open(gt_path), delimiter='\t'):
    e = event.get(r['group'])
    if e is None:
        res['no event'] += 1
        continue
    true = sum(truth[e][h] == 'P' for h in haps)
    tgt = {0: '0/0', 1: '0/1', 2: '1/1'}[true]
    if r['GT'] == './.' or int(r['GQ']) < min_gq:
        res['no call'] += 1
    else:
        ok = r['GT'] == tgt
        res['correct' if ok else 'wrong'] += 1
        conf[(tgt, r['GT'])] += 1
        if e.endswith('_rep'):
            rep['correct' if ok else 'wrong'] += 1
called = res['correct'] + res['wrong']
print(f"called {called}/{called + res['no call']}  correct {res['correct']}  wrong {res['wrong']}"
      f"  ({100 * res['correct'] / max(1, called):.1f}%)  | repeat-flank loci: correct {rep['correct']} wrong {rep['wrong']}"
      + ('  | errors: ' + ', '.join(f'{t}->{c} x{n}' for (t, c), n in sorted(conf.items()) if t != c) if res['wrong'] else ''))
