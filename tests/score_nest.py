#!/usr/bin/env python3
"""Score sine_nest.py scan output against simulation v2 truth.
Usage: score_nest.py SIMDIR PREFIX_TEMPLATE (e.g. nest_{sp}) [species ...]"""
import sys, csv, collections
S, tmpl, *sps = sys.argv[1:]
sps = sps or ['aaa', 'bbb', 'ccc', 'ddd']
nt = {r['insert']: r for r in csv.DictReader(open(f'{S}/nest_truth.tsv'), delimiter='\t')}
tot = collections.Counter()
for sp in sps:
    bed = [l.rstrip('\n').split('\t') for l in open(f'{S}/{sp}-SINEX.bed')]
    ins = {(b[0], int(b[1]), int(b[2])): b[3] for b in bed if b[3] in nt}
    sat = [(b[0], int(b[1]), int(b[2])) for b in bed if '_u' in b[3]]
    sat_ev = {b[3].rsplit('_u', 1)[0] for b in bed if '_u' in b[3]}
    pre = tmpl.format(sp=sp)
    found = set(); fp = 0; tsd_ok = 0; types_ok = 0; older = 0
    for r in csv.DictReader(open(f'{pre}.nest.tsv'), delimiter='\t'):
        s, e = map(int, r['insert'].split('-'))
        hit = [n for (c, a, b), n in ins.items() if c == r['chrom'] and min(b, e) - max(a, s) > 0.5 * (b - a)]
        if hit:
            found.add(hit[0]); tsd_ok += r['tsd_length'] != '0'
            types_ok += r['host_type'] == nt[hit[0]]['host_sub'] and r['insert_type'] == nt[hit[0]]['insert_sub']
            older += float(r['host_identity']) < float(r['insert_identity'])
        else:
            fp += 1
    sat_loci = [r for r in csv.DictReader(open(f'{pre}.loci.tsv'), delimiter='\t') if r['class'] == 'satellite']
    sat_found = sum(any(c == r['chrom'] and int(r['start']) <= a < int(r['end']) for c, a, b in sat) for r in sat_loci)
    print(f'{sp}: nested inserts found {len(found)}/{len(ins)}  false events {fp}  TSD found {tsd_ok}  '
          f'subfamilies right {types_ok}  host older by identity {older}/{len(found)} | '
          f'satellite loci {len(sat_loci)} (true {sat_found}) of {len(sat_ev)} arrays')
    tot['found'] += len(found); tot['true'] += len(ins); tot['fp'] += fp
print(f"total: {tot['found']}/{tot['true']} nested inserts, {tot['fp']} false events")
