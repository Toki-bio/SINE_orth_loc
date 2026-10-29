#!/usr/bin/env python3
"""Summary of a SINE_orth_loc pair per reference (bank mode: subfamily, or span of one), written as
summary_<sp1>-<sp2>.tsv next to the orth table.

Usage: orth_summary.py orth_<sp1>-<sp2>.tsv --copies SP1=BED SP2=BED [-o OUT]

Per reference and in total: input copies per species (BED labels), loci = rows of the orth table (C flank-mapped,
M rescued, N nested), status counts (SINE = in both, PM = in sp1 only, MP = in sp2 only), shared % = SINE / loci,
identity of the copies to the reference (median and 10-90 % range over every copy-bearing locus; divergence =
100 - identity), identity between the two copies of SINE loci, identity of the left flanks (orthology evidence).
Identities are ComPair.sh's: OneSINE/TwoSINE (locus vs reference), OneTwo (locus vs locus), FL (left flanks);
nested rows carry none. Rows made before the reference column (older runs) are summed under '(no reference)'."""
import argparse
import collections
import csv
import os
import statistics


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))] if xs else None


def fmt(xs):
    if not xs:
        return ''
    return f'{statistics.median(xs):.0f} ({q(xs, 0.1):.0f}-{q(xs, 0.9):.0f})'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('orth')
    ap.add_argument('--copies', nargs=2, required=True, metavar='SP=BED')
    ap.add_argument('-o', '--out')
    a = ap.parse_args()
    sps, copies = [], collections.defaultdict(collections.Counter)
    for spec in a.copies:
        sp, bed = spec.split('=', 1)
        sps.append(sp)
        for line in open(bed):
            f = line.rstrip('\n').split('\t')
            if len(f) >= 4:
                copies[sp][f[3]] += 1
    rows = collections.defaultdict(list)
    for r in csv.DictReader(open(a.orth), delimiter='\t'):
        if r['status'] not in ('SINE', 'PM', 'MP'):
            continue
        rows[r.get('reference') or '(no reference)'].append(r)

    def summ(ref, rs):
        st = collections.Counter(r['status'] for r in rs)
        ty = collections.Counter(r['cluster'][0] for r in rs)
        to_ref, pair, flank = [], [], []
        for r in rs:
            for k, col in (('1', 'OneSINE'), ('2', 'TwoSINE')):
                if r['sine' + k] == '1' and r[col]:
                    to_ref.append(float(r[col]))
            if r['status'] == 'SINE' and r['OneTwo']:
                pair.append(float(r['OneTwo']))
            if r['FL']:
                flank.append(float(r['FL']))
        n = len(rs)
        return [ref, copies[sps[0]][ref], copies[sps[1]][ref], n, ty['C'], ty['M'], ty['N'],
                st['SINE'], st['PM'], st['MP'], f'{100 * st["SINE"] / n:.1f}' if n else '',
                fmt(to_ref), fmt(pair), fmt(flank)]

    out = a.out or os.path.join(os.path.dirname(os.path.abspath(a.orth)),
                                os.path.basename(a.orth).replace('orth_', 'summary_'))
    with open(out, 'w') as fh:
        fh.write('\t'.join(['reference', f'copies_{sps[0]}', f'copies_{sps[1]}', 'loci', 'C', 'M', 'N',
                            'SINE', f'PM_{sps[0]}_only', f'MP_{sps[1]}_only', 'shared_pct',
                            'identity_to_reference', 'identity_between_copies', 'identity_left_flanks']) + '\n')
        for ref in sorted(rows, key=lambda x: -len(rows[x])):
            fh.write('\t'.join(map(str, summ(ref, rows[ref]))) + '\n')
        allrows = [r for rs in rows.values() for r in rs]
        tot = summ('total', allrows)
        tot[1], tot[2] = sum(copies[sps[0]].values()), sum(copies[sps[1]].values())
        fh.write('\t'.join(map(str, tot)) + '\n')
    print(f'{len(allrows)} loci, {len(rows)} references -> {out}')


if __name__ == '__main__':
    main()
