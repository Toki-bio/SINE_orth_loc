#!/usr/bin/env python3
"""Real-data labelled set for ComPair.sh calls, from read genotypes of the assembled individual.

Labels come from an independent source: reads of individual 245 (the animal behind the dvl assembly)
genotyped at the groups of a pan-SINEome built WITHOUT dvl (sine_genotype.py). For every row of an orth
table of a pair involving dvl, the other species' locus identifies the group, and the reads say whether
dvl carries the insertion there. The row's dvl-side call (sine flag) is then correct or not.

Selection (to keep label noise out):
  * read genotype homozygous (0/0 or 1/1) with GQ >= --min-gq (default 20); heterozygous loci are
    polymorphic in the individual and excluded;
  * the group must be unflagged for multicopy/inconsistent (several loci joined: the label is ambiguous);
  * the other species' locus must match exactly one group within --tol bp.
Writes a TSV with pair, alignment, status, dvl call, read genotype, label (1 = call correct), the
ComPair.sh metrics, and the dvl scaffold (for spatially held-out validation).

Usage: build_training_set.py --orth ORTH.tsv[.gz] ... --groups nodvl.groups.tsv.gz --genotypes v245.gt.tsv -o OUT.tsv
"""
import argparse, bisect, collections, csv, gzip, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from sine_registry import read_orth, insertion_anchor, parse_locus, FLANK  # noqa: E402

METRICS = ['LF', 'FL', 'LFcp', 'RFlength', 'FR', 'RFcp', 'OneTwo', 'OneSINE', 'TwoSINE', 'SL', 'SO', 'ST']
CELL = re.compile(r'^(.*):(\d+)(?:-(\d+))?\(([+-])\)$')


def opn(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--orth', nargs='+', required=True)
    ap.add_argument('--groups', required=True)
    ap.add_argument('--genotypes', required=True)
    ap.add_argument('--target', default='dvl')
    ap.add_argument('--min-gq', type=int, default=20)
    ap.add_argument('--tol', type=int, default=60)
    ap.add_argument('-o', '--output', required=True)
    a = ap.parse_args()

    gt = {r['group']: r for r in csv.DictReader(open(a.genotypes), delimiter='\t')}
    idx = collections.defaultdict(list)          # (species, chrom) -> [(pos, group)]
    gflags = {}
    with opn(a.groups) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        species = header[header.index('flags') + 1:]
        for line in fh:
            r = dict(zip(header, line.rstrip('\n').split('\t')))
            gflags[r['group']] = r['flags']
            for sp in species:
                for cell in r[sp].split(','):
                    m = CELL.match(cell)
                    if not m:
                        continue
                    ch, s, e, st = m.group(1), int(m.group(2)), m.group(3), m.group(4)
                    pts = [s] if e is None else [s, int(e)]     # a copy matches at either end
                    for p in pts:
                        idx[(sp, ch)].append((p, r['group']))
    for v in idx.values():
        v.sort()

    def group_at(sp, ch, pos):
        v = idx.get((sp, ch), [])
        i = bisect.bisect_left(v, (pos - a.tol, ''))
        hits = set()
        while i < len(v) and v[i][0] <= pos + a.tol:
            hits.add(v[i][1]); i += 1
        return hits

    stats = collections.Counter()
    with open(a.output, 'w') as out:
        out.write('\t'.join(['pair', 'alignment', 'status', 'dvl_call', 'reads_GT', 'reads_GQ', 'label',
                             'dvl_scaffold'] + METRICS) + '\n')
        for path in a.orth:
            pair = re.sub(r'^orth_|\.tsv(\.gz)?$', '', os.path.basename(path))
            rows = read_orth(path)
            lengths = collections.Counter()
            for r in rows:
                for k in ('locus1', 'locus2'):
                    _, s, e, _ = parse_locus(r[k]); lengths[e - s] += 1
            slop = lengths.most_common(1)[0][0] - FLANK if lengths else 0
            for r in rows:
                if a.target not in (r['species1'], r['species2']) or r['species1'] == r['species2']:
                    continue
                t, o = ('1', '2') if r['species1'] == a.target else ('2', '1')
                ch, s, e, st = parse_locus(r['locus' + o])
                pos = insertion_anchor(ch, s, e, st, slop)
                groups = group_at(r['species' + o], ch, pos)
                stats['rows'] += 1
                if len(groups) != 1:
                    stats['no unique group'] += 1; continue
                g = groups.pop()
                if any(f.startswith(('multicopy', 'inconsistent')) for f in gflags.get(g, '').split(',')):
                    stats['group flagged multicopy/inconsistent'] += 1; continue
                rg = gt.get(g)
                if not rg or rg['GT'] not in ('0/0', '1/1') or int(rg['GQ']) < a.min_gq:
                    stats['no confident homozygous read genotype'] += 1; continue
                call = r['sine' + t]
                label = int((call == '1') == (rg['GT'] == '1/1'))
                stats[f"labelled {r['status']} {'correct' if label else 'INCORRECT'}"] += 1
                dch = parse_locus(r['locus' + t])[0]
                out.write('\t'.join([pair, r['alignment'], r['status'], call, rg['GT'], rg['GQ'], str(label), dch]
                                    + [r.get(m, '') for m in METRICS]) + '\n')
    for k, v in sorted(stats.items()):
        print(f'{k}\t{v}', file=sys.stderr)


if __name__ == '__main__':
    main()
