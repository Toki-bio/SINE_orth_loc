#!/usr/bin/env python3
"""Combine pairwise SINE_orth_loc results into a multi-species locus registry.

  build  orth_<a>-<b>.tsv ... -o PREFIX
         Match the same genomic locus of a species across all pairwise comparisons,
         join loci through the orthologous pairs, and write
           PREFIX.loci.tsv      one row per multi-species locus: states, flags, coordinates
           PREFIX.matrix.tsv    locus x species states
           PREFIX.patterns.tsv  number of loci per state pattern
           PREFIX.nex           NEXUS binary matrix (P=1, A=0, U/X=?) of variable loci
         States: P = SINE present, A = empty site (orthologous flanks, no SINE),
                 U = no data (no validated call for this species),
                 X = ambiguous (contradicting calls, or several loci of the species
                     joined into one multi-species locus).

  orth   --stat stat_doubles_<a>-<b> stat_multi_<a>-<b> --aln DIR|BUNDLE ... \\
         --species a=a.fai b=b.fai -o orth_<a>-<b>.tsv
         Rebuild the orth table for runs made before SINE_orth_loc.bash wrote it,
         from the ComPair.sh stat files and the final alignments (loose files in a
         directory, or aln_*.aln[.gz] bundles).

Only the Python standard library is used.
"""

import argparse
import collections
import glob
import gzip
import os
import re
import sys

LOCUS_RE = re.compile(r'^(.*):(\d+)-(\d+)\(([+-])\)$')
METRICS = ['LF', 'FL', 'LFcp', 'RFlength', 'FR', 'RFcp', 'OneTwo', 'OneSINE', 'TwoSINE', 'SL', 'SO', 'ST']
ORTH_COLUMNS = ['alignment', 'cluster', 'status', 'species1', 'locus1', 'sine1',
                'species2', 'locus2', 'sine2'] + METRICS
FLANK = 300  # left flank length used by SINE_orth_loc.bash


def die(msg):
    sys.exit(f'sine_registry.py: {msg}')


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def parse_locus(text):
    m = LOCUS_RE.match(text)
    if not m:
        die(f'cannot parse locus "{text}" (expected chrom:start-end(strand))')
    return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


# ── orth: rebuild orth tables for older runs ──────────────────────────────────

def alignment_headers(sources):
    """Map alignment file name -> (seq1 locus, seq2 locus) from directories and bundles."""
    heads = {}

    def first_two(lines):
        hs = []
        for line in lines:
            if line.startswith('>'):
                hs.append(line[1:].split('::')[0].strip())
                if len(hs) == 2:
                    break
        return hs

    for src in sources:
        if os.path.isdir(src):
            for path in glob.glob(os.path.join(src, '*')):
                if re.search(r'\.(PM|MP|SINE)$', path):
                    with open(path) as fh:
                        heads[os.path.basename(path)] = first_two(fh)
        else:
            name, hs = None, []
            with open_text(src) as fh:
                for line in fh:
                    if line.startswith('##FILE '):
                        if name is not None:
                            heads[name] = hs
                        name, hs = line[7:].strip(), []
                    elif line.startswith('>') and name is not None and len(hs) < 2:
                        hs.append(line[1:].split('::')[0].strip())
            if name is not None:
                heads[name] = hs
    return heads


def cmd_orth(args):
    species = []
    chrom_species = {}
    for spec in args.species:
        name, _, fai = spec.partition('=')
        species.append(name)
        if fai:
            with open(fai) as fh:
                for line in fh:
                    chrom_species[line.split('\t', 1)[0]] = name
    if len(species) != 2:
        die('orth needs exactly two --species (NAME or NAME=genome.fai)')

    def species_of(locus):
        chrom = parse_locus(locus)[0]
        if chrom in chrom_species:
            return chrom_species[chrom]
        hits = [s for s in species if chrom.startswith(s)]
        if len(hits) != 1:
            die(f'cannot tell the species of chromosome "{chrom}": give --species NAME=genome.fai')
        return hits[0]

    heads = alignment_headers(args.aln)
    rows, missing = [], 0
    for stat in args.stat:
        with open(stat) as fh:
            for line in fh:
                f = line.split()
                if len(f) < 2 or f[1] not in ('PM', 'MP', 'SINE'):
                    continue
                bank = f[0][2:] if f[0].startswith('./') else f[0]
                aln = f'{bank}.{f[1]}'
                hs = heads.get(aln)
                if not hs or len(hs) < 2:
                    missing += 1  # alignment removed later by the pipeline (e.g. multi length filter)
                    continue
                metrics = dict(kv.split('=', 1) for kv in f[2:] if '=' in kv)
                st = f[1]
                rows.append([aln, aln.split('.')[0], st,
                             species_of(hs[0]), hs[0], '1' if st in ('PM', 'SINE') else '0',
                             species_of(hs[1]), hs[1], '1' if st in ('MP', 'SINE') else '0']
                            + [metrics.get(k, '') for k in METRICS])
    with open(args.output, 'w') as out:
        out.write('\t'.join(ORTH_COLUMNS) + '\n')
        for r in rows:
            out.write('\t'.join(r) + '\n')
    print(f'{len(rows)} locus pairs -> {args.output}'
          + (f' ({missing} stat lines without a final alignment skipped)' if missing else ''),
          file=sys.stderr)


# ── build: multi-species registry ─────────────────────────────────────────────

class UnionFind:
    def __init__(self):
        self.parent = []

    def add(self):
        self.parent.append(len(self.parent))
        return len(self.parent) - 1

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def read_orth(path):
    with open_text(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if header[:9] != ORTH_COLUMNS[:9]:
            die(f'{path} is not an orth table (run "sine_registry.py orth" for older runs)')
        return [dict(zip(header, line.rstrip('\n').split('\t'))) for line in fh if line.strip()]


def insertion_anchor(chrom, start, end, strand, slop):
    """Position of the flank/SINE junction on the 5' side of the SINE.

    Locus windows are the mapped left flank extended 3' by slop (SINE length + 300),
    so the junction sits slop bp inside the 3' end of the window."""
    return end - slop if strand == '+' else start + slop


def cmd_build(args):
    tables = []
    for path in args.orth:
        rows = read_orth(path)
        lengths = collections.Counter()
        for r in rows:
            for k in ('locus1', 'locus2'):
                _, s, e, _ = parse_locus(r[k])
                lengths[e - s] += 1
        if args.sine_length:
            slop = args.sine_length + FLANK
        elif lengths:
            slop = lengths.most_common(1)[0][0] - FLANK  # full-length windows: flank + slop
        else:
            slop = 0
        tables.append((path, rows, slop))

    uf = UnionFind()
    nodes = []            # node -> (species, chrom, strand, anchor, start, end)
    calls = []            # node -> sine flag of the row side (1/0)
    edges = []
    same_species = 0
    species_seen = []
    for path, rows, slop in tables:
        for r in rows:
            if r['species1'] == r['species2']:
                same_species += 1  # paralogous pair within one genome: not an orthology edge
                continue
            ids = []
            for sp, locus, sine in ((r['species1'], r['locus1'], r['sine1']),
                                    (r['species2'], r['locus2'], r['sine2'])):
                if sp not in species_seen:
                    species_seen.append(sp)
                chrom, s, e, strand = parse_locus(locus)
                n = uf.add()
                nodes.append((sp, chrom, strand, insertion_anchor(chrom, s, e, strand, slop), s, e))
                calls.append(sine == '1')
                ids.append(n)
            edges.append(ids)

    species = args.species.split(',') if args.species else sorted(species_seen)
    unknown = set(species_seen) - set(species)
    if unknown:
        die(f'species in the tables but not in --species: {",".join(sorted(unknown))}')

    # same genomic locus seen in several comparisons: same species, chromosome and
    # strand, insertion anchors within --tol bp (single linkage along the chromosome)
    site = UnionFind()
    for _ in nodes:
        site.add()
    order = sorted(range(len(nodes)), key=lambda i: nodes[i][:4])
    for a, b in zip(order, order[1:]):
        if nodes[a][:3] == nodes[b][:3] and nodes[b][3] - nodes[a][3] <= args.tol:
            site.union(a, b)
            uf.union(a, b)
    for a, b in edges:
        uf.union(a, b)

    comps = collections.defaultdict(list)
    for i in range(len(nodes)):
        comps[uf.find(i)].append(i)

    loci = []
    for members in comps.values():
        by_species = collections.defaultdict(lambda: collections.defaultdict(list))
        for i in members:
            by_species[nodes[i][0]][site.find(i)].append(i)
        states, coords, flags = {}, {}, set()
        for sp in species:
            sites = by_species.get(sp)
            if not sites:
                states[sp], coords[sp] = 'U', '.'
                continue
            spans = []
            for ids in sites.values():
                chrom, strand = nodes[ids[0]][1], nodes[ids[0]][2]
                spans.append(f'{chrom}:{min(nodes[i][4] for i in ids)}-{max(nodes[i][5] for i in ids)}({strand})')
            coords[sp] = ','.join(sorted(spans))
            seen = {calls[i] for i in members if nodes[i][0] == sp}
            if len(sites) > 1:
                states[sp] = 'X'
                flags.add(f'multicopy:{sp}')
            elif len(seen) > 1:
                states[sp] = 'X'
                flags.add(f'inconsistent:{sp}')
            else:
                states[sp] = 'P' if True in seen else 'A'
        first = min(members, key=lambda i: (species.index(nodes[i][0]), nodes[i][1], nodes[i][4]))
        loci.append((species.index(nodes[first][0]), nodes[first][1], nodes[first][4],
                     states, coords, sorted(flags), len(members) // 2))
    loci.sort(key=lambda x: x[:3])

    pre = args.output
    patterns = collections.Counter()
    variable = []
    with open(f'{pre}.loci.tsv', 'w') as lt, open(f'{pre}.matrix.tsv', 'w') as mt:
        lt.write('\t'.join(['locus', 'pattern', 'n_P', 'n_A', 'n_U', 'n_X', 'pairs', 'flags'] + species) + '\n')
        mt.write('\t'.join(['locus'] + species) + '\n')
        for k, (_, _, _, states, coords, flags, npairs) in enumerate(loci, 1):
            lid = f'L{k:06d}'
            pat = ''.join(states[sp] for sp in species)
            patterns[pat] += 1
            cnt = collections.Counter(pat)
            lt.write('\t'.join([lid, pat, str(cnt['P']), str(cnt['A']), str(cnt['U']), str(cnt['X']),
                                str(npairs), ','.join(flags) or '.'] + [coords[sp] for sp in species]) + '\n')
            mt.write('\t'.join([lid] + [states[sp] for sp in species]) + '\n')
            if cnt['P'] and cnt['A']:
                variable.append((lid, states))

    with open(f'{pre}.patterns.tsv', 'w') as pt:
        pt.write('pattern\t' + '\t'.join(species) + '\tloci\n')
        for pat, n in sorted(patterns.items(), key=lambda x: (-x[1], x[0])):
            pt.write(f'{pat}\t' + '\t'.join(pat) + f'\t{n}\n')

    with open(f'{pre}.nex', 'w') as nx:
        nx.write('#NEXUS\n[SINE presence/absence from sine_registry.py: 1 = SINE, 0 = empty site, ? = no data or ambiguous]\n')
        nx.write(f'BEGIN DATA;\n  DIMENSIONS NTAX={len(species)} NCHAR={len(variable)};\n'
                 '  FORMAT DATATYPE=STANDARD SYMBOLS="01" MISSING=?;\n  MATRIX\n')
        code = {'P': '1', 'A': '0', 'U': '?', 'X': '?'}
        for sp in species:
            nx.write(f'    {sp}  ' + ''.join(code[st[sp]] for _, st in variable) + '\n')
        nx.write('  ;\nEND;\n')
        nx.write('[characters: ' + ' '.join(lid for lid, _ in variable) + ']\n')

    total = len(loci)
    flagged = sum(1 for l in loci if l[5])
    print(f'{sum(len(r) for _, r, _ in tables)} locus pairs from {len(tables)} tables'
          + (f' ({same_species} same-species pairs skipped)' if same_species else ''), file=sys.stderr)
    print(f'{total} multi-species loci, {flagged} flagged, {len(variable)} variable (P and A) -> '
          f'{pre}.loci.tsv .matrix.tsv .patterns.tsv .nex', file=sys.stderr)
    print('top patterns (' + ''.join(s[0] for s in species) + ' = ' + ','.join(species) + '):', file=sys.stderr)
    for pat, n in patterns.most_common(args.top):
        print(f'  {pat}\t{n}', file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    b = sub.add_parser('build', help='combine orth tables into a multi-species registry')
    b.add_argument('orth', nargs='+', help='orth_<a>-<b>.tsv tables, one per genome pair')
    b.add_argument('-o', '--output', required=True, help='output prefix')
    b.add_argument('--species', help='comma-separated species order for the outputs (default: sorted)')
    b.add_argument('--tol', type=int, default=60,
                   help='max distance (bp) between insertion sites of the same locus in different comparisons [60]')
    b.add_argument('--sine-length', type=int,
                   help='SINE consensus length used in the runs (default: inferred from locus window sizes)')
    b.add_argument('--top', type=int, default=15, help='patterns to print [15]')
    b.set_defaults(func=cmd_build)

    o = sub.add_parser('orth', help='rebuild an orth table from stat files and alignments of an older run')
    o.add_argument('--stat', nargs='+', required=True, help='stat_doubles_<a>-<b> and stat_multi_<a>-<b>')
    o.add_argument('--aln', nargs='+', required=True, help='directories with .PM/.MP/.SINE files, or aln_*.aln[.gz] bundles')
    o.add_argument('--species', nargs=2, required=True, metavar='NAME[=FAI]',
                   help='the two species; with =genome.fai chromosomes are assigned by the index, '
                        'otherwise by chromosome name prefix')
    o.add_argument('-o', '--output', required=True)
    o.set_defaults(func=cmd_orth)

    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
