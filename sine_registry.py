#!/usr/bin/env python3
"""Combine pairwise SINE_orth_loc results into a multi-species locus registry.

  build  orth_<a>-<b>.tsv ... -o PREFIX [--copies sp=copies.bed ...] [--previous OLD]
         Build the pan-SINEome: match the same genomic locus of a species across all
         pairwise comparisons, join loci through the orthologous pairs into groups, and write
           PREFIX.groups.tsv    one row per orthologous group: family, states, flags and,
                                per species, the copy (chrom:start-end(strand)) or the
                                insertion site without SINE (chrom:pos(strand))
           PREFIX.matrix.tsv    group x species states
           PREFIX.patterns.tsv  number of groups per state pattern
           PREFIX.nex           NEXUS binary matrix (P=1, A=0, U/X=?) of variable groups
           PREFIX.evidence.tsv  the pairwise rows each group was built from
           PREFIX.aliases.tsv   group IDs of the previous build that were merged, split or retired
           PREFIX.copies.tsv    with --copies: every annotated copy, its group or why it has none
         States: P = SINE present, A = empty site (orthologous flanks, no SINE),
                 U = no data (no validated call for this species),
                 X = ambiguous (contradicting calls, or several loci of the species
                     joined into one group).
         Group IDs are kept across rebuilds with --previous; see docs/pan-sineome.md.

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


def load_copies(specs, close):
    """Annotated SINE copies per species from BED files (name column = family).

    Returns {species: [copy, ...]}; a copy is a dict with id, coordinates, family,
    junction (5' end of the copy in SINE orientation) and close (another copy within
    `close` bp, which SINE_orth_loc.bash excludes from the comparison)."""
    copies = {}
    for spec in specs:
        sp, _, path = spec.partition('=')
        if not path:
            die(f'--copies expects SPECIES=copies.bed, got "{spec}"')
        rows = []
        with open_text(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 6:
                    die(f'{path}: BED6 needed (chrom start end name score strand)')
                chrom, s, e, strand = f[0], int(f[1]), int(f[2]), f[5]
                rows.append({'id': f'{sp}:{chrom}:{s}-{e}({strand})', 'species': sp, 'chrom': chrom,
                             'start': s, 'end': e, 'strand': strand, 'family': f[3],
                             'junction': s if strand == '+' else e, 'close': False})
        rows.sort(key=lambda c: (c['chrom'], c['start']))
        # same rule as "bedtools cluster -d <close>" followed by keeping singletons
        for a, b in zip(rows, rows[1:]):
            if a['chrom'] == b['chrom'] and b['start'] - a['end'] <= close:
                a['close'] = b['close'] = True
        copies[sp] = rows
    return copies


def parse_site(cell):
    """Sites of one species cell of a groups table: chrom:start-end(strand) for a copy,
    chrom:pos(strand) for an insertion site; returns [(chrom, junction)]."""
    out = []
    if cell in ('', '.'):
        return out
    for part in cell.split(','):
        m = re.match(r'^(.*):(\d+)-(\d+)\(([+-])\)$', part)
        if m:
            out.append((m.group(1), int(m.group(2)) if m.group(4) == '+' else int(m.group(3))))
            continue
        m = re.match(r'^(.*):(\d+)\(([+-])\)$', part)
        if m:
            out.append((m.group(1), int(m.group(2))))
    return out


def assign_ids(groups, species, previous, prefix, tol):
    """Stable group IDs. Without a previous build, number groups in coordinate order.
    With one, a group keeps the ID of the previous group it shares most sites with;
    merged IDs become aliases of the kept one, and split-off groups get new IDs."""
    aliases = []
    if not previous:
        return [f'{prefix}{k:07d}' for k in range(1, len(groups) + 1)], aliases

    old_sites = collections.defaultdict(list)   # (species, chrom) -> [(junction, old id)]
    max_num = 0
    with open(f'{previous}.groups.tsv') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            row = dict(zip(header, line.rstrip('\n').split('\t')))
            gid = row['group']
            m = re.search(r'(\d+)$', gid)
            if m:
                max_num = max(max_num, int(m.group(1)))
            for sp in header[header.index('flags') + 1:]:
                for chrom, j in parse_site(row.get(sp, '.')):
                    old_sites[(sp, chrom)].append((j, gid))
    try:
        with open(f'{previous}.aliases.tsv') as fh:
            fh.readline()
            for line in fh:
                aliases.append(line.rstrip('\n').split('\t'))
                m = re.search(r'(\d+)$', aliases[-1][0])
                if m:
                    max_num = max(max_num, int(m.group(1)))
    except FileNotFoundError:
        pass
    for v in old_sites.values():
        v.sort()

    import bisect
    shared = []                                   # per new group: Counter(old id -> shared sites)
    for g in groups:
        c = collections.Counter()
        for sp, chrom, j in g['sites']:
            lst = old_sites.get((sp, chrom), [])
            i = bisect.bisect_left(lst, (j - tol, ''))
            while i < len(lst) and lst[i][0] <= j + tol:
                c[lst[i][1]] += 1
                i += 1
        shared.append(c)

    winner = {}                                   # old id -> new group index
    for gi, c in enumerate(shared):
        for oid, n in c.items():
            if oid not in winner or n > shared[winner[oid]][oid]:
                winner[oid] = gi
    won = collections.defaultdict(list)
    for oid, gi in winner.items():
        won[gi].append(oid)

    ids = []
    for gi in range(len(groups)):
        if won[gi]:
            keep = min(won[gi])
            ids.append(keep)
            aliases += [[oid, keep, 'merged'] for oid in sorted(won[gi]) if oid != keep]
        else:
            max_num += 1
            ids.append(f'{prefix}{max_num:07d}')
            aliases += [[oid, ids[-1], 'split'] for oid in sorted(shared[gi])]
    seen_old = set(winner)
    with open(f'{previous}.groups.tsv') as fh:
        fh.readline()
        for line in fh:
            gid = line.split('\t', 1)[0]
            if gid not in seen_old:
                aliases.append([gid, '.', 'retired'])
    return ids, aliases


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
    calls = []            # node -> sine flag of the row side
    edges = []            # (node1, node2, source table, orth row)
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
            edges.append((ids[0], ids[1], os.path.basename(path), r))

    copies = load_copies(args.copies or [], args.close)
    species = args.species.split(',') if args.species else sorted(set(species_seen) | set(copies))
    unknown = (set(species_seen) | set(copies)) - set(species)
    if unknown:
        die(f'species in the inputs but not in --species: {",".join(sorted(unknown))}')

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
    for a, b, _, _ in edges:
        uf.union(a, b)

    # annotated copy at each site that carries the SINE: nearest 5' junction within --tol
    by_chrom = collections.defaultdict(list)
    for sp, rows in copies.items():
        for c in rows:
            by_chrom[(sp, c['chrom'])].append(c)

    def copy_at(sp, chrom, strand, anchor):
        best = None
        for c in by_chrom.get((sp, chrom), ()):
            d = abs(c['junction'] - anchor)
            if d <= args.tol and (best is None or (c['strand'] != strand, d) < best[0]):
                best = ((c['strand'] != strand, d), c)
        return best[1] if best else None

    comps = collections.defaultdict(list)
    for i in range(len(nodes)):
        comps[uf.find(i)].append(i)

    groups = []
    for members in comps.values():
        by_species = collections.defaultdict(lambda: collections.defaultdict(list))
        for i in members:
            by_species[nodes[i][0]][site.find(i)].append(i)
        states, cells, flags, sites, fams, group_copies = {}, {}, set(), [], collections.Counter(), []
        for sp in species:
            sp_sites = by_species.get(sp)
            if not sp_sites:
                states[sp], cells[sp] = 'U', '.'
                continue
            seen = {calls[i] for ids in sp_sites.values() for i in ids}
            parts = []
            for ids in sp_sites.values():
                chrom, strand = nodes[ids[0]][1], nodes[ids[0]][2]
                anchor = sorted(nodes[i][3] for i in ids)[len(ids) // 2]
                has_sine = any(calls[i] for i in ids)
                c = copy_at(sp, chrom, strand, anchor) if has_sine and copies.get(sp) is not None else None
                if c:
                    parts.append(f"{c['chrom']}:{c['start']}-{c['end']}({c['strand']})")
                    sites.append((sp, chrom, c['junction']))
                    fams[c['family']] += 1
                    group_copies.append(c)
                else:
                    parts.append(f'{chrom}:{anchor}({strand})')
                    sites.append((sp, chrom, anchor))
                    if has_sine and sp in copies:
                        flags.add(f'unannotated:{sp}')
            cells[sp] = ','.join(sorted(parts))
            if len(sp_sites) > 1:
                states[sp] = 'X'
                flags.add(f'multicopy:{sp}')
            elif len(seen) > 1:
                states[sp] = 'X'
                flags.add(f'inconsistent:{sp}')
            else:
                states[sp] = 'P' if True in seen else 'A'
        if len(fams) > 1:
            flags.add('family_mixed:' + '/'.join(sorted(fams)))
        first = min(members, key=lambda i: (species.index(nodes[i][0]), nodes[i][1], nodes[i][3]))
        groups.append({'key': (species.index(nodes[first][0]), nodes[first][1], nodes[first][3]),
                       'members': members, 'states': states, 'cells': cells, 'flags': sorted(flags),
                       'sites': sites, 'family': fams.most_common(1)[0][0] if fams else '.',
                       'copies': group_copies, 'pairs': sum(1 for i in members) // 2})
    groups.sort(key=lambda g: g['key'])
    ids, aliases = assign_ids(groups, species, args.previous, args.id_prefix, args.tol)

    pre = args.output
    group_of_node = {}
    for gid, g in zip(ids, groups):
        g['id'] = gid
        for i in g['members']:
            group_of_node[i] = gid

    patterns = collections.Counter()
    variable = []
    with open(f'{pre}.groups.tsv', 'w') as gt, open(f'{pre}.matrix.tsv', 'w') as mt:
        gt.write('\t'.join(['group', 'family', 'pattern', 'n_P', 'n_A', 'n_U', 'n_X', 'pairs', 'flags'] + species) + '\n')
        mt.write('\t'.join(['group'] + species) + '\n')
        for g in groups:
            pat = ''.join(g['states'][sp] for sp in species)
            patterns[pat] += 1
            cnt = collections.Counter(pat)
            gt.write('\t'.join([g['id'], g['family'], pat, str(cnt['P']), str(cnt['A']), str(cnt['U']),
                                str(cnt['X']), str(g['pairs']), ','.join(g['flags']) or '.']
                               + [g['cells'][sp] for sp in species]) + '\n')
            mt.write('\t'.join([g['id']] + [g['states'][sp] for sp in species]) + '\n')
            if cnt['P'] and cnt['A']:
                variable.append((g['id'], g['states']))

    with open(f'{pre}.patterns.tsv', 'w') as pt:
        pt.write('pattern\t' + '\t'.join(species) + '\tgroups\n')
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
        nx.write('[characters: ' + ' '.join(gid for gid, _ in variable) + ']\n')

    with open(f'{pre}.evidence.tsv', 'w') as et:
        et.write('\t'.join(['group', 'source'] + ORTH_COLUMNS) + '\n')
        for a, _, src, r in edges:
            et.write('\t'.join([group_of_node[a], src] + [r.get(k, '') for k in ORTH_COLUMNS]) + '\n')

    with open(f'{pre}.aliases.tsv', 'w') as at:
        at.write('alias\tgroup\tevent\n')
        for row in aliases:
            at.write('\t'.join(row) + '\n')

    status = collections.Counter()
    if copies:
        copy_group = {}
        for g in groups:
            for c in g['copies']:
                copy_group[c['id']] = g['id']
        compared = set(species_seen)
        with open(f'{pre}.copies.tsv', 'w') as ct:
            ct.write('copy\tspecies\tchrom\tstart\tend\tstrand\tfamily\tgroup\tstatus\n')
            for sp in species:
                for c in copies.get(sp, []):
                    if c['id'] in copy_group:
                        st = 'grouped'
                    elif sp not in compared:
                        st = 'not_compared'
                    elif c['close']:
                        st = 'close_copy'
                    else:
                        st = 'no_validated_pair'
                    status[st] += 1
                    ct.write('\t'.join([c['id'], sp, c['chrom'], str(c['start']), str(c['end']), c['strand'],
                                        c['family'], copy_group.get(c['id'], '.'), st]) + '\n')

    flagged = sum(1 for g in groups if g['flags'])
    print(f'{len(edges)} locus pairs from {len(tables)} tables'
          + (f' ({same_species} same-species pairs skipped)' if same_species else ''), file=sys.stderr)
    print(f'{len(groups)} orthologous groups, {flagged} flagged, {len(variable)} variable (P and A)', file=sys.stderr)
    if copies:
        print('copies: ' + ', '.join(f'{k} {v}' for k, v in status.most_common()), file=sys.stderr)
    if aliases:
        ev = collections.Counter(a[2] for a in aliases)
        print('IDs vs previous build: ' + ', '.join(f'{k} {v}' for k, v in sorted(ev.items())), file=sys.stderr)
    print(f'-> {pre}.groups.tsv .matrix.tsv .patterns.tsv .nex .evidence.tsv .aliases.tsv'
          + (' .copies.tsv' if copies else ''), file=sys.stderr)
    print('top patterns (' + ''.join(s[0] for s in species) + ' = ' + ','.join(species) + '):', file=sys.stderr)
    for pat, n in patterns.most_common(args.top):
        print(f'  {pat}\t{n}', file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    b = sub.add_parser('build', help='combine orth tables into a pan-SINEome (orthologous groups)')
    b.add_argument('orth', nargs='+', help='orth_<a>-<b>.tsv tables, one per genome pair')
    b.add_argument('-o', '--output', required=True, help='output prefix')
    b.add_argument('--species', help='comma-separated species order for the outputs (default: sorted)')
    b.add_argument('--copies', action='append', metavar='SPECIES=BED',
                   help='all annotated SINE copies of a species (BED6, name column = family); repeat per species')
    b.add_argument('--previous', metavar='PREFIX',
                   help='previous build (PREFIX.groups.tsv, PREFIX.aliases.tsv): keep its group IDs')
    b.add_argument('--id-prefix', default='PSG', help='prefix of new group IDs [PSG]')
    b.add_argument('--tol', type=int, default=60,
                   help='max distance (bp) between insertion sites of the same locus [60]')
    b.add_argument('--close', type=int, default=300,
                   help='copies closer than this are excluded by the pipeline (bedtools cluster -d) [300]')
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
