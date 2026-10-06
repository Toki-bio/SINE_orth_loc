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
           PREFIX.edges.tsv     anchor graph: groups whose sites are neighbours along a genome,
                                with the spacer per species
           PREFIX.breakpoints.tsv  adjacencies of one genome broken in another; a_specific =
                                kept by no other genome (misjoin or lineage rearrangement)
           PREFIX.dupblocks.tsv runs of multicopy sites (duplicated or twice-assembled regions)
         States: P = SINE present, A = empty site (orthologous flanks, no SINE),
                 U = no data (no validated call for this species),
                 M = missing: the only evidence is a window clipped by a contig end (orth status MISSING),
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
import bisect
import collections
import glob
import gzip
import os
import re
import sys

LOCUS_RE = re.compile(r'^(.*):(\d+)-(\d+)\(([+-])\)$')
METRICS = ['LF', 'FL', 'LFcp', 'RFlength', 'FR', 'RFcp', 'OneTwo', 'OneSINE', 'TwoSINE', 'SL', 'SO', 'ST', 'FI']
ORTH_COLUMNS = ['alignment', 'cluster', 'status', 'species1', 'locus1', 'sine1',
                'species2', 'locus2', 'sine2'] + METRICS + ['anchor1', 'anchor2']
FLANK = 300  # left flank length used by SINE_orth_loc.bash


def die(msg):
    sys.exit(f'sine_registry.py: {msg}')


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def previous_groups(previous):
    """PREFIX.groups.tsv of a previous build, or its gzipped copy."""
    p = f'{previous}.groups.tsv'
    return p if os.path.exists(p) or not os.path.exists(p + '.gz') else p + '.gz'


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


def is_number(text):
    try:
        float(text)
        return True
    except ValueError:
        return False


def load_copies(specs, close):
    """Annotated SINE copies per species from BED files; several files per species allowed.

    sear2k output (<sp>-<FAMILY>.bed: name column = % identity, score = aligned query
    length, column 7 = bitscore) takes the family from the file name; other BED files
    give the family in the name column. Returns {species: [copy, ...]}; a copy is a dict
    with id, coordinates, family, identity, bitscore, junction (5' end of the copy in SINE
    orientation) and close (another copy of the same file within `close` bp, which
    SINE_orth_loc.bash excludes from the comparison)."""
    copies = collections.defaultdict(list)
    for spec in specs:
        sp, _, path = spec.partition('=')
        if not path:
            die(f'--copies expects SPECIES=copies.bed, got "{spec}"')
        base = re.sub(r'\.bed(\.gz)?$', '', os.path.basename(path))
        file_family = base.split('-', 1)[1] if '-' in base else base
        rows = []
        with open_text(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 6:
                    die(f'{path}: BED6 needed (chrom start end name score strand)')
                chrom, st, en, strand = f[0], int(f[1]), int(f[2]), f[5]
                sear2k = is_number(f[3])
                rows.append({'id': f'{sp}:{chrom}:{st}-{en}({strand})', 'species': sp, 'chrom': chrom,
                             'start': st, 'end': en, 'strand': strand,
                             'family': file_family if sear2k else f[3],
                             'identity': f[3] if sear2k and float(f[3]) <= 100 else '.',
                             'bitscore': f[6] if sear2k and len(f) > 6 else '.',
                             'subfamily': '.', 'subfamily_status': '.',
                             'junction': st if strand == '+' else en, 'close': False})
        rows.sort(key=lambda c: (c['chrom'], c['start']))
        # same rule as "bedtools cluster -d <close>" followed by keeping singletons
        for x, y in zip(rows, rows[1:]):
            if x['chrom'] == y['chrom'] and y['start'] - x['end'] <= close:
                x['close'] = y['close'] = True
        copies[sp] += rows
    return dict(copies)


def load_subfamilies(specs, copies):
    """Attach SINEderella step-2 assignments (assignment_full.tsv: Sequence = chrom:start-end(strand),
    Subfamily, Bitscore, Votes, Status, ...) to the copies of a species: same coordinates, or
    else the copy on the same strand that overlaps the sequence by at least half of the shorter."""
    for spec in specs:
        sp, _, path = spec.partition('=')
        if not path or sp not in copies:
            die(f'--subfamilies expects SPECIES=assignment_full.tsv for a species given with --copies, got "{spec}"')
        exact = {c['id'].split(':', 1)[1]: c for c in copies[sp]}
        by_chrom = collections.defaultdict(list)
        for c in copies[sp]:
            by_chrom[(c['chrom'], c['strand'])].append(c)
        matched = total = 0
        with open_text(path) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            for line in fh:
                r = dict(zip(header, line.rstrip('\n').split('\t')))
                seq = r.get('Sequence', '').split('|')[0]
                total += 1
                c = exact.get(seq)
                if c is None:
                    m = LOCUS_RE.match(seq)
                    if not m:
                        continue
                    chrom, st, en, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
                    best = 0
                    for cand in by_chrom.get((chrom, strand), ()):
                        ov = min(en, cand['end']) - max(st, cand['start'])
                        if ov > best and ov >= 0.5 * min(en - st, cand['end'] - cand['start']):
                            best, c = ov, cand
                if c is None:
                    continue
                matched += 1
                c['subfamily'] = r.get('Subfamily', '.') or '.'
                c['subfamily_status'] = r.get('Status', '.') or '.'
        print(f'{sp}: {matched}/{total} subfamily assignments matched to copies', file=sys.stderr)


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
    with open_text(previous_groups(previous)) as fh:
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
            near = []
            while i < len(lst) and lst[i][0] <= j + tol:
                near.append((abs(lst[i][0] - j), lst[i][1]))
                i += 1
            if near:          # only the nearest old site(s): neighbouring loci within --tol stay apart
                d = min(near)[0]
                for oid in {oid for dist, oid in near if dist == d}:
                    c[oid] += 1
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
    with open_text(previous_groups(previous)) as fh:
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
    precise = []          # node -> anchor taken from the alignment (exact) rather than the window
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
            for sp, locus, sine, anc in ((r['species1'], r['locus1'], r['sine1'], r.get('anchor1', '')),
                                         (r['species2'], r['locus2'], r['sine2'], r.get('anchor2', ''))):
                if sp not in species_seen:
                    species_seen.append(sp)
                chrom, s, e, strand = parse_locus(locus)
                n = uf.add()
                # insertion junction: from the alignment when the orth table has it (exact even for
                # merged/extended windows), else from the window and the most common window length
                a = int(anc) if anc not in ('', None) else insertion_anchor(chrom, s, e, strand, slop)
                nodes.append((sp, chrom, strand, a, s, e))
                precise.append(anc not in ('', None))
                calls.append(True if sine == '1' else False if sine == '0' else 'M' if sine == 'M' else None)
                ids.append(n)
            edges.append((ids[0], ids[1], os.path.basename(path), r))

    copies = load_copies(args.copies or [], args.close)
    load_subfamilies(args.subfamilies or [], copies)
    # copies flagged by sine_nest.py scan (role in column 5): satellite units are not independent
    # insertions (their groups become X); nested inserts and close/dimer copies keep their state
    caution = collections.defaultdict(list)          # (species, chrom) -> [(start, end, role)]
    for spec in args.caution or []:
        sp, path = spec.split('=', 1)
        with open_text(path) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) >= 5:
                    caution[(sp, f[0])].append((int(f[1]), int(f[2]), f[4]))
    def caution_role(sp, c):
        for s0, e0, role in caution.get((sp, c['chrom']), ()):
            if min(e0, c['end']) - max(s0, c['start']) > 0.5 * (c['end'] - c['start']):
                return role
        return None
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
        tol = args.precise_tol if precise[a] and precise[b] else args.tol
        if nodes[a][:3] == nodes[b][:3] and nodes[b][3] - nodes[a][3] <= tol:
            site.union(a, b)
            uf.union(a, b)
    for a, b, _, _ in edges:
        uf.union(a, b)

    # annotated copy at each site that carries the SINE. A site is usually the copy's 5' end,
    # but where two assemblies' annotations disagree on the copy's strand the flank lands at
    # its 3' end, so either end within --tol counts (5' end and same strand preferred).
    ends = collections.defaultdict(list)            # (species, chrom) -> sorted [(position, n, copy)]
    for sp, rows in copies.items():
        for n, c in enumerate(rows):
            ends[(sp, c['chrom'])].append((c['junction'], n, c))
            other_end = c['end'] if c['strand'] == '+' else c['start']
            ends[(sp, c['chrom'])].append((other_end, n, c))
    for v in ends.values():
        v.sort(key=lambda x: x[0])
    ends_pos = {k: [x[0] for x in v] for k, v in ends.items()}

    def copy_at(sp, chrom, strand, anchor, tol):
        v = ends.get((sp, chrom))
        if not v:
            return None
        pos = ends_pos[(sp, chrom)]
        best = None
        for i in range(bisect.bisect_left(pos, anchor - tol), bisect.bisect_right(pos, anchor + tol)):
            p, _, c = v[i]
            key = (p != c['junction'], c['strand'] != strand, abs(p - anchor))
            if best is None or key < best[0]:
                best = (key, c)
        return best[1] if best else None

    # sites of one species that land on the same copy are one locus
    node_copy = {}
    first_node_of_copy = {}
    for i, (sp, chrom, strand, anchor, _, _) in enumerate(nodes):
        if calls[i] is True and sp in copies:
            c = copy_at(sp, chrom, strand, anchor, args.precise_tol if precise[i] else args.tol)
            if c:
                node_copy[i] = c
                j = first_node_of_copy.setdefault(c['id'], i)
                if j != i:
                    site.union(i, j)
                    uf.union(i, j)

    comps = collections.defaultdict(list)
    for i in range(len(nodes)):
        comps[uf.find(i)].append(i)

    groups = []
    for members in comps.values():
        by_species = collections.defaultdict(lambda: collections.defaultdict(list))
        for i in members:
            by_species[nodes[i][0]][site.find(i)].append(i)
        states, cells, flags, sites, fams, group_copies = {}, {}, set(), [], collections.Counter(), []
        loci = []                                   # (species, chrom, position, P|A) of every site
        subfams = collections.Counter()
        for sp in species:
            sp_sites = by_species.get(sp)
            if not sp_sites:
                states[sp], cells[sp] = 'U', '.'
                continue
            # sites with a presence call; sites known only from a clipped window (M) or a link (None)
            # count only when the species has no called site
            called = {k: ids for k, ids in sp_sites.items() if any(calls[i] in (True, False) for i in ids)}
            if not called:
                ids = [i for v in sp_sites.values() for i in v]
                states[sp] = 'M' if any(calls[i] == 'M' for i in ids) else 'U'
                cells[sp] = ','.join(sorted({f'{nodes[i][1]}:{nodes[i][3]}({nodes[i][2]})' for i in ids})) if states[sp] == 'M' else '.'
                continue
            sp_sites = called
            seen = {calls[i] for ids in sp_sites.values() for i in ids if calls[i] in (True, False)}
            parts = []
            for ids in sp_sites.values():
                chrom, strand = nodes[ids[0]][1], nodes[ids[0]][2]
                anchor = sorted(nodes[i][3] for i in ids)[len(ids) // 2]
                has_sine = any(calls[i] is True for i in ids)
                found = collections.Counter(node_copy[i]['id'] for i in ids if i in node_copy)
                c = next(node_copy[i] for i in ids if i in node_copy and node_copy[i]['id'] == found.most_common(1)[0][0]) if found else None
                if c:
                    parts.append(f"{c['chrom']}:{c['start']}-{c['end']}({c['strand']})")
                    sites.append((sp, chrom, c['junction']))
                    loci.append((sp, chrom, c['junction'], 'P'))
                    fams[c['family']] += 1
                    if c['subfamily_status'] == 'assigned':
                        subfams[c['subfamily']] += 1
                    group_copies.append(c)
                else:
                    parts.append(f'{chrom}:{anchor}({strand})')
                    sites.append((sp, chrom, anchor))
                    loci.append((sp, chrom, anchor, 'P' if has_sine else 'A'))
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
        for c in group_copies:
            role = caution_role(c['species'], c) if caution else None
            if role == 'unit':
                flags.add(f"satellite:{c['species']}"); states[c['species']] = 'X'
            elif role in ('insert', 'close', 'dimer'):
                flags.add(f"{'nested' if role == 'insert' else role}:{c['species']}")
        if len(fams) > 1:
            flags.add('family_mixed:' + '/'.join(sorted(fams)))
        if len(subfams) > 1:
            flags.add('subfamily_mixed:' + '/'.join(sorted(subfams)))
        first = min(members, key=lambda i: (species.index(nodes[i][0]), nodes[i][1], nodes[i][3]))
        groups.append({'key': (species.index(nodes[first][0]), nodes[first][1], nodes[first][3]),
                       'members': members, 'states': states, 'cells': cells, 'flags': sorted(flags),
                       'sites': sites, 'loci': loci, 'family': fams.most_common(1)[0][0] if fams else '.',
                       'subfamily': subfams.most_common(1)[0][0] if subfams else '.',
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
        gt.write('\t'.join(['group', 'family', 'subfamily', 'pattern', 'n_P', 'n_A', 'n_U', 'n_X', 'n_M', 'pairs', 'flags'] + species) + '\n')
        mt.write('\t'.join(['group'] + species) + '\n')
        for g in groups:
            pat = ''.join(g['states'][sp] for sp in species)
            patterns[pat] += 1
            cnt = collections.Counter(pat)
            gt.write('\t'.join([g['id'], g['family'], g['subfamily'], pat, str(cnt['P']), str(cnt['A']), str(cnt['U']),
                                str(cnt['X']), str(cnt['M']), str(g['pairs']), ','.join(g['flags']) or '.']
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
        code = {'P': '1', 'A': '0', 'U': '?', 'X': '?', 'M': '?'}
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

    write_graph(pre, groups, species, args.max_spacer, args.min_dup_run)

    status = collections.Counter()
    if copies:
        copy_group = {}
        for g in groups:
            for c in g['copies']:
                copy_group[c['id']] = g['id']
        compared = set(species_seen)
        with open(f'{pre}.copies.tsv', 'w') as ct:
            ct.write('copy\tspecies\tchrom\tstart\tend\tstrand\tfamily\tidentity\tbitscore\t'
                     'subfamily\tsubfamily_status\tgroup\tstatus\n')
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
                                        c['family'], c['identity'], c['bitscore'], c['subfamily'],
                                        c['subfamily_status'], copy_group.get(c['id'], '.'), st]) + '\n')

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
          ' .edges.tsv .breakpoints.tsv .dupblocks.tsv'
          + (' .copies.tsv' if copies else ''), file=sys.stderr)
    print('top patterns (' + ''.join(s[0] for s in species) + ' = ' + ','.join(species) + '):', file=sys.stderr)
    for pat, n in patterns.most_common(args.top):
        print(f'  {pat}\t{n}', file=sys.stderr)


def write_graph(pre, groups, species, max_spacer, min_dup_run):
    """The anchor graph: groups are nodes, neighbouring sites along each genome are edges.

    PREFIX.edges.tsv        one row per pair of groups whose sites are neighbours in at least one
                            genome; per species the spacer (bp between the two sites) or '.'
    PREFIX.breakpoints.tsv  for each ordered species pair (a, b): neighbours in a, among groups with a
                            single site in both, that are not neighbours in b (other chromosome or
                            moved); a_specific = no other informative species keeps the adjacency
                            (a misjoin of assembly a, or a rearrangement on its lineage)
    PREFIX.dupblocks.tsv    runs of consecutive sites of multicopy groups in one genome: duplicated
                            regions, or haplotypes assembled twice
    """
    track = {}                                    # species -> chrom -> sorted [(pos, group, state)]
    count = {sp: collections.Counter() for sp in species}
    for g in groups:
        for sp, chrom, pos, state in g['loci']:
            track.setdefault(sp, collections.defaultdict(list))[chrom].append((pos, g['id'], state))
            count[sp][g['id']] += 1
    for chroms in track.values():
        for v in chroms.values():
            v.sort()

    edges = {}                                    # (g1, g2) -> {species: spacer}
    for sp, chroms in track.items():
        for v in chroms.values():
            for (p1, g1, _), (p2, g2, _) in zip(v, v[1:]):
                if g1 == g2 or (max_spacer and p2 - p1 > max_spacer):
                    continue
                key = (g1, g2) if g1 < g2 else (g2, g1)
                edges.setdefault(key, {})[sp] = p2 - p1
    with open(f'{pre}.edges.tsv', 'w') as fh:
        fh.write('\t'.join(['group1', 'group2', 'n_species'] + species) + '\n')
        for (g1, g2), sp_len in sorted(edges.items()):
            fh.write('\t'.join([g1, g2, str(len(sp_len))] + [str(sp_len.get(sp, '.')) for sp in species]) + '\n')

    single = {sp: {gid for gid, n in count[sp].items() if n == 1} for sp in species}
    loc = {sp: {gid: (chrom, pos) for chrom, v in track.get(sp, {}).items() for pos, gid, _ in v} for sp in species}

    def order(sp, keep):
        """per chromosome, the groups in keep in coordinate order; and group -> (chrom, rank)"""
        rank, runs = {}, []
        for chrom, v in track.get(sp, {}).items():
            run = [gid for _, gid, _ in v if gid in keep]
            runs.append((chrom, run))
            for k, gid in enumerate(run):
                rank[gid] = (chrom, k)
        return runs, rank

    calls = collections.defaultdict(dict)         # (a, g1, g2) -> {b: conserved?}
    rows = []
    summary = []
    for a in species:
        for b in species:
            if a == b:
                continue
            shared = single[a] & single[b]
            runs_a, _ = order(a, shared)
            _, rank_b = order(b, shared)
            n = kept = 0
            for chrom, run in runs_a:
                for g1, g2 in zip(run, run[1:]):
                    n += 1
                    (c1, k1), (c2, k2) = rank_b[g1], rank_b[g2]
                    ok = c1 == c2 and abs(k1 - k2) == 1
                    calls[(a, g1, g2)][b] = ok
                    if ok:
                        kept += 1
                    else:
                        rows.append([a, b, g1, g2, chrom, loc[a][g1][1], loc[a][g2][1],
                                     f'{c1}:{loc[b][g1][1]}', f'{c2}:{loc[b][g2][1]}',
                                     'other_chrom' if c1 != c2 else f'moved:{abs(loc[b][g2][1] - loc[b][g1][1])}'])
            summary.append((a, b, len(shared), n, kept))
    specific = collections.Counter()
    with open(f'{pre}.breakpoints.tsv', 'w') as fh:
        fh.write('species_a\tspecies_b\tgroup1\tgroup2\tchrom_a\tpos1_a\tpos2_a\tlocus1_b\tlocus2_b\tkind\ta_specific\n')
        for r in rows:
            c = calls[(r[0], r[2], r[3])]
            spec = len(c) >= 2 and not any(c.values())
            fh.write('\t'.join(map(str, r)) + '\t' + ('yes' if spec else 'no') + '\n')
    for (a, _, _), c in calls.items():
        if len(c) >= 2 and not any(c.values()):
            specific[a] += 1

    blocks = []
    for sp, chroms in track.items():
        for chrom, v in chroms.items():
            run = []
            for pos, gid, _ in v + [(None, None, None)]:
                if gid is not None and count[sp][gid] > 1:
                    run.append((pos, gid))
                    continue
                if len(run) >= min_dup_run:
                    blocks.append((sp, chrom, run[0][0], run[-1][0], len(run), ','.join(g for _, g in run)))
                run = []
    with open(f'{pre}.dupblocks.tsv', 'w') as fh:
        fh.write('species\tchrom\tstart\tend\tsites\tgroups\n')
        for blk in blocks:
            fh.write('\t'.join(map(str, blk)) + '\n')

    print(f'anchor graph: {len(edges)} edges; adjacency kept between genomes (groups with one site in both):',
          file=sys.stderr)
    for a, b, sh, n, kept in summary:
        if a < b:
            back = next(x for x in summary if x[0] == b and x[1] == a)
            print(f'  {a}-{b}: {sh} shared groups, {kept}/{n} adjacencies of {a} kept in {b}, '
                  f'{back[4]}/{back[3]} of {b} kept in {a}', file=sys.stderr)
    print('species-specific breakpoints (adjacency kept by no other species): '
          + ', '.join(f'{sp} {specific[sp]}' for sp in species), file=sys.stderr)
    dup = collections.Counter(blk[0] for blk in blocks)
    print(f'duplicated blocks (>= {min_dup_run} consecutive multicopy sites): '
          + ', '.join(f'{sp} {dup[sp]}' for sp in species), file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    b = sub.add_parser('build', help='combine orth tables into a pan-SINEome (orthologous groups)')
    b.add_argument('orth', nargs='+', help='orth_<a>-<b>.tsv tables, one per genome pair')
    b.add_argument('-o', '--output', required=True, help='output prefix')
    b.add_argument('--species', help='comma-separated species order for the outputs (default: sorted)')
    b.add_argument('--copies', action='append', metavar='SPECIES=BED',
                   help='annotated SINE copies of a species: sear2k <sp>-<FAMILY>.bed (family from the file '
                        'name) or BED6 with the family in the name column; repeat per species and family')
    b.add_argument('--caution', action='append', metavar='SPECIES=BED',
                   help='nest_<sp>.caution.bed from sine_nest.py scan: satellite units -> X, flags nested/close/dimer')
    b.add_argument('--subfamilies', action='append', metavar='SPECIES=TSV',
                   help='SINEderella step-2 assignment_full.tsv of a species; repeat per species')
    b.add_argument('--previous', metavar='PREFIX',
                   help='previous build (PREFIX.groups.tsv, PREFIX.aliases.tsv): keep its group IDs')
    b.add_argument('--id-prefix', default='PSG', help='prefix of new group IDs [PSG]')
    b.add_argument('--tol', type=int, default=60,
                   help='max distance (bp) between insertion sites of the same locus [60]')
    b.add_argument('--precise-tol', type=int, default=20,
                   help='max distance between insertion sites when both anchors come from alignments '
                        '(orth anchor1/anchor2 columns): keeps a nested insert apart from its host [20]')
    b.add_argument('--close', type=int, default=300,
                   help='copies closer than this are excluded by the pipeline (bedtools cluster -d) [300]')
    b.add_argument('--sine-length', type=int,
                   help='SINE consensus length used in the runs (default: inferred from locus window sizes)')
    b.add_argument('--top', type=int, default=15, help='patterns to print [15]')
    b.add_argument('--max-spacer', type=int, default=0,
                   help='do not link sites further apart than this (bp) in the edge table [0 = no limit]')
    b.add_argument('--min-dup-run', type=int, default=3,
                   help='consecutive multicopy sites that make a duplicated block [3]')
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
