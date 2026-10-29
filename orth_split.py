#!/usr/bin/env python3
"""Split the alignment bundles of a SINE_orth_loc pair by reference (bank mode: subfamily, or span of one), so that
one subfamily's SINE, PM or MP alignments open on their own in sine_loci_browser.html.

Usage: orth_split.py orth_<sp1>-<sp2>.tsv [--bundles aln_<sp1>-<sp2>_SINE.aln.gz ...] [-o DIR] [--min-loci 20]

Writes DIR/<reference>/aln_<sp1>-<sp2>_<status>.aln.gz (status = SINE, PM, MP) and DIR/index.tsv (reference, dir,
status, alignments). References with fewer than --min-loci rows go together under 'other' (same cut as orth_summary's
table on the Tal page); reference names are made directory-safe. Default bundles: the three next to the orth table;
default DIR: by_reference next to it. Alignments whose name is not in the orth table (none expected) are counted
and reported, not written."""
import argparse
import collections
import csv
import gzip
import os
import re


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('orth')
    ap.add_argument('--bundles', nargs='+')
    ap.add_argument('-o', '--out')
    ap.add_argument('--min-loci', type=int, default=20)
    a = ap.parse_args()
    d = os.path.dirname(os.path.abspath(a.orth))
    pair = re.sub(r'^orth_|\.tsv$', '', os.path.basename(a.orth))
    bundles = a.bundles or [os.path.join(d, f'aln_{pair}_{s}.aln.gz') for s in ('SINE', 'PM', 'MP')]
    out = a.out or os.path.join(d, 'by_reference')

    ref_of, n_ref = {}, collections.Counter()
    for r in csv.DictReader(open(a.orth), delimiter='\t'):
        if r['status'] in ('SINE', 'PM', 'MP') and r['alignment']:
            ref_of[r['alignment']] = r.get('reference') or '(no reference)'
            n_ref[ref_of[r['alignment']]] += 1
    dirs = {ref: re.sub(r'[^A-Za-z0-9._-]+', '_', ref) if n >= a.min_loci else 'other' for ref, n in n_ref.items()}
    os.makedirs(out, exist_ok=True)

    handles, counts, unlisted = {}, collections.Counter(), collections.Counter()
    for b in bundles:
        name, lines = None, []

        def flush():
            if name is None:
                return
            status = name.rsplit('.', 1)[-1]
            ref = ref_of.get(name)
            if ref is None:
                unlisted[status] += 1
                return
            key = (ref, status)
            if key not in handles:
                os.makedirs(os.path.join(out, dirs[ref]), exist_ok=True)
                handles[key] = gzip.open(os.path.join(out, dirs[ref], f'aln_{pair}_{status}.aln.gz'), 'at')
            handles[key].write(f'##FILE {name}\n' + ''.join(lines))
            counts[key] += 1

        with gzip.open(b, 'rt') as fh:
            for line in fh:
                if line.startswith('##FILE '):
                    flush()
                    name, lines = line[7:].strip(), []
                elif name is not None:
                    lines.append(line)
            flush()
    for h in handles.values():
        h.close()
    with open(os.path.join(out, 'index.tsv'), 'w') as fh:
        fh.write('reference\tdir\tstatus\talignments\n')
        for (ref, status), n in sorted(counts.items(), key=lambda kv: (-n_ref[kv[0][0]], kv[0])):
            fh.write(f'{ref}\t{dirs[ref]}\t{status}\t{n}\n')
    print(f'{sum(counts.values())} alignments, {len(set(dirs.values()))} directories -> {out}'
          + (f'; not in the orth table: {dict(unlisted)}' if unlisted else ''))


if __name__ == '__main__':
    main()
