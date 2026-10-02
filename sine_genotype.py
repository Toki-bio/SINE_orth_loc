#!/usr/bin/env python3
"""Genotype pan-SINEome loci (orthologous SINE insertions) directly from sequencing reads.

  templates  For every orthologous group of a pan-SINEome (sine_registry.py build), write two
             template sequences: the allele with the SINE (left flank + copy + right flank,
             from a genome where the copy is annotated) and the empty allele (flank + flank,
             from a genome with an empty site, or the SINE allele with the copy cut out).
             Writes PREFIX.fa and PREFIX.tsv (junction positions in each template).

  Map the reads of a sample to PREFIX.fa, e.g.
      bwa index PREFIX.fa; bwa mem -t 8 PREFIX.fa reads.fq > sample.sam          (short reads)
      minimap2 -ax map-hifi -t 8 PREFIX.fa reads.fq > sample.sam                  (HiFi)

  call       (= count + genotype for one SAM; use count on each SAM and genotype on all count
             files to combine several read files of one sample.)
             Decide for every confidently placed read (MAPQ >= --min-mapq, primary or
             supplementary) which allele it shows at a group: a junction crossed with at least
             --anchor bp on both sides, or the SINE missing from the SINE template (deletion) /
             a SINE-sized insertion in the empty template. One vote per read and group. Junctions
             where most reads are ambiguously placed (repetitive flank) are not used. The
             genotype (0/0 absent, 0/1 heterozygous, 1/1 present; --ploidy 1: 0 or 1) maximises
             a binomial likelihood that accounts for the SINE allele offering two junctions to
             short reads.

Only the Python standard library is used.
"""

import argparse
import collections
import math
import re
import sys

LOCUS = re.compile(r'^(.*):(\d+)-(\d+)\(([+-])\)$')
SITE = re.compile(r'^(.*):(\d+)\(([+-])\)$')
CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


def die(msg):
    sys.exit(f'sine_genotype.py: {msg}')


def open_text(path):
    import gzip
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


class Fasta:
    """Random access to a FASTA file through its samtools .fai index."""

    def __init__(self, path):
        self.fh = open(path, 'rb')
        self.idx = {}
        with open(path + '.fai') as fai:
            for line in fai:
                name, length, offset, lb, lw = line.split('\t')[:5]
                self.idx[name] = (int(length), int(offset), int(lb), int(lw))

    def fetch(self, chrom, start, end):
        length, offset, lb, lw = self.idx[chrom]
        start, end = max(0, start), min(length, end)
        if end <= start:
            return ''
        a = offset + (start // lb) * lw + start % lb
        b = offset + ((end - 1) // lb) * lw + (end - 1) % lb + 1
        self.fh.seek(a)
        return self.fh.read(b - a).decode().replace('\n', '').replace('\r', '').upper()


# ── templates ─────────────────────────────────────────────────────────────────

def cmd_templates(args):
    genomes = {}
    for spec in args.genome:
        sp, _, path = spec.partition('=')
        genomes[sp] = Fasta(path)
    F = args.flank
    made = collections.Counter()
    with open_text(args.groups) as fh, open(f'{args.out}.fa', 'w') as fa, open(f'{args.out}.tsv', 'w') as tv:
        header = fh.readline().rstrip('\n').split('\t')
        species = header[header.index('flags') + 1:]
        tv.write('group\tallele\tsource\tlength\tjunctions\n')
        for line in fh:
            g = dict(zip(header, line.rstrip('\n').split('\t')))
            pattern = g['pattern']
            p_src = a_src = None
            for sp, st in zip(species, pattern):
                cell = g[sp]
                if sp not in genomes or ',' in cell:
                    continue
                if st == 'P' and p_src is None and LOCUS.match(cell):
                    p_src = (sp, LOCUS.match(cell))
                elif st == 'A' and a_src is None and SITE.match(cell):
                    a_src = (sp, SITE.match(cell))
            if p_src is None:
                made['skipped: no genome with an annotated copy'] += 1
                continue
            sp, m = p_src
            chrom, cs, ce = m.group(1), int(m.group(2)), int(m.group(3))
            seq = genomes[sp].fetch(chrom, cs - F, ce + F)
            left = min(F, cs)
            fa.write(f'>{g["group"]}|P\n{seq}\n')
            tv.write(f'{g["group"]}\tP\t{sp}:{chrom}:{cs}-{ce}\t{len(seq)}\t{left},{left + ce - cs}\n')
            if a_src:
                sp, m = a_src
                chrom, pos = m.group(1), int(m.group(2))
                seq = genomes[sp].fetch(chrom, pos - F, pos + F)
                tv.write(f'{g["group"]}\tA\t{sp}:{chrom}:{pos}\t{len(seq)}\t{min(F, pos)}\n')
                made['empty allele from a genome'] += 1
            else:                                   # no genome without the SINE: cut the copy out
                seq = genomes[sp].fetch(chrom, cs - F, cs) + genomes[sp].fetch(chrom, ce, ce + F)
                tv.write(f'{g["group"]}\tA\tsynthetic:{sp}:{chrom}:{cs}-{ce}\t{len(seq)}\t{left}\n')
                made['empty allele synthetic (copy cut out)'] += 1
            fa.write(f'>{g["group"]}|A\n{seq}\n')
    for k, v in made.most_common():
        print(f'{v}\t{k}', file=sys.stderr)


# ── call ──────────────────────────────────────────────────────────────────────

def parse_alignment(cigar, start, max_del):
    """Reference spans (split at deletions > max_del), large deletions and insertions."""
    ref = start
    spans, dels, ins = [], [], []
    cur = None
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op in 'M=X':
            if cur is None:
                cur = [ref, ref + n]
            else:
                cur[1] = ref + n
            ref += n
        elif op in 'DN':
            if n > max_del and cur is not None:
                spans.append(cur)
                dels.append((ref, ref + n))
                cur = None
            elif cur is not None:
                cur[1] = ref + n
            ref += n
        elif op == 'I':
            ins.append((ref, n))
    if cur is not None:
        spans.append(cur)
    return spans, dels, ins


def load_templates(path):
    tpl, groups = {}, []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            g, allele, _, _, js = line.rstrip('\n').split('\t')
            tpl[f'{g}|{allele}'] = [int(x) for x in js.split(',')]
            if allele == 'P':
                groups.append(g)
    return tpl, groups


COUNT_COLS = ['P_reads', 'A_reads', 'hq_P0', 'hq_P1', 'hq_A', 'lq_P0', 'lq_P1', 'lq_A']


def count_sam(args, tpl):
    """Per group: reads showing each allele, and confident / ambiguous reads per junction."""
    W, T = args.anchor, args.tol
    read_ev = collections.defaultdict(set)
    hq, lq = collections.Counter(), collections.Counter()
    n_len = sum_len = 0
    with open_text(args.sam) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            f = line.split('\t', 11)
            flag = int(f[1])
            if flag & (4 | 256):
                continue
            if not flag & 2048 and f[9] != '*' and n_len < 100000:
                n_len += 1
                sum_len += len(f[9])
            tname = f[2]
            js = tpl.get(tname)
            if js is None:
                continue
            g, allele = tname.rsplit('|', 1)
            spans, dels, ins = parse_alignment(f[5], int(f[3]) - 1, args.max_del)
            cover = lambda a, b: any(s <= a and e >= b for s, e in spans)
            found = []
            if allele == 'P':
                j1, j2 = js
                for k, j in enumerate(js):
                    if cover(j - W, j + W):
                        found.append((f'P{k}', 'P'))
                for ds, de in dels:
                    if abs(ds - j1) <= T and abs(de - j2) <= T and cover(j1 - W, ds) and cover(de, j2 + W):
                        found.append(('A', 'A'))
            else:
                j = js[0]
                if cover(j - W, j + W):
                    found.append(('A', 'A'))
                for pos, n in ins:
                    if abs(pos - j) <= T and n >= args.min_ins and cover(j - W, pos) and cover(pos, j + W):
                        found.append(('P0', 'P'))
            confident = int(f[4]) >= args.min_mapq
            for jid, al in found:
                if confident:
                    hq[(g, jid)] += 1
                    read_ev[(f[0], g)].add(al)
                else:
                    lq[(g, jid)] += 1
    counts = collections.defaultdict(lambda: [0] * len(COUNT_COLS))
    for (read, g), al in read_ev.items():
        if len(al) == 1:
            counts[g][0 if 'P' in al else 1] += 1
    for (g, jid), n in hq.items():
        counts[g][{'P0': 2, 'P1': 3, 'A': 4}[jid]] += n
    for (g, jid), n in lq.items():
        counts[g][{'P0': 5, 'P1': 6, 'A': 7}[jid]] += n
    return counts, n_len, sum_len


def genotype(args, tpl, groups, counts, mean_len, out_path):
    W = args.anchor
    e = args.error
    fs = [0.0, 0.5, 1.0] if args.ploidy == 2 else [0.0, 1.0]
    labels = ['0/0', '0/1', '1/1'] if args.ploidy == 2 else ['0', '1']
    nocall = './.' if args.ploidy == 2 else '.'
    tally = collections.Counter()
    w = max(1.0, mean_len - 2 * W)
    with open(out_path, 'w') as out:
        out.write('group\tP_reads\tA_reads\tP_junctions_informative\tGT\tGQ\n')
        for g in groups:
            k, a, hq0, hq1, hqa, lq0, lq1, lqa = counts.get(g, [0] * len(COUNT_COLS))
            j1, j2 = tpl[f'{g}|P']
            # a junction is informative unless most reads at it are ambiguously placed
            m = (lq0 <= hq0) + (lq1 <= hq1)
            a_inf = lqa <= hqa
            # reads a SINE-carrying haplotype yields relative to an empty one: one read window of
            # length L - 2W per informative junction, overlapping when reads are longer than the SINE
            ratio = 1.0 if m < 2 else (w + min(j2 - j1, w)) / w
            if m == 0 or not a_inf or k + a < args.min_reads:
                gt, gq = nocall, 0
            else:
                ll = []
                for fr in fs:
                    pp = ratio * fr / (ratio * fr + (1 - fr))
                    pp = min(1 - e, max(e, pp))
                    ll.append(k * math.log(pp) + a * math.log(1 - pp))
                best = max(range(len(ll)), key=lambda i: ll[i])
                second = max(ll[i] for i in range(len(ll)) if i != best)
                gt = labels[best]
                gq = min(99, int(round(10 * (ll[best] - second) / math.log(10))))
            tally[gt] += 1
            out.write(f'{g}\t{k}\t{a}\t{m}\t{gt}\t{gq}\n')
    print(f'{len(groups)} groups (mean read length {mean_len:.0f}): ' + ', '.join(f'{k} {v}' for k, v in sorted(tally.items())),
          file=sys.stderr)


def cmd_call(args):
    tpl, groups = load_templates(args.templates)
    counts, n, s = count_sam(args, tpl)
    genotype(args, tpl, groups, counts, s / n if n else 150, args.out)


def cmd_count(args):
    tpl, _ = load_templates(args.templates)
    counts, n, s = count_sam(args, tpl)
    with open(args.out, 'w') as out:
        out.write(f'#reads_sampled={n}\tbases_sampled={s}\n')
        out.write('group\t' + '\t'.join(COUNT_COLS) + '\n')
        for g, c in sorted(counts.items()):
            out.write(g + '\t' + '\t'.join(map(str, c)) + '\n')
    print(f'{len(counts)} groups with evidence -> {args.out}', file=sys.stderr)


def cmd_genotype(args):
    tpl, groups = load_templates(args.templates)
    counts = collections.defaultdict(lambda: [0] * len(COUNT_COLS))
    n = s = 0
    for path in args.counts:
        with open_text(path) as fh:
            meta = dict(kv.split('=') for kv in fh.readline()[1:].strip().split('\t'))
            n += int(meta['reads_sampled'])
            s += int(meta['bases_sampled'])
            fh.readline()
            for line in fh:
                f = line.rstrip('\n').split('\t')
                c = counts[f[0]]
                for i, v in enumerate(f[1:]):
                    c[i] += int(v)
    genotype(args, tpl, groups, counts, s / n if n else 150, args.out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    t = sub.add_parser('templates', help='write allele templates for every group')
    t.add_argument('--groups', required=True, help='PREFIX.groups.tsv of a pan-SINEome build')
    t.add_argument('--genome', nargs='+', required=True, metavar='SP=FASTA', help='genomes of the build (with .fai)')
    t.add_argument('--flank', type=int, default=300, help='flank length on each side [300]')
    t.add_argument('-o', '--out', required=True, help='output prefix')
    t.set_defaults(func=cmd_templates)
    def common(p, sam=True):
        p.add_argument('--templates', required=True, help='PREFIX.tsv from templates')
        if sam:
            p.add_argument('--sam', required=True, help='reads of one sample mapped to PREFIX.fa (SAM, may be gzipped)')
            p.add_argument('--anchor', type=int, default=30, help='min aligned bp on each side of a junction [30]')
            p.add_argument('--min-mapq', type=int, default=20, help='min mapping quality [20]')
            p.add_argument('--max-del', type=int, default=30, help='deletions up to this length do not break a span [30]')
            p.add_argument('--tol', type=int, default=30, help='max distance of a SINE-sized indel from the junction [30]')
            p.add_argument('--min-ins', type=int, default=100, help='min insertion length counted as a SINE [100]')
        else:
            p.add_argument('--anchor', type=int, default=30, help='as used for counting [30]')
        p.add_argument('-o', '--out', required=True)

    def gt_opts(p):
        p.add_argument('--min-reads', type=int, default=2, help='min reads to call a genotype [2]')
        p.add_argument('--error', type=float, default=0.02, help='per-read error rate in the likelihood [0.02]')
        p.add_argument('--ploidy', type=int, choices=(1, 2), default=2)

    c = sub.add_parser('call', help='genotype a sample from reads mapped to the templates (count + genotype)')
    common(c); gt_opts(c); c.set_defaults(func=cmd_call)
    k = sub.add_parser('count', help='per-group read counts of one SAM (to combine several SAMs of a sample)')
    common(k); k.set_defaults(func=cmd_count)
    y = sub.add_parser('genotype', help='genotype a sample from one or more count files')
    common(y, sam=False); gt_opts(y)
    y.add_argument('counts', nargs='+', help='count files of one sample')
    y.set_defaults(func=cmd_genotype)
    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
