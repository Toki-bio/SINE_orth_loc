#!/usr/bin/env python3
"""Simulate a diploid individual and its reads from the test genomes.

simulate_reads.py --hap aaa.bnk --hap ccc.bnk --mut 0.005 --cov 10 --type short|hifi --seed 1 -o out.fq
Each haplotype is a copy of the given genome with --mut independent substitutions (a new
individual, not the assembled one); reads are sampled uniformly from both haplotypes
(--cov is the total depth), on random strands, with substitution errors.
"""
import argparse
import random

B = 'ACGT'


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip().upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def mutate(s, r, rng):
    s = list(s)
    for i in range(len(s)):
        if rng.random() < r:
            s[i] = rng.choice([b for b in B if b != s[i]])
    return ''.join(s)


def rc(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hap', action='append', required=True)
    ap.add_argument('--mut', type=float, default=0.005)
    ap.add_argument('--cov', type=float, required=True)
    ap.add_argument('--type', choices=('short', 'hifi'), required=True)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('-o', '--out', required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    err = 0.005 if a.type == 'short' else 0.002
    n = 0
    with open(a.out, 'w') as out:
        for h, path in enumerate(a.hap):
            chroms = {k: mutate(v, a.mut, rng) for k, v in read_fasta(path).items()}
            for name, seq in chroms.items():
                L = len(seq)
                mean = 150 if a.type == 'short' else 15000
                nreads = int(a.cov / len(a.hap) * L / mean)
                for _ in range(nreads):
                    rl = 150 if a.type == 'short' else max(3000, int(rng.gauss(15000, 3000)))
                    if rl >= L:
                        continue
                    s = rng.randrange(0, L - rl)
                    r = seq[s:s + rl]
                    if rng.random() < 0.5:
                        r = rc(r)
                    r = ''.join(rng.choice([b for b in B if b != c]) if rng.random() < err else c for c in r)
                    n += 1
                    out.write(f'@r{n}_h{h}_{name}_{s}\n{r}\n+\n{"I" * len(r)}\n')


if __name__ == '__main__':
    main()
