#!/usr/bin/env python3
"""Check sine_nest.py tint on a simulated insertion process with known activity periods.
Types T1..T6 insert at times ~ N(mu, sd) into a growing genome; an insertion landing inside an existing
element is a nesting event. The fitted chronology must recover the order of mu."""
import bisect, math, os, random, subprocess, sys, tempfile
def simulate(mus, sds, counts, seed):
    rng = random.Random(seed)
    ev = sorted((rng.gauss(m, s), t) for t, (m, s, c) in enumerate(zip(mus, sds, counts)) for _ in range(c))
    G = 5_000_000; L = 300
    starts = []; types = []; nest = []
    for _, t in ev:
        p = rng.randrange(G)
        i = bisect.bisect_right(starts, p) - 1
        if i >= 0 and p < starts[i] + L:
            nest.append((t, types[i]))
        # the new element shifts nothing in this simple model; it occupies [p, p+L)
        j = bisect.bisect_left(starts, p); starts.insert(j, p); types.insert(j, t)
    return nest
def run(mus, sds, counts, seed, label):
    nest = simulate(mus, sds, counts, seed)
    d = tempfile.mkdtemp()
    with open(f'{d}/n.tsv', 'w') as fh:
        fh.write('locus\tchrom\thost_type\tinsert_type\n')
        for a, b in nest:
            fh.write(f'x\tc\tT{b + 1}\tT{a + 1}\n')
    subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), '..', 'sine_nest.py'), 'tint', f'{d}/n.tsv',
                    '-o', f'{d}/t', '--boot', '50'], check=True, capture_output=True)
    rows = [l.rstrip('\n').split('\t') for l in open(f'{d}/t.chronology.tsv')][1:]
    est = {r[0]: float(r[1]) for r in rows}
    order = [r[0] for r in rows]
    truth = [f'T{i + 1}' for i in sorted(range(len(mus)), key=lambda i: mus[i])]
    # Spearman between fitted mu and true mu
    t_rank = {t: i for i, t in enumerate(truth)}; e_rank = {t: i for i, t in enumerate(order)}
    n = len(truth); rho = 1 - 6 * sum((t_rank[t] - e_rank[t]) ** 2 for t in truth) / (n * (n * n - 1))
    print(f'{label}: {len(nest)} nesting events; true order {" < ".join(truth)}; fitted {" < ".join(order)}; '
          f'Spearman {rho:.2f}; mu ' + ' '.join(f'{t}={est[t]:.2f}' for t in truth))
    return rho
ok = True
ok &= run([0, 1, 2, 3, 4, 5], [1] * 6, [6000] * 6, 1, 'equal widths, equal copy numbers') > 0.99
ok &= run([0, 0.5, 2, 2.2, 4, 6], [1] * 6, [3000, 8000, 4000, 2000, 6000, 3000], 2, 'overlapping, unequal copy numbers') > 0.9
ok &= run([0, 1, 2, 3, 4, 5], [0.5, 2, 1, 0.7, 1.5, 1], [5000] * 6, 3, 'unequal activity widths (model violated)') > 0.8
ok &= run([0, 1, 2, 3, 4, 5], [1] * 6, [400] * 6, 4, 'sparse: ~100 nesting events') > 0.7
sys.exit(0 if ok else 1)
