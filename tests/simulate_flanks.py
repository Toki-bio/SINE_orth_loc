#!/usr/bin/env python3
"""Toy for the flank-mapping test (mapper, flank length, one or both flanks): two genomes aaa, bbb with the
bank of simulate_bank.py and hard flanks, each event with truth.

Usage: simulate_flanks.py [--div 0.03] [--events 450] [--seed 11]
  --div    substitutions per lineage in every sequence (flank divergence between aaa and bbb ~ 2 x div)
Event classes (truth.tsv `class`):
  plain      a SINE (full copy of A-D, or a partial piece of D) with unique flanks
  repflank   the 600 bp next to the SINE's 5' end are a copy of one young repeat (2 % from its consensus,
             ~70 copies in each genome): the left flank is not unique
  segdup     bbb only: 3 kb around the site duplicated (1 % apart) on contig bbb_dup; the true pair is the
             original, a pair with the duplicate is a paralog
  ctgend     bbb only: the contig ends 120-250 bp right of the site (3' side in SINE orientation)
  rfindel    bbb only: 200 bp of unrelated sequence inserted 60 bp into the right flank
where: both (insertion before the split), aaa only, bbb only.
Writes aaa.bnk, bbb.bnk, aaa-bank.bed, bbb-bank.bed, bank.fa, truth.tsv (event, class, label, aaa, bbb,
aaa_site, bbb_site, bbb_dup_site; site = contig:insertion junction = copy start on +, end on -)."""
import argparse, random
ap = argparse.ArgumentParser()
ap.add_argument('--div', type=float, default=0.03)
ap.add_argument('--events', type=int, default=450)
ap.add_argument('--seed', type=int, default=11)
a = ap.parse_args()
random.seed(a.seed)
B = "ACGT"
def rnd(n): return "".join(random.choice(B) for _ in range(n))
def mut(s, r): return "".join(random.choice([b for b in B if b != c]) if random.random() < r else c for c in s)
def rc(s): return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

random.seed(7)                                         # the bank of simulate_bank.py
A = rnd(110); Bc = rnd(250); C = A + Bc; D = rnd(300) + Bc + rnd(50)
random.seed(a.seed)
BANK = {"A": A, "B": Bc, "C": C, "D": D}
open("bank.fa", "w", newline="\n").write("".join(f">{k}\n{v}\n" for k, v in BANK.items()))
PIECES = [(1, 90), (120, 230), (451, 580)]
REP = rnd(600)

ev = []
for i in range(a.events):
    cls = random.choices(["plain", "repflank", "segdup", "ctgend", "rfindel"], [60, 15, 10, 10, 5])[0]
    where = random.choice(["both", "both", "aaa", "bbb"])
    if random.random() < 0.15:
        x, y = random.choice(PIECES); lab = f"D:{x}-{y}"; seq = D[x - 1:y]
    else:
        lab = random.choice("ABCD"); seq = BANK[lab]
    e = dict(id=f"e{i}", cls=cls, where=where, lab=lab, ins=mut(seq, 0.10 if where == "both" else 0.05),
             strand=random.choice("+-"), tsd=rnd(random.randint(8, 14)), pre=rnd(1250), post=rnd(1250))
    if cls == "repflank":                              # 5' side of the SINE = left of + copies, right of - copies
        r = mut(REP, 0.02)
        if e["strand"] == "+": e["pre"] = e["pre"][:-600] + r
        else: e["post"] = rc(r) + e["post"][600:]
    e["junk"] = rnd(200)
    ev.append(e)

N_CHR = 3
per = (a.events + N_CHR - 1) // N_CHR
out = {"aaa": {}, "bbb": {}}; bed = {"aaa": [], "bbb": []}; site = {}; dup = {}
for c in range(N_CHR):
    evs = ev[c * per:(c + 1) * per]
    for sp in ("aaa", "bbb"):
        seq, pos, cuts = [], 0, []
        for e in evs:
            pre, post = mut(e["pre"], a.div), mut(e["post"], a.div)
            if sp == "bbb" and e["cls"] == "rfindel":    # right flank = 3' side: post on +, pre on -
                if e["strand"] == "+": post = post[:60] + e["junk"] + post[60:]
                else: pre = pre[:-60] + e["junk"] + pre[-60:]
            t = mut(e["tsd"], a.div)
            if e["where"] in ("both", sp):
                ins = mut(e["ins"], a.div); ins = ins if e["strand"] == "+" else rc(ins)
                mid = t + ins + t
                st = pos + len(pre) + len(t); en = st + len(ins)
                bed[sp].append([c, st, en, e["lab"], e["strand"]])
                j = st if e["strand"] == "+" else en
            else:
                mid = t; j = pos + len(pre) + len(t)
            site[(e["id"], sp)] = [c, j]
            if sp == "bbb" and e["cls"] == "ctgend":
                k = random.randint(120, 250)
                cuts.append(pos + len(pre) + len(mid) + k if e["strand"] == "+" else pos + len(pre) - k)
            seq += [pre, mid, post]; pos += len(pre) + len(mid) + len(post)
        out[sp][c] = ("".join(seq), sorted(cuts))

# contigs: split at the cuts (100 bp of the sequence after a cut are lost, as in an assembly gap)
contigs = {"aaa": {}, "bbb": {}}; where_is = {}
for sp in ("aaa", "bbb"):
    for c, (s, cuts) in out[sp].items():
        bounds = [0] + [x for x in cuts] + [len(s)]
        for k in range(len(bounds) - 1):
            lo = bounds[k] + (100 if k else 0); hi = bounds[k + 1]
            name = f"{sp}_chr{c + 1}" + (f"_c{k + 1}" if len(bounds) > 2 else "")
            contigs[sp][name] = s[lo:hi]; where_is[(sp, c, k)] = (name, lo, hi)
def place(sp, c, x):
    for (s2, c2, k), (name, lo, hi) in where_is.items():
        if s2 == sp and c2 == c and lo <= x < hi: return name, x - lo
    return None, None

# segmental duplicates in bbb
dupseq = []
for e in ev:
    if e["cls"] != "segdup": continue
    c, j = site[(e["id"], "bbb")]; s = out["bbb"][c][0]
    lo, hi = max(0, j - 1500), min(len(s), j + 1500)
    off = sum(len(x) for x in dupseq) + 2000
    dupseq += [rnd(2000), mut(s[lo:hi], 0.01)]
    dup[e["id"]] = ("bbb_dup", off + (j - lo))
    for b in list(bed["bbb"]):
        if b[0] == c and lo <= b[1] and b[2] <= hi:
            bed["bbb"].append(["dup", off + b[1] - lo, off + b[2] - lo, b[3], b[4]])
if dupseq: contigs["bbb"]["bbb_dup"] = "".join(dupseq) + rnd(2000)

for sp in ("aaa", "bbb"):
    with open(f"{sp}.bnk", "w", newline="\n") as fh:
        for n, q in contigs[sp].items():
            fh.write(f">{n}\n" + "\n".join(q[i:i + 80] for i in range(0, len(q), 80)) + "\n")
    rows = []
    for c, st, en, lab, strand in bed[sp]:
        if c == "dup": n, s0, e0 = "bbb_dup", st, en
        else:
            n, s0 = place(sp, c, st); n2, e0 = place(sp, c, en - 1)
            if n is None or n != n2: continue
            e0 += 1
        rows.append((n, s0, e0, lab, strand))
    with open(f"{sp}-bank.bed", "w", newline="\n") as fh:
        for r in sorted(rows):
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t0\t{r[4]}\n")
with open("truth.tsv", "w", newline="\n") as fh:
    fh.write("event\tclass\tlabel\taaa\tbbb\taaa_site\tbbb_site\tbbb_dup_site\n")
    for e in ev:
        s = []
        for sp in ("aaa", "bbb"):
            n, x = place(sp, *site[(e["id"], sp)]); s.append(f"{n}:{x}" if n else "-")
        d = f"{dup[e['id']][0]}:{dup[e['id']][1]}" if e["id"] in dup else "-"
        fh.write(f"{e['id']}\t{e['cls']}\t{e['lab']}\t{int(e['where'] in ('both','aaa'))}\t{int(e['where'] in ('both','bbb'))}"
                 f"\t{s[0]}\t{s[1]}\t{d}\n")
