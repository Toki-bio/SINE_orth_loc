#!/usr/bin/env python3
"""Bank-mode toy: two genomes (aaa, bbb) with copies of four consensi of different length, as in a SINEderella bank.

Bank (bank.fa): A 110 bp, B 250 bp, C = A + B (composite, 360 bp), D 600 bp (B at D positions 301-550).
Events (truth.tsv), one every 2.5 kb:
  class   full      a whole copy of A, B, C or D
          partial   a piece of D (positions 1-90, 120-230 or 451-580), labelled D:a-b as in SINEderella's
                    family_view/partial_copies.tsv
          swap      a C copy that bbb's annotation labels B (an assignment that differs between species)
  where   both (inserted before the split), aaa only, bbb only
Shared copies are 10 % from their consensus before the split; every sequence then diverges 3 % per lineage.
Writes aaa.bnk, bbb.bnk, aaa-bank.bed, bbb-bank.bed, bank.fa, truth.tsv (event, class, label, aaa, bbb, aaa_site,
bbb_site: chrom:insertion junction = BED start of the copy on +, end on -, or the empty site), and the single-consensus
inputs D.fa, aaa-D.bed, bbb-D.bed for the comparison run."""
import random
random.seed(7)
B = "ACGT"
def rnd(n): return "".join(random.choice(B) for _ in range(n))
def mut(s, r): return "".join(random.choice([b for b in B if b != c]) if random.random() < r else c for c in s)
def rc(s): return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

A = rnd(110); Bc = rnd(250); C = A + Bc; D = rnd(300) + Bc + rnd(50)
BANK = {"A": A, "B": Bc, "C": C, "D": D}
open("bank.fa", "w", newline="\n").write("".join(f">{k}\n{v}\n" for k, v in BANK.items()))
open("D.fa", "w", newline="\n").write(f">D\n{D}\n")
PIECES = [(1, 90), (120, 230), (451, 580)]

events = []
for i in range(3 * 110):
    cls = random.choice(["full", "full", "full", "partial", "swap"])
    where = random.choice(["both", "both", "aaa", "bbb"])
    if cls == "full":
        lab = random.choice("ABCD"); seq = BANK[lab]
    elif cls == "partial":
        a, b = random.choice(PIECES); lab = f"D:{a}-{b}"; seq = D[a - 1:b]
    else:
        lab = "C"; seq = C
    events.append(dict(id=f"e{i}", cls=cls, where=where, lab=lab, seq=seq, strand=random.choice("+-"),
                       tsd=rnd(random.randint(8, 14))))
anc_copy = {e["id"]: mut(e["seq"], 0.10) for e in events}

genome = {"aaa": {}, "bbb": {}}; bed = {"aaa": [], "bbb": []}; bedD = {"aaa": [], "bbb": []}; site = {}
for chrom in range(1, 4):
    flanks = [rnd(2500) for _ in range(111)]
    evs = events[(chrom - 1) * 110: chrom * 110]
    for sp in ("aaa", "bbb"):
        s = [mut(flanks[0], 0.03)]; pos = len(s[0])
        for k, e in enumerate(evs):
            has = e["where"] in ("both", sp)
            if has:
                ins = anc_copy[e["id"]] if e["where"] == "both" else mut(e["seq"], 0.05)
                ins = mut(ins, 0.03)
                if e["strand"] == "-": ins = rc(ins)
                t = mut(e["tsd"], 0.03)
                piece = t + ins + t
                st, en = pos + len(t), pos + len(t) + len(ins)
                lab = "B" if (e["cls"] == "swap" and sp == "bbb") else e["lab"]
                bed[sp].append((f"{sp}_chr{chrom}", st, en, lab, e["strand"]))
                bedD[sp].append((f"{sp}_chr{chrom}", st, en, "D", e["strand"]))
                site[(e["id"], sp)] = f"{sp}_chr{chrom}:" + str(st if e["strand"] == "+" else en)
            else:
                piece = mut(e["tsd"], 0.03)
                site[(e["id"], sp)] = f"{sp}_chr{chrom}:{pos + len(piece)}"
            s.append(piece); pos += len(piece)
            f = mut(flanks[k + 1], 0.03); s.append(f); pos += len(f)
        genome[sp][f"{sp}_chr{chrom}"] = "".join(s)

for sp in ("aaa", "bbb"):
    with open(f"{sp}.bnk", "w", newline="\n") as fh:
        for c, q in genome[sp].items():
            fh.write(f">{c}\n" + "\n".join(q[i:i + 80] for i in range(0, len(q), 80)) + "\n")
    for name, rows in ((f"{sp}-bank.bed", bed[sp]), (f"{sp}-D.bed", bedD[sp])):
        with open(name, "w", newline="\n") as fh:
            for c, st, en, lab, strand in rows:
                fh.write(f"{c}\t{st}\t{en}\t{lab}\t0\t{strand}\n")
with open("truth.tsv", "w", newline="\n") as fh:
    fh.write("event\tclass\tlabel\taaa\tbbb\taaa_site\tbbb_site\n")
    for e in events:
        fh.write(f"{e['id']}\t{e['cls']}\t{e['lab']}\t{int(e['where'] in ('both','aaa'))}\t{int(e['where'] in ('both','bbb'))}"
                 f"\t{site[(e['id'],'aaa')]}\t{site[(e['id'],'bbb')]}\n")
