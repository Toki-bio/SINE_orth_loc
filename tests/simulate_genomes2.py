#!/usr/bin/env python3
"""Simulation v2: the v1 tree ((aaa,bbb),(ccc,ddd)) plus the hard cases, each with truth.

Event classes (truth.tsv column `class`):
  plain       SINE insertion on one branch
  repflank    left flank (5' in SINE orientation) is a young repeat (as in v1)
  flankindel  plain SINE + an unrelated 80-300 bp insertion 30-200 bp into its right flank on another branch
  host        old SINE (root) that receives a nested insertion
  nested      young SINE inserted into a host SINE on a later branch (TSD duplicated inside the host)
  close       second SINE inserted 40-250 bp downstream of a plain one, on its own branch
  satellite   tandem array of SINE+spacer units; unit number differs between lineages
Assembly artefacts:
  ddd: chromosomes cut into contigs 120-250 bp right of 18 insertion sites (contig_ends.tsv)
  ccc: a 30 kb region of chromosome 1 assembled twice (ccc_chr9, 0.4% divergent; haplotig.tsv)
Subfamilies S1..S4 (5% apart) are used by branch age (root S1/S2, ab/cd S2/S3, tips S3/S4); hosts are S1,
nested inserts S3/S4, so the true chronology is S1 < S2 < S3 < S4.
Annotation (<sp>-SINEX.bed) mimics sear2k: a SINE piece is reported only if >= 80% of the consensus length.
Writes <sp>.bnk, <sp>-SINEX.bed, SINEX.fa, subfamilies.fa, truth.tsv, <sp>.sites.tsv, nest_truth.tsv,
contig_ends.tsv, haplotig.tsv."""
import random
random.seed(22)
B = "ACGT"
def rnd(n): return "".join(random.choice(B) for _ in range(n))
def mut(s, r): return "".join(random.choice([b for b in B if b != c]) if random.random() < r else c for c in s)
def rc(s): return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

SP = ["aaa", "bbb", "ccc", "ddd"]
CLADES = {"root": SP, "ab": ["aaa", "bbb"], "cd": ["ccc", "ddd"], "a": ["aaa"], "b": ["bbb"], "c": ["ccc"], "d": ["ddd"]}
AGE = {"root": 0, "ab": 1, "cd": 1, "a": 2, "b": 2, "c": 2, "d": 2}
LATER = {"root": ["ab", "cd", "a", "b", "c", "d"], "ab": ["a", "b"], "cd": ["c", "d"]}
BR = list(CLADES)

base = rnd(260)
SUB = {f"S{i}": mut(base, 0.05) for i in range(1, 5)}
open("SINEX.fa", "w").write(">SINEX\n" + base + "\n")
open("subfamilies.fa", "w").write("".join(f">{k}\n{v}\n" for k, v in SUB.items()))
SUB_BY_AGE = {0: ["S1", "S2"], 1: ["S2", "S3"], 2: ["S3", "S4"]}
REPX = rnd(300)
MIN_REPORT = int(0.8 * len(base))

def sine_copy(sub, div): return mut(SUB[sub], div)

# ------------------------------------------------------------------ events on the ancestral genome
anc = {}; events = []; eid = 0
CLASSES = ["plain", "repflank", "flankindel", "plain", "host", "close", "plain", "satellite", "flankindel", "host"]
for c in range(1, 4):
    L = 420000; a = list(rnd(L))
    for p in range(4000, L - 4000, 2800):
        cls = CLASSES[eid % len(CLASSES)]
        if cls == "satellite" and c != 2:
            cls = "plain"                                    # satellites only on chromosome 2
        strand = "+" if cls == "satellite" else random.choice("+-")
        ev = dict(c=c, p=p, cls=cls, strand=strand, tsd=random.randint(8, 14), id=eid)
        if cls == "host":
            ev["br"] = "root"; ev["sub"] = "S1"; ev["ins"] = sine_copy("S1", 0.15)
            nb = random.choice(LATER["root"]); ny = random.choice(SUB_BY_AGE[AGE[nb]][1:] + ["S4"])
            k = random.randint(40, 220)                      # breakpoint inside the host (host orientation)
            ev["nest"] = dict(br=nb, sub=ny, ins=sine_copy(ny, 0.03), k=k, tsd=random.randint(8, 14),
                              strand=random.choice("+-"))
        elif cls == "satellite":
            ev["br"] = "root"; ev["sub"] = "S2"
            unit = sine_copy("S2", 0.08) + rnd(45)
            ev["units"] = {sp: 4 + (2 if sp in ("aaa", "bbb") else 0) + (1 if sp == "ddd" else 0) for sp in SP}
            ev["unit"] = unit
        else:
            br = BR[(eid * 7) % len(BR)]
            ev["br"] = br; ev["sub"] = random.choice(SUB_BY_AGE[AGE[br]]); ev["ins"] = sine_copy(ev["sub"], 0.10 - 0.03 * AGE[br])
        if cls == "repflank":                                 # left flank (SINE orientation) = REPX copy
            r = mut(REPX, 0.01)
            if strand == "+": a[p - 300:p] = list(r)
            else: a[p:p + 300] = list(rc(r))
        if cls == "flankindel":                               # unrelated insertion in the right flank
            ev["indel"] = dict(br=random.choice([b for b in BR if b != ev["br"]]), seq=rnd(random.randint(80, 300)),
                               d=random.randint(30, 200))
        if cls == "close":
            ev["close"] = dict(br=random.choice([b for b in BR if b != ev["br"]]), sub="S3", ins=sine_copy("S3", 0.04),
                               d=random.randint(40, 250), strand=random.choice("+-"), tsd=random.randint(8, 14))
        events.append(ev); eid += 1
    anc[c] = "".join(a)

def name(ev, part=""):
    return f"EV{ev['id']}{part}_{ev['br'] if not part else ev[part.strip('_')]['br']}"

# truth rows: one per element (insertion event)
truth = open("truth.tsv", "w"); truth.write("event\tbranch\tclass\tsubfamily\t" + "\t".join(SP) + "\n")
def trow(n, br, cls, sub, states=None):
    st = states or ["P" if s in CLADES[br] else "A" for s in SP]
    truth.write(f"{n}\t{br}\t{cls}\t{sub}\t" + "\t".join(st) + "\n")
nest_truth = open("nest_truth.tsv", "w"); nest_truth.write("host\tinsert\thost_sub\tinsert_sub\tinsert_branch\tbreakpoint\n")
for ev in events:
    n = f"EV{ev['id']}_{ev['br']}" + ("_rep" if ev["cls"] == "repflank" else "")
    ev["name"] = n
    if ev["cls"] == "satellite":
        trow(n, "root", "satellite", ev["sub"], [f"x{ev['units'][s]}" for s in SP]); continue
    trow(n, ev["br"], ev["cls"], ev["sub"])
    if ev["cls"] == "host":
        nn = f"EV{ev['id']}n_{ev['nest']['br']}"; ev["nest"]["name"] = nn
        trow(nn, ev["nest"]["br"], "nested", ev["nest"]["sub"])
        nest_truth.write(f"{n}\t{nn}\tS1\t{ev['nest']['sub']}\t{ev['nest']['br']}\t{ev['nest']['k']}\n")
    if ev["cls"] == "close":
        nn = f"EV{ev['id']}c_{ev['close']['br']}"; ev["close"]["name"] = nn
        trow(nn, ev["close"]["br"], "close", "S3")
    if ev["cls"] == "flankindel":
        trow(f"EV{ev['id']}i_{ev['indel']['br']}", ev["indel"]["br"], "otherTE", "-")

# ------------------------------------------------------------------ per-species genomes
def oriented(seq, strand): return seq if strand == "+" else rc(seq)

contig_cut = {}                                               # ddd: (chrom) -> list of cut positions
for sp in SP:
    seqs = {}; bed = []; sites = []
    for c in range(1, 4):
        a = anc[c]; out = []; last = 0
        nonlocal_cur = [0]                                   # running coordinate in this species
        # build a list of genomic actions sorted by ancestral position
        acts = []
        for ev in events:
            if ev["c"] != c: continue
            acts.append((ev["p"], "main", ev))
            if "indel" in ev:
                d = ev["indel"]["d"]; acts.append((ev["p"] + d if ev["strand"] == "+" else ev["p"] - d, "indel", ev))
            if "close" in ev:
                d = ev["close"]["d"]; acts.append((ev["p"] + d if ev["strand"] == "+" else ev["p"] - d, "close", ev))
        acts.sort(key=lambda x: x[0])
        for p, kind, ev in acts:
            out.append(a[last:p]); nonlocal_cur[0] += p - last; last = p
            pos = nonlocal_cur[0]
            if kind == "main":
                sites.append((ev["name"], f"{sp}_chr{c}", pos, ev["strand"]))
                if ev["cls"] == "satellite":
                    n = ev["units"][sp]; arr = ""
                    for u in range(n):
                        unit = mut(ev["unit"], 0.02)
                        bed.append((f"{sp}_chr{c}", pos + len(arr), pos + len(arr) + len(SUB['S2']), f"{ev['name']}_u{u}", ev["strand"]))
                        arr += unit
                    piece = arr + a[p - ev["tsd"]:p]       # satellites are built on the + strand
                    out.append(piece); nonlocal_cur[0] += len(piece); continue
                if sp not in CLADES[ev["br"]]:
                    continue
                host = ev["ins"]
                if ev["cls"] == "host":
                    nst = ev["nest"]; k = nst["k"]
                    if sp in CLADES[nst["br"]]:
                        y = oriented(nst["ins"], nst["strand"])
                        hpart = host[:k] + y + host[k - nst["tsd"]:k] + host[k:]
                        # host halves (host orientation) and the insert, in genome coordinates
                        h1 = (0, k); yi = (k, k + len(y)); h2 = (k + len(y), len(hpart))
                        segs = [("h1", h1), ("y", yi), ("h2", h2)]
                    else:
                        hpart = host; segs = [("h", (0, len(host)))]
                        bp = pos + k if ev["strand"] == "+" else pos + len(host) - k    # empty site of the insert
                        sites.append((nst["name"], f"{sp}_chr{c}", bp, nst["strand"] if ev["strand"] == "+" else ("-" if nst["strand"] == "+" else "+")))
                    g = oriented(hpart, ev["strand"]); L = len(g)
                    for tag, (s0, e0) in segs:
                        gs, ge = (pos + s0, pos + e0) if ev["strand"] == "+" else (pos + L - e0, pos + L - s0)
                        if tag == "y":
                            ystrand = nst["strand"] if ev["strand"] == "+" else ("-" if nst["strand"] == "+" else "+")
                            bed.append((f"{sp}_chr{c}", gs, ge, nst["name"], ystrand))
                        elif e0 - s0 >= MIN_REPORT:
                            bed.append((f"{sp}_chr{c}", gs, ge, ev["name"] + ("" if tag == "h" else f"_{tag}"), ev["strand"]))
                    piece = g + a[p - ev["tsd"]:p]
                else:
                    g = oriented(host, ev["strand"])
                    bed.append((f"{sp}_chr{c}", pos, pos + len(g), ev["name"], ev["strand"]))
                    piece = g + a[p - ev["tsd"]:p]
                out.append(piece); nonlocal_cur[0] += len(piece)
            elif kind == "indel":
                if sp in CLADES[ev["indel"]["br"]]:
                    out.append(ev["indel"]["seq"]); nonlocal_cur[0] += len(ev["indel"]["seq"])
            elif kind == "close":
                cl = ev["close"]
                sites.append((cl["name"], f"{sp}_chr{c}", pos, cl["strand"]))
                if sp in CLADES[cl["br"]]:
                    g = oriented(cl["ins"], cl["strand"])
                    bed.append((f"{sp}_chr{c}", pos, pos + len(g), cl["name"], cl["strand"]))
                    piece = g + a[p - cl["tsd"]:p]; out.append(piece); nonlocal_cur[0] += len(piece)
        out.append(a[last:]); seqs[f"{sp}_chr{c}"] = mut("".join(out), 0.02)
    # ---- assembly artefacts
    if sp == "ccc":                                          # haplotig: 30 kb of chr1 assembled twice
        s0, e0 = 100000, 130000
        seqs["ccc_chr9"] = rnd(2000) + mut(seqs["ccc_chr1"][s0:e0], 0.004) + rnd(2000)
        for b in list(bed):
            if b[0] == "ccc_chr1" and s0 <= b[1] and b[2] <= e0:
                bed.append(("ccc_chr9", b[1] - s0 + 2000, b[2] - s0 + 2000, b[3] + "_hap", b[4]))
        open("haplotig.tsv", "w").write(f"ccc\tccc_chr1\t{s0}\t{e0}\tccc_chr9\t2000\n")
    if sp == "ddd":                                          # contig breaks next to 18 sites
        plain_sites = [s for s in sites if not s[0].endswith("_rep")][5::12][:18]
        cuts = {}
        for n, ch, pos, st in plain_sites:
            off = random.randint(120, 250)
            cuts.setdefault(ch, []).append(pos + 360 + off if st == "+" else pos - off)
        new = {}; remap = {}
        for ch, s in seqs.items():
            cs = sorted(x for x in cuts.get(ch, []) if 0 < x < len(s)); starts = [0] + cs; ends = cs + [len(s)]
            for i, (x, y) in enumerate(zip(starts, ends)):
                nm = ch if i == 0 else f"{ch}_ctg{i}"
                new[nm] = s[x:y]; remap.setdefault(ch, []).append((x, y, nm))
        def rm(ch, pos):
            for x, y, nm in remap.get(ch, [(0, 10**12, ch)]):
                if x <= pos < y: return nm, pos - x
            return ch, pos
        bed2 = []
        for ch, s, e, n, st in bed:
            nm, ns = rm(ch, s); bed2.append((nm, ns, ns + (e - s), n, st))
        bed = bed2
        sites = [(n, *rm(ch, pos), st) for n, ch, pos, st in sites]
        seqs = new
        with open("contig_ends.tsv", "w") as fh:
            fh.write("event\n" + "".join(f"{n}\n" for n, *_ in plain_sites))
    with open(f"{sp}.bnk", "w") as fa:
        for ch, g in seqs.items():
            fa.write(f">{ch}\n" + "\n".join(g[i:i + 80] for i in range(0, len(g), 80)) + "\n")
    with open(f"{sp}-SINEX.bed", "w") as fh:
        for ch, s, e, n, st in sorted(bed, key=lambda x: (x[0], x[1])):
            fh.write(f"{ch}\t{s}\t{e}\t{n}\t0\t{st}\n")
    with open(f"{sp}.sites.tsv", "w") as fh:
        for n, ch, pos, st in sites:
            fh.write(f"{n}\t{ch}\t{pos}\t{st}\n")
cnt = {}
for ev in events: cnt[ev["cls"]] = cnt.get(ev["cls"], 0) + 1
print(len(events), "primary events:", cnt)
