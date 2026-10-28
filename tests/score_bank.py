#!/usr/bin/env python3
"""Score orth_aaa-bbb.tsv of the bank toy (simulate_bank.py) against truth.tsv: both loci of a row must be the
same event (anchors within 40 bp of its sites) and the status must match the truth (SINE, PM or MP).
Usage: score_bank.py SIMDIR RUNDIR [refs_aaa-bbb.tsv]"""
import collections, csv, sys
S, run = sys.argv[1], sys.argv[2]
truth = list(csv.DictReader(open(f"{S}/truth.tsv"), delimiter="\t"))
sites = collections.defaultdict(list)
for t in truth:
    for sp in ("aaa", "bbb"):
        c, p = t[f"{sp}_site"].rsplit(":", 1); sites[c].append((int(p), t))
def ev(locus, anchor):
    if not anchor: return None
    c = locus.rsplit(":", 1)[0]; a = int(anchor)
    best = min(sites[c], key=lambda x: abs(x[0] - a), default=None)
    return best[1] if best and abs(best[0] - a) <= 40 else None
res = collections.Counter(); got = {}
for r in csv.DictReader(open(f"{run}/orth_aaa-bbb.tsv"), delimiter="\t"):
    if r["status"] not in ("SINE", "PM", "MP"): continue
    e1, e2 = ev(r["locus1"], r["anchor1"]), ev(r["locus2"], r["anchor2"])
    if not e1 or not e2: res["unmatched"] += 1; continue
    if e1 is not e2: res["wrong_pair"] += 1; continue
    p1 = e1["aaa"] if r["species1"] == "aaa" else e1["bbb"]; p2 = e1["bbb"] if r["species1"] == "aaa" else e1["aaa"]
    want = {("1", "1"): "SINE", ("1", "0"): "PM", ("0", "1"): "MP"}[(p1, p2)]
    ok = r["status"] == want; res["right_status" if ok else "wrong_status"] += 1
    got[e1["event"]] = got.get(e1["event"], False) or ok
print("rows:", dict(res))
by = collections.defaultdict(collections.Counter)
for t in truth:
    k = (t["class"], "both" if t["aaa"] == t["bbb"] == "1" else "one")
    by[k]["events"] += 1; by[k]["right"] += got.get(t["event"]) is True; by[k]["wrong_status"] += got.get(t["event"]) is False
for k in sorted(by): print(f"{k[0]:8s} {k[1]:5s} events {by[k]['events']:4d}  right {by[k]['right']:4d}  only wrong status {by[k]['wrong_status']:3d}")
