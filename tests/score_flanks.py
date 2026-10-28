#!/usr/bin/env python3
"""Score orth_aaa-bbb.tsv of the flank toy (simulate_flanks.py) against truth.tsv.
A row is matched to events by its anchors (within 40 bp of the event's site in each genome; in bbb also the
segmental duplicate). Row outcomes: right (same event, right status), wrong_status, wrong_pair (different
events), paralog (paired with the duplicate), unmatched. Per class: events, right, wrong (any wrong row).
Usage: score_flanks.py SIMDIR RUNDIR [LABEL]   (prints one tab-separated summary line per class)"""
import collections, csv, sys
S, run = sys.argv[1], sys.argv[2]
label = sys.argv[3] if len(sys.argv) > 3 else run
truth = list(csv.DictReader(open(f"{S}/truth.tsv"), delimiter="\t"))
sites = collections.defaultdict(list)
for t in truth:
    for k, kind in (("aaa_site", "aaa"), ("bbb_site", "bbb"), ("bbb_dup_site", "dup")):
        if t[k] != "-":
            c, p = t[k].rsplit(":", 1); sites[c].append((int(p), t, kind))
def ev(locus, anchor):
    if not anchor: return None
    c = locus.rsplit(":", 1)[0]; x = int(anchor)
    b = min(sites[c], key=lambda s: abs(s[0] - x), default=None)
    return b if b and abs(b[0] - x) <= 40 else None
rows = collections.Counter(); right = set(); wrong = set()
for r in csv.DictReader(open(f"{run}/orth_aaa-bbb.tsv"), delimiter="\t"):
    if r["status"] not in ("SINE", "PM", "MP"): continue
    h1, h2 = ev(r["locus1"], r["anchor1"]), ev(r["locus2"], r["anchor2"])
    if not h1 or not h2: rows["unmatched"] += 1; continue
    e1, e2 = h1[1], h2[1]
    if e1 is not e2:
        rows["wrong_pair"] += 1; wrong.add(e1["event"]); wrong.add(e2["event"]); continue
    if "dup" in (h1[2], h2[2]): rows["paralog"] += 1; wrong.add(e1["event"]); continue
    p1, p2 = (e1["aaa"], e1["bbb"]) if r["species1"] == "aaa" else (e1["bbb"], e1["aaa"])
    want = {("1", "1"): "SINE", ("1", "0"): "PM", ("0", "1"): "MP"}[(p1, p2)]
    if r["status"] == want: rows["right"] += 1; right.add(e1["event"])
    else: rows["wrong_status"] += 1; wrong.add(e1["event"])
by = collections.defaultdict(collections.Counter)
for t in truth:
    for k in (t["class"], "ALL"):
        by[k]["events"] += 1; by[k]["right"] += t["event"] in right; by[k]["wrong"] += t["event"] in wrong
print(label, "rows", " ".join(f"{k}={v}" for k, v in sorted(rows.items())), sep="\t")
for k in sorted(by):
    print(label, k, f"events={by[k]['events']}", f"right={by[k]['right']}", f"wrong={by[k]['wrong']}", sep="\t")
