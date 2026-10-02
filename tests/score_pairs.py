#!/usr/bin/env python3
"""Score orth_<a>-<b>.tsv against the simulation truth: both loci of a row must be the same event.
Usage: score_pairs.py SIMDIR RUNDIR sp1 sp2"""
import sys, re, collections, csv
S, run, a, b = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
CL={"root":"abcd","ab":"ab","cd":"cd","a":"a","b":"b","c":"c","d":"d"}
copies=collections.defaultdict(list); sites=collections.defaultdict(list); present={}
for sp in (a,b):
    for l in open(f"{S}/{sp}-SINEX.bed"):
        ch,s,e,n,_,st=l.rstrip().split("\t"); s,e=int(s),int(e)
        copies[(sp,ch)].append((s if st=="+" else e, n))
    for l in open(f"{S}/{sp}.sites.tsv"):
        n,ch,cur,st=l.rstrip().split("\t")
        if sp[0] not in CL[n.split("_")[1]]: sites[(sp,ch)].append((int(cur),n))
rows=list(csv.DictReader(open(f"{run}/orth_{a}-{b}.tsv"),delimiter="\t"))
L=collections.Counter(int(m.group(3))-int(m.group(2)) for r in rows for k in ("locus1","locus2") for m in [re.match(r"(.*):(\d+)-(\d+)\((.)\)",r[k])])
slop=L.most_common(1)[0][0]-300
def ev(sp,locus,sine):
    m=re.match(r"(.*):(\d+)-(\d+)\((.)\)",locus); ch,s,e,st=m.group(1),int(m.group(2)),int(m.group(3)),m.group(4)
    anc=e-slop if st=="+" else s+slop
    pool=copies[(sp,ch)] if sine=="1" else sites[(sp,ch)]
    best=min(pool,key=lambda x:abs(x[0]-anc),default=None)
    return best[1] if best and abs(best[0]-anc)<=80 else None
res=collections.Counter(); found=set()
for r in rows:
    e1=ev(r["species1"],r["locus1"],r["sine1"]); e2=ev(r["species2"],r["locus2"],r["sine2"])
    if e1 and e2 and e1.replace("_dup","")==e2.replace("_dup",""):
        if "_dup" in e1+e2: res["paralog(dup)"]+=1
        else: res["correct"]+=1; found.add(e1)
    elif e1 is None or e2 is None: res["unmatched"]+=1
    else: res["wrong"]+=1
truth=[l.split("\t")[0] for l in open(f"{S}/truth.tsv")][1:]
det=[n for n in truth if a[0] in CL[n.split("_")[1]] or b[0] in CL[n.split("_")[1]]]
rep=[n for n in det if n.endswith("_rep")]
print(f"{a}-{b}: rows={len(rows)} correct={res['correct']} wrong={res['wrong']} paralog={res['paralog(dup)']} unmatched={res['unmatched']} | recall {len(found)}/{len(det)}  (repeat-flank {sum(n in found for n in rep)}/{len(rep)}, other {sum(n in found for n in det if n not in rep)}/{len(det)-len(rep)})")
