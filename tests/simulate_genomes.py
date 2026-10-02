#!/usr/bin/env python3
"""Simulate 4 genomes on the tree ((aaa,bbb),(ccc,ddd)) with SINE insertions on every branch,
a young repeat as left flank of every 4th locus, and 10 segmental duplications in aaa.
Writes <sp>.bnk, <sp>-SINEX.bed, SINEX.fa, truth.tsv and <sp>.sites.tsv in the current directory."""
import random
random.seed(21)
B="ACGT"
def rnd(n): return "".join(random.choice(B) for _ in range(n))
def mut(s,r): return "".join(random.choice([b for b in B if b!=c]) if random.random()<r else c for c in s)
def rc(s): return s[::-1].translate(str.maketrans("ACGT","TGCA"))
SP=["aaa","bbb","ccc","ddd"]
CLADES={"root":SP,"ab":["aaa","bbb"],"cd":["ccc","ddd"],"a":["aaa"],"b":["bbb"],"c":["ccc"],"d":["ddd"]}
BR=list(CLADES)
sine=rnd(260); open("SINEX.fa","w").write(">SINEX\n"+sine+"\n")
REPX=rnd(300)                       # young repeat occupying whole left flanks of some loci
anc={}; events=[]; eid=0
for c in range(1,4):
    L=400000; a=list(rnd(L))
    for p in range(4000,L-4000,2600):
        br=BR[eid%len(BR)]; ins=mut(sine,0.10); strand=random.choice("+-")
        flankrep = (eid%4==1)
        if flankrep:                # left flank (upstream in SINE orientation) = REPX copy
            r=mut(REPX,0.01)
            if strand=="+": a[p-300:p]=list(r)
            else:           a[p:p+300]=list(rc(r))
        if strand=="-": ins=rc(ins)
        events.append(dict(c=c,p=p,br=br,ins=ins,strand=strand,name=f"EV{eid}_{br}"+("_rep" if flankrep else ""),tsd=random.randint(8,14)))
        eid+=1
    anc[c]="".join(a)
SDS=[e for e in events if "aaa" in CLADES[e["br"]]][3:60:6]   # duplicated in aaa lineage
truth=open("truth.tsv","w"); truth.write("event\tbranch\t"+"\t".join(SP)+"\n")
for e in events: truth.write(f"{e['name']}\t{e['br']}\t"+"\t".join("P" if s in CLADES[e['br']] else "A" for s in SP)+"\n")
for sp in SP:
    fa=open(f"{sp}.bnk","w"); bed=open(f"{sp}-SINEX.bed","w"); seqs={}; coords={}; st=open(f"{sp}.sites.tsv","w")
    for c in range(1,4):
        a=anc[c]; out=[]; last=0; cur=0
        for e in events:
            if e["c"]!=c: continue
            p=e["p"]; out.append(a[last:p]); cur+=p-last; last=p
            st.write(f"{e['name']}\t{sp}_chr{c}\t{cur}\t{e['strand']}\n")
            if sp in CLADES[e["br"]]:
                bed.write(f"{sp}_chr{c}\t{cur}\t{cur+len(e['ins'])}\t{e['name']}\t0\t{e['strand']}\n")
                coords[e["name"]]=(c,cur); piece=e["ins"]+a[p-e["tsd"]:p]
                out.append(piece); cur+=len(piece)
        out.append(a[last:]); seqs[c]=mut("".join(out),0.02)
    if sp=="aaa":                   # segmental duplications onto a new contig
        contig=""; 
        for e in SDS:
            c,pos=coords[e["name"]]; seg=seqs[c][pos-1500:pos+len(e["ins"])+1500]
            base=len(contig)+3000; contig+=rnd(3000)+mut(seg,0.005)
            bed.write(f"aaa_chr9\t{base+1500}\t{base+1500+len(e['ins'])}\t{e['name']}_dup\t0\t{e['strand']}\n")
        seqs[9]=contig+rnd(3000)
    for c,g in seqs.items(): fa.write(f">{sp}_chr{c}\n"+"\n".join(g[i:i+80] for i in range(0,len(g),80))+"\n")
print(len(events),"events,",sum(1 for e in events if e["name"].endswith("_rep")),"with repeat flanks,",len(SDS),"segmental duplications in aaa")
