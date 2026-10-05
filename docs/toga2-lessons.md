# What SINE_orth_loc can borrow from TOGA2 — and where it is better

Comparison of SINE_orth_loc (branch `claude/wonderful-tesla-x39s12`, da26c47) with TOGA2 v2.0.10
(source and wiki; see `.claude/skills/toga/`). Numbers marked *measured* come from the dva–mix rerun
(`darevskia-sineome-data/dva-mix-v2`, code a53c4eb) and dar7 (`pan/`).

The two tools solve the same core problem — *is this locus in genome B the orthologous position of
this feature in genome A?* — for different features (coding genes vs SINE insertions) with opposite
designs: TOGA2 projects a reference through whole-genome alignment chains and decides orthology from
**long-range context** (10 kb flanks, introns, neighbouring genes); SINE_orth_loc maps **300 bp of
5′ flank** of every copy in both genomes and decides from **local identity** of the aligned flanks
and the insertion site.

## Status (v2.2)

Implemented and tested (see README "Nested, close, satellite and edge cases"): **2** Missing state
(`contig_end` → `MISSING` rows → registry state M), **3** flank-indel masking (`FI=1`), **6** compound loci
(`sine_nest.py scan/orth/supersede`, extended to nested SINE-in-SINE with TSD-verified host reassembly
and satellites), exact alignment anchors in the orth table (`--precise-tol`), and a TinT-style
chronology. **4** was examined with a real read-labelled set instead of TOGA-style weak labels:
ComPair calls are 99.2–99.8% consistent with reads, and a metric-based model ranks the rest with
AUC 0.70 — useful for triage, not as a replacement classifier. Open: **1** neighbour-anchor score,
**5** weighted graph splitting, **7** split-contig sites, **8** shared flank index.

## Where SINE_orth_loc is objectively better

| | SINE_orth_loc | TOGA2 | why it matters |
|---|---|---|---|
| Direction | symmetric: copies of **both** genomes are queried in every pair; pan-SINEome has no reference | reference → query projection only | presence/absence characters need both directions. *Measured*: 12,786 dva−/mix+ loci (MP) in dva–mix are, by construction, invisible to a dva-referenced projection |
| Evidence for absence | an empty site requires both flanks aligned (≥ 150 columns, > 65% id, right flank re-tested after trimming) and the element absent between them | a "Deleted" exon = not found under weak criteria (≥ 45% id / BLOSUM ≥ 20 / chain overlap / SpliceAI) | absence is the informative state for insertion phylogenetics; SINE_orth_loc requires positive evidence of the empty site, TOGA infers absence from failure to find |
| Prerequisite | genome FASTA + element BED; BWA + MAFFT | pairwise LASTZ chains (`make_lastz_chains`) | whole-genome alignment is usually the dominant compute; *measured*: full dva–mix run 12,624 s on one server |
| Repeats | built for TE junctions | coding genes; chains seeded on soft-masked genomes | the object of study lives where whole-genome chains are weakest |
| Identifiers | stable `PSG` IDs with merged/split/retired aliases; *measured*: identical rebuild of dar7 keeps every ID (4640840) | projection names embed chain IDs, which change with every re-alignment | a growing, cited resource needs stable IDs |
| Individuals without assemblies | `sine_genotype.py` genotypes groups from reads; *measured*: 99.6% agreement (343/92,620 contradictions) on D. valentini 245 vs its own assembly | needs an assembly per individual | population sampling, parthenogen haplotypes |
| Multi-genome structure | pan-SINEome groups + anchor graph (edges, breakpoints, duplicated blocks) across all genomes | synteny is a classifier feature; orthogroups are a union-find afterthought | rearrangement / haplotig detection comes out directly |
| Transparency | every metric per locus, thresholds tunable live in `sine_loci_browser.html`, every copy and cluster accounted (`copies.tsv`, `clusters_*.tsv`, `rescue_*.tsv`) | fixed XGBoost model, rule trees; rejection logs | users can see *why* each call was made and change the criterion |
| Footprint | ~2.2 k lines bash/awk/Python stdlib | ~50 k lines Python/Rust/Cython/C/Nextflow/HDF5, IQ-TREE, PRANK, SpliceAI | install, audit, maintenance |

## Where TOGA2 is objectively better (and what to borrow)

Ranked by expected gain for the Darevskia data.

### 1. Long-range context for ambiguous placements (TOGA: `synt`, `flank_cov`)
TOGA's decisive features are *context*: how much of a 10 kb neighbourhood aligns and how many
neighbouring genes the same chain covers. SINE_orth_loc sees 300 bp (two-flank rescue: two × 300 bp).
*Measured*: 5,273 copies end `ambiguous` in the dva–mix rescue (3,365 mix copies with exactly two equal
placements in dva), 10,242 `no_concordant_pair`.

Borrow — a **neighbour-anchor score**, computable from data we already have: for a candidate placement
P of copy c, count validated pan-SINEome anchors (groups already paired between the two genomes) within
±50 kb of c in genome A whose partners lie within ±50 kb of P in genome B, in the same order. The
orthologous placement inherits its neighbours; a paralogous one does not. This is TOGA's synteny feature
with SINE anchors instead of genes.
Caveat: it cannot break ties for **uncollapsed haplotigs** (both copies carry the neighbours —
*measured*: 55% of tested dva double placements are ≥ 95% identical over ≥ 80% of 10 kb); those should
become an explicit 1:many state (below), not a forced choice.

### 2. Missing ≠ absent (TOGA: Missing vs Lost/Deleted)
TOGA never calls a loss where the sequence is simply not there (assembly gap, chain/contig end).
SINE_orth_loc reports such loci as badLF/badRF/shortRF and they end as `U` or `no_validated_pair`.
*Measured* in dva–mix rejected alignments: 1,975 of 2,242 badLF, 1,503 of 16,354 badRF and 205 of
2,003 shortRF have a window shorter than the full flank + element + flank although the 300 bp flank mapped
full-length — most likely the window hits a contig boundary (to confirm with the `.fai` lengths, not available in the data repo). (No window contained ≥ 10 N; HiFi assemblies.)

Borrow — classify every locus that fails for lack of sequence as **M (missing)** with a reason:
`contig_end` (window clipped by the `.fai` length), `assembly_gap` (≥ 10 N in the window, N treated as
missing data rather than mismatch in identity), `no_hit`. Carry M into the registry as a state distinct
from `U` (no comparison) so that the matrix says *why* a cell is empty.

### 3. Mask large indels in the flank instead of failing (TOGA: BIG_INS/BIG_DEL masked)
*Measured*: 9,800 of 16,354 badRF loci (60%) contain a one-sided gap of ≥ 50 bp in the right flank
(a second insertion or a deletion next to the site). Re-scoring the right flank with those gap blocks
masked lets 2,145 pass the existing criteria (1,751 SINE, 394 PM/MP) — about +2.7% loci for dva–mix.

Borrow — in `ComPair.sh`, compute FR/FRcut on columns outside one-sided gap runs ≥ 50 bp, keep the
original metrics, and add a flag `flank_indel` so the call can be filtered.

### 4. A probability instead of a cliff (TOGA: XGBoost orthology probability)
ComPair.sh is a hard decision tree on 7 thresholds. TOGA turns features into a probability that later
steps use as edge weights. Borrow the idea, **not** TOGA's training recipe (its labels are weak: positives
= genes with a single chain, negatives = "non-top chains").

Better labels are available here: simulation truth (`tests/simulate_genomes.py`), conspecific pairs
(dva–dvl: true PM/MP should be rare), and read-validated loci (`sine_genotype.py` vs assembly). Features:
the ComPair metrics, MAPQ/AS/XS of both flank hits, number of alternative hits, neighbour-anchor score (1),
flank-indel flag (3). A logistic model is enough; keep the thresholds as the transparent default.

### 5. Weighted orthology graph with conservative splitting (TOGA: `split_graph`, second-best eviction)
`sine_registry.py` joins sites by union-find: one wrong pair merges two groups (→ `multicopy`/`X`).
*Measured* in dar7: 2,880 multicopy and 1,626 inconsistent groups. TOGA keeps edge weights (orthology
probability) and only cuts a many:many component when every removed edge is < 0.75 or < 0.9 × the
weakest kept edge, without separating strongly connected references; when several references claim one
locus it keeps the best ("second-best projection" eviction).

Borrow — weight registry edges by pair confidence (4); split multicopy components with TOGA's two rules;
when a site is claimed by two groups keep the stronger; report 1:many explicitly (state `P2`, or
`multicopy` with copy list) instead of collapsing to `X`.

### 6. Nested / close copies (TOGA: nested genes collapsed, not dropped)
SINE_orth_loc drops every copy within 300 bp of another (`bedtools cluster -d 300`, `uniq -u`).
*Measured*: 31,529 copies (4.4% of 714,365 in dar7) are `close_copy`. SINE-in-SINE and dimeric insertions
are phylogenetically informative (relative age within one locus).

Borrow — treat a cluster of close copies as one **compound locus**: flanks outside the cluster, one
alignment, per-element presence read from the alignment columns of each element.

### 7. Assemblies split across contigs (TOGA: fragment stitching, default on)
Left and right flank of one site on two contigs currently fail as shortRF/contig-end. With the two-flank
machinery in `rescue_multi.py`, allow the left flank on one contig end and the right flank on another
contig end → state `split_site` (the empty/occupied status is then unknown, but orthology is supported).

### 8. Linear scaling in genome number (TOGA: one run per query against a reference)
SINE_orth_loc runs all pairs, N(N−1)/2 (21 for 7 genomes, 190 for 20). TOGA scales as N−1 but pays with
reference-centrism. Keep symmetry, borrow the shape: one shared **index of representative flanks** of all
pan-SINEome groups; each new genome is mapped once against it and its own copies are added (already in
`docs/pan-sineome.md` "planned").

### 9. Engineering
- Resume from a step and re-run only failed batches (TOGA2 `--resume_from`, failed-batch lists); the bash
  pipeline restarts whole stages.
- Per-genome QC z-scores (TOGA `toga2orthogroups` OrthoZ/FamZ) on pan-SINEome outputs: fraction U/X/M,
  genome-specific breakpoints, duplicated blocks; *measured* outliers already visible (arm, nai breakpoints;
  dva duplicated blocks).
- Browser export as a UCSC track hub (BigBed per genome with group IDs and states), next to the HTML viewer.
- Benchmarks against external truth, as TOGA does with Ensembl/BUSCO: simulation, conspecific pairs,
  read-validated loci — reported with every release.

## What not to borrow

- CESAR2/SpliceAI/codon logic — coding-gene specific.
- Whole-genome chains as the *only* locus source: they would remove the cheap-prerequisite advantage.
  As an optional second locus source where chains already exist (Hiller-lab alignments of >3,000 pairs),
  projecting SINE junctions through chains would give TOGA-strength context for free.
- TOGA's 0.5-threshold-on-a-human–rat-model approach to labels (see 4).

## Suggested order

2 (M state) → 3 (indel masking) → 6 (compound loci) → 1 (neighbour-anchor score) → 5 (weighted graph)
→ 4 (probability) → 7 → 8. Items 2, 3 and 6 are local changes to `ComPair.sh`/the clustering step with
measured gains; 1 and 5 need the registry; 4 needs the label sets assembled first.
