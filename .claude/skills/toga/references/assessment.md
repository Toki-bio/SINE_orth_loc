# TOGA / TOGA2 — assessment

Based on reading the source (TOGA @861a837, TOGA2 v2.0.10) and wikis; paper-level claims from abstracts.

## Strengths

1. **Orthology from genomic context, not similarity.** The classifier asks "does this chain align the gene's
   introns, flanks and neighbouring genes?" — the property that separates the orthologous locus from
   paralogs, retrocopies and processed pseudogenes. This is why TOGA beats similarity-based orthology on
   recently duplicated and retrocopied genes, where sequence similarity is nearly uninformative.
2. **Annotation and orthology in one step, with explicit uncertainty.** Every query gene comes with its
   reference ortholog, an orthology probability, a loss status, and the exact mutations. Missing (assembly)
   is separated from Lost (biology) — a distinction most annotation pipelines don't make.
3. **Codon/splice-aware alignment.** CESAR2 aligns reading frames and splice sites jointly, so exon
   boundaries, frameshifts, stops and splice-site mutations are called consistently; mutation masking rules
   (terminal 10%, compensation, U12, Sec stops) encode years of curation experience.
4. **Scale and homogeneity.** Same method, same reference, same rules across hundreds-thousands of genomes
   (TOGA1: 488 placental mammals + 308 birds; TOGA2: 2,162 vertebrate assemblies). For comparative work the
   consistency matters more than per-genome optimality. Paper: TOGA1 BUSCO completeness higher than Ensembl in
   97% and NCBI in 91.5% of species; TOGA2 513× less memory, 6.1× faster (exon-wise CESAR).
5. **Assembly QC for free.** % intact ancestral genes is a gene-level assembly-quality measure that correlates
   with BUSCO but uses far more genes.
6. **TOGA2 additions are the right ones:** exon-wise alignment (memory), SpliceAI to tolerate intron gain/loss and
   splice-site shifts (fewer false losses), FI status, UTRs, gene-tree step for many:many, multi-reference
   integration, fragmented-assembly stitching by default, orthogroups for CAFE5 with species QC, Rust speed-ups.
7. **Inspectability.** UCSC tracks with per-projection mutation plots, exon meta tables, rejection reasons.

## Weaknesses and failure modes

1. **Projection-only.** Genes absent from the reference annotation (lineage-specific genes, novel families,
   de novo genes) cannot be found; quality is capped by the reference annotation. Multi-reference `integrate`
   only partly mitigates this.
2. **Everything rides on the chains.** Alignment sensitivity and chaining define which loci are even
   considered. Chains from other pipelines (or from HAL) are under-documented ("To be expanded"); chain score
   cut-offs (5,000 / 15,000) are vertebrate-tuned. Fragmented or low-scoring chains silently lose orthologs.
   Chaining artefacts (read-through chimeras) need heuristic splitting (5×, 500 kb rules).
3. **The orthology classifier is small, old and trained on one pair.** XGBoost, 50 depth-3 trees on 4–5
   features, trained once on human–rat; the same `train.tsv` ships unchanged in TOGA2. Labels are weak:
   positives are single-chain 1:1 genes (easy cases), negatives are "non-top chains over 1:1 genes", not
   curated paralogs. Probabilities are therefore not calibrated for other clades or distances, and the
   0.5 threshold is a convention. The long-distance model is optional and was labelled experimental.
   Features are alignment-coverage fractions whose distributions shift with divergence and repeat content
   (e.g. TE-rich introns/flanks in squamates align less → lower `flank_cov`/`intr_perc` for true orthologs).
4. **Gene-loss calls are rule-based and assembly-sensitive.** Fixed thresholds (2 exons / 20% exons / 40%
   CDS exon / middle 80%) with no probability; a single assembly base error can produce UL, two can produce
   L. Short-read or old assemblies inflate losses; "Missing" depends on detecting N-gaps. UL is counted as
   functional by default, which hides early pseudogenisation but also hides real losses. Status is the best
   orthologous projection — in a 1:many case one intact copy makes the gene "intact" even if the syntenic
   copy is lost (paralogous projections only count when no orthologous chain exists).
5. **SpliceAI is human-trained**; the paper reports good transfer across vertebrates, but thresholds
   (0.02 probability, 600 bp deviation — the latter marked TODO in code) are heuristics, and splice-site
   "correction" can manufacture intactness. The correction-mode default differs between `run` (3) and
   config/defaults (0) — runs launched differently are not equivalent.
6. **Gene-tree step is thin.** Two-species trees (reference + query copies only), midpoint-rooted, no
   species tree or outgroup, PRANK + IQ-TREE on cliques ≤ 50; duplication/speciation labelling uses a lookup
   table the authors themselves question (S+S → S). Fine for clean 2:2 cases; not a real reconciliation for
   large tandem families (ZNFs, olfactory receptors, KRABs), which are exactly where many:many arises.
7. **Pairwise, reference-centric orthology.** Orthology is reference→query; multi-species orthogroups are
   built afterwards by union-find over shared query orthologs, which inflates families when one species has
   chimeric or spanning genes (hence the OrthoZ/FamZ QC).
8. **Software maturity.** Large, fast-moving codebase (~50 kLOC, Python/Rust/Cython/C/Nextflow/HDF5),
   heavy dependency stack, no unit tests (CI = build + one end-to-end test), wiki pipeline pages still stubs,
   behaviour changing between patch releases (2.0.5→2.0.10 changed naming, fragment stitching, paralog
   handling). Results must be tied to an exact version.
9. **Coding genes only.** No ncRNA, no non-coding conserved elements (UTRs are projected, not modelled).

## When to use / not

Use: many genomes of a clade with one excellent reference annotation; chromosome-level or good scaffold
assemblies; questions about orthologs, gene loss, duplication, selection (codon alignments), assembly QC.
Prefer de novo/evidence-based (BRAKER/EGAPx/Helixer) or combine (TOGA as hints/CAT-like) when: lineage-specific
genes matter, reference is distant (rule of thumb, not from the papers: beyond ~100–150 My in vertebrates, or non-model clades), or you need ncRNA.
For gene-family evolution, combine TOGA2 orthogroups with a species-tree-aware reconciliation if families are large.

## Checks to run on any TOGA result

1. Version, chain source, chain score distribution, reference annotation version recorded.
2. % FI+I of reference genes vs a closely related well-assembled species (outliers = assembly or alignment problem).
3. Orthology probability distribution — bimodal near 0/1 is good; mass near 0.5 means the classifier is out of its training domain.
4. Gene losses: inspect a sample in the browser (mutation plots); check reads support frameshifts/stops; check shared mutations across related species/assemblies.
5. many2many and `#paralog` counts; fragmented projections (`$` names) concentrated on few contigs → assembly breaks.
6. For orthogroups: OrthoZ/FamZ flags before CAFE5.

## Relevance for squamate / Darevskia work

- Reference: a lacertid with good annotation (e.g. *Podarcis*, *Lacerta agilis*) rather than human/chicken;
  divergence within Lacertidae is close to the training regime; human/chicken → lizard is the far end.
- Expect lower `flank_cov`/`intr_perc` in TE-rich intergenic regions; check the probability distribution and
  consider retraining with lizard 1:1 orthologs (`create_train.ipynb` logic) if it smears.
- Gene-loss calls in the parthenogen haplotypes must be cross-checked with the parental species and reads:
  haplotype assemblies of hybrids are exactly where assembly errors look like losses.
- TOGA's principle (orthology from flanking context and synteny, not similarity) is the same one SINE_orth_loc
  uses for SINE loci; TOGA's edge/graph resolution (many2many splitting) and the pan-SINEome anchor graph
  (`sine_registry.py` edges/breakpoints) are complementary — gene and TE anchors on the same genomes.
