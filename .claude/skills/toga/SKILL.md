---
name: toga
description: Expert knowledge of TOGA and TOGA2 (Hiller lab, Tool to infer Orthologs from Genome Alignments) — how they annotate genes, classify orthologs vs paralogs vs processed pseudogenes, call gene loss, and resolve 1:1/1:many/many:many orthology from pairwise whole-genome alignment chains. Use when planning, running, debugging or interpreting a TOGA/TOGA2 run, choosing a reference/query setup, reading its output files or loss statuses, judging whether its results can be trusted for a species pair, or comparing it with other annotation/orthology tools (BUSCO, Liftoff, OrthoFinder, CAT, BRAKER).
---

# TOGA / TOGA2

Source-level knowledge, verified against code: TOGA `hillerlab/TOGA` @861a837 (final commit, Nov 2025:
"switch to TOGA2") and TOGA2 `hillerlab/TOGA2` v2.0.10 @c4e1150 (Sep 2026), plus both GitHub wikis.
Papers: Kirilenko et al. 2023 *Science* 380:eabn3107 (TOGA1); Malovichko et al. 2026 bioRxiv
10.64898/2026.06.30.735536 (TOGA2). The paper full texts were not readable when this skill was
written; paper numbers below are from abstracts/search snippets — say so when quoting them.

Deeper material:
- `references/internals.md` — step-by-step pipeline, features, thresholds, decision trees, file:function pointers.
- `references/assessment.md` — strengths, weaknesses, failure modes, when (not) to use, checks to run.

## The idea in one paragraph

Instead of searching for genes in a new genome, TOGA *projects* a well-annotated reference
transcript through a pairwise genome-alignment **chain** (UCSC chain format, ideally from the Hiller-lab
`make_lastz_chains` pipeline) into the query. Each transcript×chain pair is a candidate
**projection**. Orthology is decided from **alignment context, not sequence similarity**: an
orthologous chain aligns the gene's introns and 10 kb flanks and spans neighbouring genes (synteny);
a paralogous chain aligns mostly exons; a processed pseudogene aligns exons with the introns missing.
A small XGBoost model turns these features into an orthology probability. The query locus is then
annotated exon-by-exon with **CESAR2** (codon- and splice-site-aware HMM aligner), inactivating
mutations are called, each projection gets a **loss status**, and a gene graph (+ in TOGA2, gene
trees) yields 1:1 / 1:many / many:1 / many:many orthology.

## Pipeline (TOGA2 step names)

1. `setup`/`prepare_data` — filter chains (score ≥ 5,000 kept; ≥ 15,000 required for orthologous use), index chains/BED.
2. `feature_extraction` (Rust) — per transcript×chain: `gl_exo`, `loc_exo`, `flank_cov`, `synt`, `intr_perc`, `exon_perc`, `exon_qlen`, + TOGA2 `clipped_exon_qlen`, `clipped_intr_cover`.
3. `classification` — XGBoost SE (single-exon) / ME (multi-exon) models → probability; ≥ 0.5 ORTH, else PARA; spanning chains −1; processed pseudogenes −2 (rules, not ML).
4. fragmented-projection recovery (default ON in TOGA2) — stitch exons of one transcript across chains on *different* contigs.
5. `preprocessing` — exon search spaces from the chain; exons grouped by overlapping search space (exon-wise CESAR jobs, memory estimate, limit 15 GB default).
6. `alignment` — CESAR2 per exon group; SpliceAI-guided splice-site correction and intron-gain search; exon presence (Present/Missing/Deleted); mutation calling; projection loss status.
7. `gene_inference` — query genes = same-strand projections sharing ≥ 1 coding base; precedence ortholog > paralog > processed pseudogene.
8. `loss_summary` — transcript/gene status = best projection status: `FI > I > PI > UL > L > M > (PG, PP) > N`.
9. `orthology_resolution` — bipartite gene graph → components; many:many split by edge-weight rules; remaining cliques ≤ 50 seqs → PRANK + IQ-TREE2 (5,000 UFBoot), midpoint root, 1:1 pairs at UFBoot ≥ 90.
10. `finalize` — UTR projection, renaming (`_a/_b` copies, `lost_`/`missing_`/`paralog_`/`retro_` prefixes), GTF, UCSC BigBed with mutation SVGs.

## Loss statuses (TOGA2)

| status | meaning (projection level) |
|---|---|
| FI Fully Intact | no inactivating *or masked* mutations, all exons present (new in TOGA2) |
| I Intact | middle 80% of CDS present, no inactivating mutation there |
| PI Partially Intact | ≥ 50% CDS present, middle 80% unmutated, some missing sequence |
| UL Uncertain Loss | ≥ 1 inactivating mutation in middle 80%, below Lost criteria |
| L Lost | ≥ 2 mutated exons (≤ 10 exons) or ≥ 20% mutated exons (> 10); or a ≥ 40%-CDS exon deleted / with ≥ 2 mutations; or ≥ 65% deleted; or longest intact stretch < 20% |
| M Missing | < 50% of CDS present, no evidence of loss (gaps, contig ends, no chain) |
| PG / PP / N | paralogous-only, processed pseudogene, technical reject |

Masked (non-inactivating) by default: first 10% of CDS *if* an in-frame ATG follows within it (all of it with `-m10m`), last 10%, compensated frameshifts (no stop in shifted frame), frame-preserving exon deletions, stops aligned to reference in-frame stops (Sec/readthrough), U12 and non-canonical-U2 splice sites, precise intron deletions, SpliceAI-supported new introns/site shifts, big (≥ 50 bp) in-frame indels, missing start/stop.
TOGA2 change: when ≥ 1 real inactivating mutation exists, masked terminal mutations and compensating frameshifts *count* toward Lost.
Default "functional" set: `--accepted_loss_symbols FI,I,PI,UL` — note UL counts as functional.

## Running — what matters in practice

- Inputs: reference + query `.2bit`; chain file (blank-line separated, no `#` lines); reference BED12 (CDS length % 3 == 0, names `[A-Za-z0-9.-#]` only, recommended `transcriptID#geneID`); isoform file (gene\ttranscript, **no header**; mandatory unless `--no_isoform_file`); U12/non-canonical intron file; optional SpliceAI output for the query (`toga2.py spliceai`).
- `toga2.py prepare-input` validates/prepares the reference; `toga2.py test` is the install check; `cookbook` lists example commands; `from-config` runs from a TSV of arguments.
- Parallelism via Nextflow configs (`nextflow_configs/`), Apptainer image available.
- Multi-reference: run per reference, then `toga2.py integrate` (priority-ordered references → one query annotation). Cross-species orthogroups for CAFE5: `toga2orthogroups.py` (union-find over shared query orthologs, OrthoZ/FamZ species QC).
- Defaults worth knowing: orthology threshold 0.5; long-distance model OFF (`--use_long_distance_model`); `run` CLI default `--spliceai_correction_mode 3` but `defaults.py`/config default 0 (help text says 0) — set it explicitly; tree step on unless `--skip_gene_trees`; max clique 50.
- Pin the version: v2.0.6→2.0.10 each changed gene inference/orthology behaviour (e.g. v2.0.10 changed fragment stitching on same contigs and fixed paralogs being split from orthologs).

## Reading results

- Projection names `transcript#gene#chain[#paralog|#retro]`, fragments `…#chain1,chain2$k`.
- Query gene names follow orthology: 1:1 = reference gene id; 1:many `_a,_b`; many:1 comma list (or `id+` if > 3); prefixes show the *query gene's* best status.
- Key files: `query_annotation.bed` (+ `.with_utrs.bed`), `query_genes.tsv/bed`, `orthology_classification.tsv`, `loss_summary.tsv`, `inactivating_mutations.tsv`, `processed_pseudogenes.bed`, `meta/exon_meta.tsv.gz`, rejection logs.
- A transcript is Lost only if **all** its orthologous projections are lost; retrogenes never rescue a loss.
- Spanning chain only (no CDS aligned): shortest chain gap wins; assembly gap inside → Missing, else Lost.

## Judgement rules (see references/assessment.md)

- Trust TOGA most for close-to-moderate divergence with chromosome-level assemblies and a high-quality reference (human/mouse-like annotation); it is a projection method — it cannot find genes absent from the reference.
- Treat single-assembly gene-loss calls as hypotheses: confirm with reads (frameshift/stop supported?), a second assembly, or related species sharing the mutation.
- The orthology classifier was trained once on human–rat (identical `train.tsv` in TOGA1 and TOGA2); retrain (`train_model.py` + `create_train.ipynb`) or validate on a known pair for very distant pairs or non-mammals.
- Many:many and young tandem families are where both chaining and the classifier are weakest; check `#paralog`, many2many rows and fragmented projections by eye in the browser.
