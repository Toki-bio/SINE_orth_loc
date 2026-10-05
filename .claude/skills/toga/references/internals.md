# TOGA / TOGA2 internals (verified in source)

Versions: TOGA @861a837, TOGA2 v2.0.10 @c4e1150. Paths are repo-relative.

## 1. Orthology features (per reference transcript × chain)

TOGA1 `chain_runner.py:get_features`, TOGA2 `src/rust/src/bin/feature_extraction.rs` (same definitions,
TOGA2 adds two). Computed with an overlapSelect-style intersection of chain aligned blocks (reference
side) with tracks built from the BED12: exons, CDS, gene body ("grange"), gene ± 10 kb flanks.

| feature | definition | intuition |
|---|---|---|
| `synt` (`synt_log` = log10) | number of reference genes whose CDS the chain covers by ≥ 1 base | orthologous chains span neighbours |
| `gl_exo` | CDS bases aligned by the chain / all aligned bases of the chain | paralog/retrocopy chains are exon-dominated |
| `loc_exo` | CDS aligned / (gene body aligned − UTR exons aligned) | ~1 when introns don't align (processed pseudogene) |
| `flank_cov` | aligned bases in the two 10 kb flanks / 20,000 | orthologous context aligns |
| `intr_perc` | intron bases aligned / intron length | |
| `exon_perc` | exon bases aligned / exon length | |
| `exon_qlen` | CDS aligned / (chain query span − UTR aligned); 0 if ≤ 1 CDS exon covered | short chain made of exons only |
| `gl_score`, `chain_len` | chain score, reference span | |
| TOGA2 `clipped_exon_qlen`, `clipped_intr_cover` | same ideas restricted to the chain clipped to the CDS | processed-pseudogene detection |

Models (`classify_chains.py`): ME features `gl_exo, loc_exo, flank_cov, synt_log, intr_perc`;
SE features `gl_exo, flank_cov, exon_perc, synt_log`. XGBoost, 50 trees, depth 3, lr 0.1
(TOGA2 `train_model.py` hardcodes metaparameters; LD model depth 4).

Training set (`models/train.tsv`, 40,440 rows, human–rat; **byte-identical in TOGA1 and TOGA2**),
built by `models/create_train.ipynb`:
- longest isoform per gene; only Ensembl `ortholog_one2one` human–rat genes; genes overlapped by ≤ 10 chains; ≥ 10% exon coverage;
- positives (y=1): genes covered by exactly **one** such chain;
- negatives (y=0): for genes covered by several chains, every chain except the lowest chain ID (= highest score in UCSC chain files).
So "paralog" in training means "non-top chain over a 1:1 gene" — partly real paralogs/retrocopies, partly repeats and fragments; positives are the easy single-chain cases.

Long-distance (LD) model (TOGA2 `models/train_long_distance.tsv`, 29,374 rows, 11,359 pos / 18,015 neg;
features = the above + first-model `score` + `single_exon`): stacked on the first model, **off by default**,
TOGA1 logged it as "experimental, not recommended for research purposes yet".

Rule-based classes:
- spanning chain: `exon_cover == 0 & synt > 1` → score −1 (no CDS aligned; used for Missing/Lost of whole genes).
- processed pseudogene (−2): ME, `synt == 1 & exon_qlen > 0.95 & pred < thr & exon_perc > 0.65`;
  TOGA2 additionally: any chain below threshold, or a non-best orthologous chain of a transcript with
  several orthologous chains, with `clipped_exon_qlen > 0.3 & 0 ≤ clipped_intr_cover < 0.1`.
- chains with score < `min_orthologous_chain_score` (15,000) are not used for orthologous annotation (only PP/retro).

(TOGA1 bug, cosmetic: the log counts spanning/PP as `pred == 1.0 / 2.0` instead of −1/−2.)

## 2. Locus definition and alignment

TOGA1 (`CESAR_wrapper.py`): one CESAR2 run per projection over the whole chain-projected locus;
memory ≈ CESAR's own formula on exon lengths × locus length (`memory_check`), jobs bucketed by memory.
Exons inside assembly gaps detected from N runs (`--gap size`); `--fragmented_genome` optional stitching.

TOGA2 (`cesar_preprocess.py`, `modules/cesar_wrapper_executables.py`):
- each exon gets a chain-defined expected locus + search space (flank 50 bp, `EXTRA_FLANK` 10%);
- `group_defined_exons`: exons whose search spaces intersect form one CESAR group; exons absent from the
  chain are attached to neighbouring groups and searched in the inter-exon space (`--max_space_size`);
  spanning-chain projections are one group;
- memory estimated per group (`cesar_memory_check`), limit 15 GB default → paper: 513× less memory, 6.1× faster;
- chimeric ("read-through") projections split when two exons' search spaces are separated by > 5× the
  reference intron and > 500 kb, or an exon projects to > 5× its length; the side with more aligned
  sequence is kept (`MAX_CHAIN_INTRON_LEN` 500,000; `MAX_CHAIN_GAP_SIZE` 1,000,000);
- fragmented projections (exons on different contigs) recovered by `stitch_fragments.py` by default;
  same-contig fragments only with `--enable_same_contigs` (since v2.0.10).

Exon presence (TOGA2): **Present** if it overlaps its chain locus by ≥ 1 bp, OR nucleotide identity ≥ 45%
and BLOSUM score ≥ 20% (`MIN_ID_THRESHOLD`, `MIN_BLOSUM_THRESHOLD`), OR both splice sites have SpliceAI
probability ≥ 0.02 (one site for terminal exons). Otherwise **Missing** if outside the chain span or the
search space contains an assembly gap (or ≥ 90% N), else **Deleted**.

SpliceAI (`spliceai_manager.py`, `predict_with_spliceai.py`, `intron_gain_check.py`): genome-wide donor/
acceptor probabilities for the query (human-trained models). Uses: rescue exon presence; correct exon
boundaries (`--spliceai_correction_mode` 0–7; mode 3 = correct all canonical U2 sites when a better-supported
alternative exists); split a CESAR exon by a new query intron (min intron 10 bp) when that removes
frameshifts/stops ("intron gain", masked); `MAX_DEV_FROM_SPLICEAI = 600` (marked TODO in code).

## 3. Mutations and loss status

TOGA1 `modules/inact_mut_check.py` + `modules/gene_losses_summary.py`; TOGA2 in the alignment step.
Mutation classes: FS_INS/FS_DEL, STOP, SSMD (donor ≠ GT/GC), SSMA (acceptor ≠ AG), BIG_INS/BIG_DEL (≥ 50 bp,
masked), START_MISSING (masked), STOP_MISSING (masked), Missing exon, Deleted exon, COMPENSATION, intron deletion.
Compensation: ≥ 2 frameshifts summing to 0 mod 3 with no stop in the shifted-frame segment.

TOGA1 decision tree (`get_projection_classes`), thresholds: intact if no mutations and %intact(missing ignored) > 0.6;
%intact codons < 0.35 → L (or M if > 65% of frame outside the chain); < 0.49 → UL; "middle 80%" test;
missing < 0.5 → PI else M/PM; %intact(missing as intact) < 0.2 → L; multi-exon L if affected exons ≥
max(2, 20% of exons) or a > 40%-CDS exon deleted or with ≥ 2 mutations; single-exon L if < 60% intact and ≥ 2 mutations.
TOGA1 classes: I, PI, UL, L, M, PM (partially missing), PG, N. TOGA2 drops PM, adds FI and PP, and lets masked
terminal mutations / compensations count once a real inactivating mutation exists. Full TOGA2 criteria: wiki
"Loss status classification" and "Computing CDS integrity statistics".

## 4. Genes and orthology

Query genes (`modules/infer_query_genes.py`): projections on the same strand sharing ≥ 1 CDS base; orthologs
evict overlapping paralog projections; both evict processed pseudogenes.

Initial resolution (TOGA1 `modules/orthology_type_map.py`, TOGA2 `modules/initial_orthology_resolver.py`):
bipartite graph reference genes ↔ query genes (edges from orthologous, non-lost projections, weight = orthology
probability) → connected components → one2one / one2many / many2one / many2many / one2zero.
many2many splitting (`split_graph`): cut leaf-adjacent edges only if no isolated nodes appear, no "strongly
connected" reference genes (> 1 shared neighbour) are separated, and every removed edge < 0.75 or
< 0.9 × the weakest preserved edge; otherwise the clique stays.

Fine resolution (TOGA2 `fine_orthology_resolver.py`, `modules/tree_analysis.py`): for unresolved cliques
≤ 50 sequences, protein sequences → PRANK MSA → IQ-TREE2 (`--alrt 5000 -B 5000`, model set restricted) or RAxML
(100 bootstraps) → `.contree` → **midpoint root** → each internal node labelled from children (R/Q leaves; R+Q = S
speciation, R+R = Dr, Q+Q = Dq, …; lookup table `CAT_DDICT`, which labels S+S as S — the code itself asks
"So ancestor for S-S is also S?", and a 4-leaf `special_case` patches the obvious instance) → leaf pairs under an S
node become 1:1 if no "P" node lies on the path to the root and all supports on that path ≥ 90 (UFBoot).
Only reference and query sequences are in each tree (no outgroup).

Second-best projections: when several reference genes project onto one query locus, lower-probability ones are
discarded as likely wrong orthologs.

UTRs (`src/rust/src/bin/utr_projector.rs`): reference UTR exons projected through the chain, capped by absolute
(3,000 bp) and relative (2.5× reference) length thresholds, optional boundary extrapolation.

## 5. Multi-reference and multi-species

`toga2.py integrate` (`modules/integrate.py`): reads several TOGA2 results for one query (JSON of references with
priority), takes loci from all references, infers genes, picks best isoforms by priority/intactness/chain support,
writes merged BED/GTF and a reference-support matrix.

`toga2orthogroups.py`: union-find over reference genes that share a query ortholog in any species (optional
PANTHER merging) → CAFE5 matrix; `--one-to-one` per-reference-gene single-copy matrix; species QC: OrthoZ (rate
of query genes spanning several reference genes) and FamZ (copy-number outliers), flag at z > 3.

## 6. Engineering

TOGA1: Python + C helpers + Nextflow/para; CESAR2.0 submodule; postoga. TOGA2: Python ≥ 3.10 (README says
3.13), Rust (feature extraction, UTRs, plots, filters), Cython, HDF5 stores, Nextflow, IQ-TREE2/PRANK, SpliceAI;
Apptainer recipes. Testing = CI build + `toga2.py test` end-to-end on bundled hg38/mm10 data; no unit tests.
Changelog shows behaviour changes in every 2.0.x release — record the version with results.
