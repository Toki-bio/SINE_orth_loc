# SINE_orth_loc
Search for the orthologous SINE-containing loci in two genome assemblies

## Versions

| Version | Commit | What it is |
| :- | :- | :- |
| [`v1.0-genes2023`](https://github.com/Toki-bio/SINE_orth_loc/tree/v1.0-genes2023) | `5c6c456` | The workflow used in Kosushkin et al. 2023, *Genes* 14(11):2089 ([doi:10.3390/genes14112089](https://doi.org/10.3390/genes14112089)) |
| [`v2.0`](https://github.com/Toki-bio/SINE_orth_loc/tree/v2.0) | | Current version (this branch) |

To reproduce the published analysis, use `v1.0-genes2023`:
`git clone --branch v1.0-genes2023 https://github.com/Toki-bio/SINE_orth_loc`.

Changes in v2.0, with the same classification logic (`ComPair.sh` thresholds unchanged):

- Alignments are written to one bundle per category instead of one file per locus;
  batches run in scratch space (see below). Added `aln_bundle.sh` and bundle mode in
  `sine_loci_browser.html`.
- `orth_<sp1>-<sp2>.tsv` per pair and `sine_registry.py` for multi-species locus matrices.
- `SINE_orth_loc_flexible.sh` wrapper (input checks, `--threads`, `--scratch`, `--force`).
- Fixes that can change results of some runs:
  - chromosome names containing dots (e.g. `NC_000001.11`) produced empty coordinates in
    `statbed_*`;
  - parallel `ComPair.sh` jobs appended to one shared `stat` file, which can lose or
    interleave lines on network file systems;
  - `bwa mem -t=N` was not read as a thread count (mapping ran single-threaded);
  - the tool check did not detect missing tools, and mawk instead of GNU awk silently
    produced no loci.
- Removed `SINE_orth_loc.sh`, an outdated copy of `SINE_orth_loc.bash` that called
  `ComPair.sh` by its former name `script10` (still available in `v1.0-genes2023`).

On a simulated genome pair, v2.0 gives the same `statbed_*`, `MP_PM_SINE_*` and
alignments as the previous version.

If you use SINE_orth_loc, please cite the paper (see `CITATION.cff`).

## Usage

    SINE_orth_loc_flexible.sh -g1 sp1.fa -g2 sp2.fa -s SINE.fa -b1 sp1_SINE.bed -b2 sp2_SINE.bed \
        -n1 sp1 -n2 sp2 -o results_sp1-sp2 -t 16 --scratch "$TMPDIR"

`SINE_orth_loc_flexible.sh --help` lists all options. Requires mafft, esl-alipid, seqkit,
bedtools, samtools, bwa, sam2bed (BEDOPS) and GNU awk.

On a cluster, set `-t` to the number of allocated cores and `--scratch` to node-local
storage (defaults to `$SLURM_TMPDIR` or `$TMPDIR` when set). Use `--force` to reuse an
existing output directory from a batch job.

## Multi-copy loci (multimappers)

When a flank matches several places, loci join into clusters of 3–10 (multi) or more
(poly). The multi stage keeps at most one pair per cluster; poly clusters were never
analysed. Since v2.1, every copy of a multi-copy cluster that is not yet in a final
alignment goes through a **two-flank rescue** (`rescue_multi.py`): left and right flanks
are mapped with `bwa mem -a`, the best hits of each flank serve as anchors, and the other
flank is aligned in the window next to each anchor where it must lie (same strand, at
the distance of an empty site or of a SINE). The placement whose two flanks score clearly
best (`--margin`, default 20) is accepted and checked with ComPair.sh like a double
(clusters `M<n>R`); ties are reported as ambiguous instead of guessed.

- `clusters_<sp1>-<sp2>.tsv` — every cluster: class, number of loci, fate (resolved and
  how, rejected by ComPair.sh, unresolved) and the rescue outcome of its copies
- `rescue_<sp1>-<sp2>.tsv` — per copy: resolved (and the new cluster), ambiguous, no hit
- `RESCUE=0` (or `--no-rescue` in the wrapper) reproduces the v2.0 behaviour

### Running a pair in pieces (`STAGE`, `SHARD`)

Where jobs have a time limit (CI runners, array jobs), a pair can be split; the alignment
stage, which takes most of the time, runs in independent shards that need only the `.fai`
files, the BED files and the consensus:

    STAGE=prep   SINE_orth_loc.bash a.bnk b.bnk SINE.fa   # mapping, clustering, batch files
    STAGE=align SHARD=3/8 SINE_orth_loc.bash a.bnk b.bnk SINE.fa   # batches 3, 11, 19, ... -> shard_3/
    STAGE=finish SINE_orth_loc.bash a.bnk b.bnk SINE.fa   # merge shard_*/, rescue, tables, nesting

The shards can run on different machines (copy `shard_<k>/` back before `finish`). The result
is identical to `STAGE=all` (default): same orth table, cluster table and alignment bundles.

On a simulated 4-species set where a quarter of the loci have a young repeat (~1%
divergence) as left flank, plus 10 lineage-specific segmental duplications
(`tests/`), recall per pair rose from 207–283 to 257–324 of 258–324 loci (repeat-flank
loci: from 17–40 to 63–81 of 64–81), with no wrong pairs; the duplicated loci are
reported as ambiguous and come out as `multicopy` in the registry.

The multi stage used to tell the two species apart by the first 3 characters of the
sequence names, which works only for genomes renamed with species prefixes (`dva_chr1`).
With GenBank-style names (`CM0…`, `JAW…`) every candidate pair looked like one species
and all multi clusters were dropped silently (empty `stat_multi_*`). It now uses the
chromosome list of genome 1.

## Nested, close, satellite and edge cases (v2.2)

**SINE-in-SINE.** A SINE inserted into an older SINE splits it in two; the annotation (sear2k:
≥ 80% length) reports the young copy but not the host halves, so the young copy's flanks *are*
SINE sequence, and in a window holding host and insert the consensus aligns to whichever is most
similar — a flank-based pair can join an old host in one genome with a young insert in the other.
`sine_nest.py` (`NEST=1`, default; needs `nhmmer` from HMMER):

- `scan` searches every copy ± 400 bp with nhmmer for SINE pieces ≥ 20 nt (TinT's minimum) and
  classifies each locus: **nested** (two host pieces whose consensus coordinates join up around an
  insert, the target site duplication verified in the sequence — reassembled hosts in
  `nest_<sp>.reassembled.bed`), **split** (one copy annotated as two pieces), **dimer**, **satellite**
  (≥ 3 periodic units sharing their spacers: not independent insertions), **close**, single;
- `orth` compares compound loci through the flanks *outside* the whole locus (unique sequence) and
  tests each element by its junctions against its empty-site junction with one TSD removed — the
  host as one element from its outer ends, the insert inside it;
- `supersede` moves flank-based rows with an anchor inside a compound locus to
  `orth_<a>-<b>.superseded.tsv`; the compound rows (`N<k>R`) replace them;
- `tint` orders (sub)families in time from nesting events (who inserted into whom), after
  Kriegs et al. 2007 / Churakov et al. 2010 (TinT): one normally distributed activity period per type,
  no target preference ⇒ P(A inserted into B | A–B nesting) = Φ((μA − μB)/√2); μ by maximum likelihood,
  bootstrap intervals, types compared one way only reported as bounds. This is a reconstruction from the
  published model assumptions, not the original TinT code. Within a nesting event the host is the
  older element by definition; host and insert identity to the consensus are reported as a check.

**Flank indels.** A second insertion or deletion of ≥ 50 bp next to the site used to fail the right
flank (`badRF`); ComPair.sh now re-tests the right flank with such one-sided gap runs masked
(`CP_INDEL_MIN`, 0 = off) and marks passing calls `FI=1`.

**Contig ends.** Loci whose window is clipped by a sequence end are `contig_end` (the sequence is
missing, not diverged) and enter the orth table as `MISSING` rows; the registry state is **M**.

**Exact anchors.** ComPair.sh reports where the SINE starts in each aligned locus; the orth table
carries `anchor1/anchor2`, and the registry joins sites with these exact anchors within
`--precise-tol` (20 bp) instead of estimating them from window lengths (60 bp), which kept a nested
insert's empty site apart from its host's junction.

Tests (`tests/simulate_genomes2.py`: nested SINEs with TSDs, close insertions, flank indels, satellite
arrays, contig ends, a haplotig, subfamilies; `tests/score_pairs2.py` scores with the alignment
anchors against the full element truth):

| class (4 pairs) | flank-based only (`NEST=0`) | with `sine_nest.py` (`NEST=1`) |
| :- | :- | :- |
| nested inserts | 18/19/24/27 of 37/54/56/51, 2–4 wrong per pair | 34/54/56/48, 0 wrong |
| hosts | 69/65/66/64 of 88, 2–6 wrong | 79/76/76/72, 0 wrong |
| close copies | 34/31/31/25 of 59/75/71/73 | 54/63/58/55 |
| plain / repeat flank | 100%, 0 wrong | 100%, 0 wrong |

Flank indels (89 per pair): 74/41/43/63 before masking, 80/54/56/72 with it, 0 wrong either way.

`scan` finds 106/110 nested inserts (0 false, TSD in all, host the more diverged in all) and 14/14
satellites; `tint` recovers the true order of 6 types (Spearman 1.00) from simulated insertion
histories, also with unequal activity widths and ~180 events. Registry over the 4 simulated genomes:
informative groups 38 → 130, groups joining different events (`family_mixed`) 37 → 19. On the real
dva–mix rejected alignments, 198 of 3,000 `badRF` loci (6.6%) pass after indel masking.

Real genomes (D. valentini dva, D. mixta mix; `darevskia-sineome-data`, branch `actions/nest`): 47 and 32
nestings (host the more diverged element in 41/47 and 30/32, TSD found in 36/47 and 25/32, median 14 bp),
10 and 7 satellite arrays, 481 and 373 dimers. Of the dva–mix flank-based rows, 736 fall in compound loci;
the compound test agrees with 702 of them, differs in 18 (15 times flank "SINE" vs compound PM/MP — the
host/insert pattern) and has no call for 16 nested ones. It adds 3,630 calls at elements without any
flank-based row (2,661 SINE, 969 PM/MP), mostly copies that were excluded as closer than 300 bp.

Real-data precision of ComPair.sh calls, measured with reads (`tests/build_training_set.py`): in the six
pairs with D. valentini dvl, the dvl-side call of 265,145 rows was checked against read genotypes of the
assembled individual (GQ ≥ 20, homozygous, unflagged groups): 0.22% of SINE, 0.77% of PM and 0.45% of MP
calls are contradicted — about the genotyper's own error rate (0.37%). A logistic model on the ComPair
metrics (`tests/fit_error_model.py`, held out by scaffold) ranks the contradicted calls with AUC 0.70:
the 1% highest-risk calls are contradicted 5.2% of the time (18× the average), so the score is useful to
pick calls for inspection, not to replace the thresholds.

## Alignment bundles

Clusters are processed in batches of 100 in the scratch directory, and the resulting
alignments are appended to one bundle per category instead of one file per locus:

    aln_<sp1>-<sp2>_PM.aln.gz  aln_<sp1>-<sp2>_MP.aln.gz  aln_<sp1>-<sp2>_SINE.aln.gz
    aln_<sp1>-<sp2>_rejected.aln.gz   (alignments that failed ComPair.sh QC)

A bundle is plain text: a `##FILE <name>` line followed by that alignment in FASTA,
repeated. Names are the same as the former individual files (e.g. `C12R.cl.dbl.PM`).

- `sine_loci_browser.html` → **Open Bundle** opens one or more bundles directly
  (**Open Folder** still works for loose files).
- **⚙ Thresholds** in the browser re-classifies the loaded alignments with the ComPair.sh
  rules (an exact JavaScript port) under thresholds you set: counts per class before and
  after, the changed alignments marked (e.g. `PM→badRF`, filter **changed**), and for each
  alignment its metrics and the test that decides it. Load the `rejected` bundle as well to
  see what relaxed thresholds would accept. **Copy as env** gives the setting for the next run:
  `CP_FLANK_LEN=150 CP_FLANK_NID=65 CP_FLANK_PID=65 CP_SINE_PID=65 CP_SINE_NID=100
  CP_CUT_PID=70 CP_CUT_NID=120 SINE_orth_loc_flexible.sh ...` (these are the defaults).
- `aln_bundle.sh list|count|get|extract` lists, prints or restores individual files, e.g.
  `aln_bundle.sh extract -p '^C12R\.' out/ aln_sp1-sp2_*.aln.gz`.
- `aln_bundle.sh pack [--remove] DIR PREFIX` bundles the loose `.PM/.MP/.SINE` files of
  runs made with earlier versions.

## Multi-species comparison (pan-SINEome)

Each pairwise run also writes `orth_<sp1>-<sp2>.tsv`: one row per validated alignment with
both loci, their species and whether each carries the SINE (plus the ComPair.sh metrics).
`sine_registry.py` (Python 3, standard library only) combines the pairs into orthologous
groups — one per ancestral insertion, with the copy or the empty site in every genome:

    sine_registry.py build results_*/results/orth_*.tsv --species lag,dva,pmu,pra \
        --copies lag=lag-Squam1.bed --subfamilies lag=lag/assignment_full.tsv \
        --copies dva=dva-Squam1.bed --subfamilies dva=dva/assignment_full.tsv ... -o lacertids

`--copies` takes the [sear2k](https://github.com/Toki-bio/sear2k) BED of each species and
family (`<sp>-<FAMILY>.bed`, family from the file name); `--subfamilies` takes the
[SINEderella](https://github.com/Toki-bio/SINEderella) step-2 `assignment_full.tsv`.

It matches the same locus of a species across comparisons (same chromosome and strand,
insertion site within `--tol` bp, default 60), joins loci through the orthologous pairs and
writes

- `lacertids.groups.tsv` — one row per group: family, subfamily, state pattern, flags, and per species
  the copy (`chrom:start-end(strand)`) or the empty site (`chrom:pos(strand)`)
- `lacertids.copies.tsv` — with `--copies`: every annotated copy with identity, bitscore and
  subfamily, its group, or why it has
  none (`close_copy`, `no_validated_pair`, `not_compared`)
- `lacertids.evidence.tsv` — the pairwise rows behind each group
- `lacertids.matrix.tsv`, `.patterns.tsv`, `.nex` — states, pattern counts, NEXUS 0/1 matrix
- `lacertids.aliases.tsv` — IDs of the previous build that were merged, split or retired
- `lacertids.edges.tsv`, `.breakpoints.tsv`, `.dupblocks.tsv` — the anchor graph: groups
  that are neighbours along each genome with their spacers, adjacencies broken between genomes
  (`a_specific`: kept by no other genome — a misjoin or a lineage rearrangement), and runs of
  multicopy sites (duplicated or twice-assembled regions)

States: **P** SINE present, **A** empty site with orthologous flanks, **U** no data,
**X** ambiguous (contradicting calls, flag `inconsistent:<sp>`, or several loci of one
species joined into one group, flag `multicopy:<sp>`). To add a genome, run its
comparisons with the existing ones and rebuild with `--previous <last build>`: groups keep
their IDs. The format is described in [`docs/pan-sineome.md`](docs/pan-sineome.md).

## Genotyping pan-SINEome loci from reads (prototype)

`sine_genotype.py` genotypes the orthologous groups of a pan-SINEome in any sequenced
individual, without assembling it: two templates per group (the allele with the SINE and
the empty allele, from the assembled genomes), reads mapped to them, and per read the allele
it shows (a junction crossed with ≥30 bp on both sides, or the SINE missing from / inserted
into a template); repetitive junctions are not used. Output: 0/0, 0/1, 1/1 with read counts
and a genotype quality (GQ).

    sine_genotype.py templates --groups lacertids.groups.tsv --genome lag=lag.fa dva=dva.fa ... -o tpl
    bwa index tpl.fa && bwa mem -t 16 tpl.fa sample.fq > sample.sam          # HiFi: minimap2 -ax map-hifi
    sine_genotype.py call --templates tpl.tsv --sam sample.sam -o sample.gt.tsv

Simulation (`tests/`: pan-SINEome of 3 species; individuals with 0.5% private mutations),
% correct of called loci (calls / 389 groups):

| individual | short 5× | short 10× | short 20× | HiFi 10× |
| :- | :- | :- | :- | :- |
| new individual of a species in the pan-SINEome | 100 (336) | 100 (380) | 100 (389) | 100 (386) |
| individual of a species **not** in the pan-SINEome | 100 (323) | 100 (381) | 100 (386) | 99.7 (383) |
| hybrid of two species (heterozygous loci) | 83 (353) | 93 (388) | 98 (387) | 97 (387) |
| hybrid, only GQ ≥ 20 | 100 (51) | 100 (231) | 100 (351) | 99.7 (314) |

Heterozygous loci need depth: at 10× each haplotype gives ~3 junction reads and the empty
allele is sometimes missed by chance; such calls get a low GQ. For hybrid parthenogens use
≥20× short reads or HiFi, and filter on GQ.

For runs made before `orth_*.tsv` existed, rebuild the table from the stat files and the
alignments (a directory of loose files or the bundles):

    sine_registry.py orth --stat stat_doubles_aaa-bbb stat_multi_aaa-bbb --aln work_dir/ \
        --species aaa=aaa.bnk.fai bbb=bbb.bnk.fai -o orth_aaa-bbb.tsv
