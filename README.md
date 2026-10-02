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

## Alignment bundles

Clusters are processed in batches of 100 in the scratch directory, and the resulting
alignments are appended to one bundle per category instead of one file per locus:

    aln_<sp1>-<sp2>_PM.aln.gz  aln_<sp1>-<sp2>_MP.aln.gz  aln_<sp1>-<sp2>_SINE.aln.gz
    aln_<sp1>-<sp2>_rejected.aln.gz   (alignments that failed ComPair.sh QC)

A bundle is plain text: a `##FILE <name>` line followed by that alignment in FASTA,
repeated. Names are the same as the former individual files (e.g. `C12R.cl.dbl.PM`).

- `sine_loci_browser.html` → **Open Bundle** opens one or more bundles directly
  (**Open Folder** still works for loose files).
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

States: **P** SINE present, **A** empty site with orthologous flanks, **U** no data,
**X** ambiguous (contradicting calls, flag `inconsistent:<sp>`, or several loci of one
species joined into one group, flag `multicopy:<sp>`). To add a genome, run its
comparisons with the existing ones and rebuild with `--previous <last build>`: groups keep
their IDs. The format is described in [`docs/pan-sineome.md`](docs/pan-sineome.md).

For runs made before `orth_*.tsv` existed, rebuild the table from the stat files and the
alignments (a directory of loose files or the bundles):

    sine_registry.py orth --stat stat_doubles_aaa-bbb stat_multi_aaa-bbb --aln work_dir/ \
        --species aaa=aaa.bnk.fai bbb=bbb.bnk.fai -o orth_aaa-bbb.tsv
