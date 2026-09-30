# SINE_orth_loc
Search for the orthologous SINE-containing loci in two genome assemblies

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

## Multi-species comparison

Each pairwise run also writes `orth_<sp1>-<sp2>.tsv`: one row per validated alignment with
both loci, their species and whether each carries the SINE (plus the ComPair.sh metrics).
`sine_registry.py` (Python 3, standard library only) combines the pairs:

    sine_registry.py build results_*/results/orth_*.tsv --species aaa,bbb,ccc,ddd -o lacertids

It matches the same locus of a species across comparisons (same chromosome and strand,
insertion site within `--tol` bp, default 60), joins loci through the orthologous pairs and
writes

- `lacertids.loci.tsv` — one row per multi-species locus: state pattern, flags, coordinates
- `lacertids.matrix.tsv` — locus × species states
- `lacertids.patterns.tsv` — number of loci per pattern (e.g. `PPAA`)
- `lacertids.nex` — NEXUS 0/1 matrix of variable loci for PAUP*/MrBayes

States: **P** SINE present, **A** empty site with orthologous flanks, **U** no data,
**X** ambiguous (contradicting calls, flag `inconsistent:<sp>`, or several loci of one
species joined into one locus, flag `multicopy:<sp>`). A new genome only needs its
comparisons with the existing ones; rebuild the registry with the added tables.

For runs made before `orth_*.tsv` existed, rebuild the table from the stat files and the
alignments (a directory of loose files or the bundles):

    sine_registry.py orth --stat stat_doubles_aaa-bbb stat_multi_aaa-bbb --aln work_dir/ \
        --species aaa=aaa.bnk.fai bbb=bbb.bnk.fai -o orth_aaa-bbb.tsv
