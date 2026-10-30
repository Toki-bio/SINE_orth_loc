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

## Bank mode: SINE loci from a SINEderella run

For genomes with many subfamilies and composite SINEs (e.g. *Rhinolophus*, 21 consensi of
105-608 bp), the SINE loci are taken from SINEderella's per-copy list, not searched again, and
one run covers the whole bank:

    sinederella_bed.py RUN_A -o aaa-BANK.bed     # RUN_*: SINEderella result directories
    sinederella_bed.py RUN_B -o bbb-BANK.bed     # (loci_verdict.tsv, family_view/)
    SINE_orth_loc.bash aaa.bnk bbb.bnk BANK.fa   # BANK.fa: the bank both runs used

Bank mode starts when the consensus file has more than one sequence. It changes three things:

- **Copies.** `sinederella_bed.py` keeps copies with verdict `SINE` (assigned to a subfamily, flanks
  unique) and partial copies (`family_view/partial_copies.tsv`, status `partial`; label
  `<consensus>:<a>-<b>`, the span they cover). Grey and other not-SINE copies are left out.
  Overlapping copies become one locus.
- **Reference row.** Every alignment gets the consensus of its own copy (or the covered span of a
  partial copy) as reference row, named `>BANK <subfamily>` (or `>BANK <subfamily>:<a>-<b>` for a
  span), so the subfamily is on the row in every alignment and viewer. A copy at least 10 % longer
  than its labelled consensus gets the longer bank consensus that contains it and comes closest to
  its length (a composite copy whose votes went to a unit inside it). When the two loci carry
  different labels, the longer reference is used. The `reference` column of `orth_<sp1>-<sp2>.tsv`
  names the reference row of every locus (single mode: the consensus name; nested rows: the
  element); `refs_<sp1>-<sp2>.tsv` adds the labels of each cluster's copies. Windows are sized for
  the longest consensus.
- **Insertion junction (`anchor1`, `anchor2`).** The junction of a locus whose annotated copy sits at
  its window's junction is that copy's start (`<cluster>.off`, written with the reference); the other
  locus gets the alignment column where the first reaches it. Without an annotated copy (single
  mode, rescue windows) the junction is where the reference row's first residue aligns, measured
  after ComPair's short-left-flank fix. A reference longer than the copy (partial copy, composite)
  starts anywhere in the left flank; this put ~5 % of rsi–rfe anchors >50 bp off before 2026-10-29.
- **Summary.** `orth_summary.py` (run at the end) writes `summary_<sp1>-<sp2>.tsv`: per reference,
  input copies per genome, loci by row type and status, shared %, and the median (10–90 %) identity
  of copies to the reference, between the two copies of shared loci, and of the left flanks.
- **Bundles per reference.** `orth_split.py` (run after the summary) splits the SINE, PM and MP
  bundles into `by_reference/<reference>/aln_<sp1>-<sp2>_<status>.aln.gz` (references with fewer
  than 20 loci under `other`; `by_reference/index.tsv` lists them), so one subfamily's alignments
  of one status open alone in `sine_loci_browser.html`.
  The browser lists the subfamilies it finds in the loaded alignments (reference row) as a
  drop-down next to the SINE / PM / MP / rejected chips; `?ref=SUBFAMILY` and `?type=PM` in its
  address preset them, as `?bundle=` names the bundles to load.
- **SINE identity threshold.** `ComPair.sh` needs > 100 identical nt between the two SINEs and with
  the reference (`CP_SINE_NID`). With `CP_SINE_FRAC` (bank mode default 50) the threshold is lowered
  to that percentage of the reference length when that is smaller, e.g. 52 nt for a 105 bp
  consensus or a 105 bp partial span. Consensi of 200 bp or more keep 100 nt.

Toy test (`tests/simulate_bank.py`: consensi of 110, 250, 360 = 110 + 250 and 600 bp; full,
partial and relabelled copies, shared or in one genome; 330 events):

    python3 tests/simulate_bank.py
    bash SINE_orth_loc.bash aaa.bnk bbb.bnk bank.fa && python3 tests/score_bank.py . .

Bank mode: 330/330 events with the right status, no wrong pairs. The same toy run with the 600 bp
consensus alone (`D.fa`, `aaa-D.bed`, `bbb-D.bed`): 70/330.

## Flank length, mapper, one or both flanks (testing)

    FLANK=500 MAPPER=minimap2 MM2_OPTS="-x sr -k 11 -w 3" BOTHFLANKS=1 SINE_orth_loc.bash ...

- `FLANK` (300): flank length; copies closer than this form one locus; passed to the rescue and the
  nested-locus steps.
- `MAPPER` (bwa) / `MM2_OPTS` (`-x sr -k 11 -w 3 -p 0.5 -N 50`): the mapper of the flanks
  (`minimap2 -a --secondary=no`), of the rescue and of the nested-locus test (`minimap2 -a --secondary=yes
  -N 50`, the counterpart of `bwa mem -a`). With minimap2 no bwa index is needed: a `.mmi` index is built
  once per genome in the run directory (`<genome>.mmi`, tagged with `MM2_OPTS` in `<genome>.mmi.opts` and
  rebuilt when the options change; about 9 GB and 1-2 min for a 2 Gbp genome).
- `BOTHFLANKS` (1): maps the right flanks too; a left-flank hit is kept only where the right flank of
  the same copy maps on the same target sequence and strand, -30 to SINElength + 100 bp downstream,
  or where that cannot be tested: the copy's own right flank is shorter than `MINRF` (100 bp; a copy
  at a contig end of its genome, tagged `RF=short` in the hit's MAPQ field) or the hit lies closer
  than SINElength + 100 + `MINRF` bp to the target contig end on the right-flank side (`RF=edge`).
  `BOTHFLANKS=0` maps the left flanks only (the behaviour before 2026-10-30).
- `ADAPT` (`600 1000`) / `ADAPT_MAPQ` (20): adaptive flank extension. A flank left unresolved after
  mapping (no hit, or MAPQ below `ADAPT_MAPQ`; typically a flank inside a repeat) is cut again from its
  genome at each length of `ADAPT` in turn (shorter at a contig end) and mapped again with the same mapper
  and index; a hit with MAPQ >= `ADAPT_MAPQ` (and above the old one) replaces the old hit. The new hit is
  cut back to `FLANK` bp at its far side, keeping the SINE-side end, and tagged `EXT=<length>` in the
  MAPQ field, so `both_flanks` and every later step see `FLANK`-bp hits as before. Applied to left and
  right flanks in both directions before the two-flank check; files `<hits>.ext<L>.{in.bed,bnk,sam,bed}`.
  The extension stops at the nearest other annotated copy of the species (`<sp>-<SINE>.bed`, clustered
  copies included), so an extended flank never contains another SINE copy (flanks left at `FLANK` + 50 bp
  or shorter are not mapped again), and `EXT=<n>` gives the length n actually reached. `ADAPT=""` turns it off.

Toy (`tests/simulate_flanks.py`, 450 events: 266 plain, 63 with a repeat as left flank, 55 segmental
duplicates in bbb, 40 contig ends, 26 right-flank indels; `RESCUE=0 NEST=0`), events with the right
status, `tests/score_flanks.py`. No setting gave a wrong status or a wrong pair; one-flank runs at 300 bp
gave 1-2 paralog pairs (segmental duplicate), both-flank runs none.

| flank divergence | setting | 300, one | 300, both | 500, one | 500, both |
|---|---|---|---|---|---|
| ~6 % (`--div 0.03`) | bwa | 375 | 381 | 389 | 388 |
| | minimap2 `-x sr` | 368 | 375 | 389 | 388 |
| | minimap2 `-x sr -k 11 -w 3` | 375 | 379 | 389 | 388 |
| | minimap2 `-x asm20` | 304 | 226 | 383 | 356 |
| ~16 % (`--div 0.08`) | bwa | 302 | 269 | 335 | 333 |
| | minimap2 `-x sr` | 101 | 22 | 166 | 69 |
| | minimap2 `-x sr -k 11 -w 3` | 336 | 345 | 349 | 361 |
| | minimap2 `-x asm20` | 0 | 0 | 15 | 0 |

Per class: 500 bp flanks recover the right-flank indels (20 -> 26); both flanks recover repeat-flank loci
(bwa 500 bp, 16 %: 14 -> 37) and, before the contig-end rule, lost contig-end loci (26 -> 17); segmental
duplicates stay mostly unresolved without the rescue (~10/55).

With the contig-end rule (2026-10-30; the table above is without it) the both-flank runs keep every
contig-end locus the one-flank runs find and lose nothing else: contig ends 24 -> 28 of 40 at ~6 %,
17 -> 23 (bwa 300), 17 -> 25 (bwa 500), 21 -> 26 (minimap2 -k 11, 300) at ~16 %; all events right
385 / 392 (bwa 300 / 500, ~6 %) and 275 / 341 / 350 (bwa 300 / bwa 500 / minimap2 -k 11 300, ~16 %),
no wrong status or pair. Both flanks is therefore the default.

minimap2 tuning for repeat flanks (2026-10-30, both flanks, 300 bp, `-x sr -k 11 -w 3` plus the option
named; contig ends / all events right): lifting the minimizer frequency cutoff (`-f 0 -U 1000000`) changes
nothing (the toy repeat, ~70 copies, never reached the cutoff of 1000); `-w 1` adds at most one locus;
`-k 13 -w 5` equals bwa at ~6 % (385, repeat flanks 62) but collapses at ~16 % (323, 26);
`-p 0.5 -N 50` is the best: repeat flanks 62 (= bwa) and all 384 (bwa 385) at ~6 %, 35 and 352 at ~16 %
(bwa 26 and 275). Segmental duplicates stay bwa's: 14 vs 7 at ~16 %. No time cost (208-224 s per run for
every variant). So `MAPPER=minimap2 MM2_OPTS="-x sr -k 11 -w 3 -p 0.5 -N 50"` is the minimap2 setting
(the default `MM2_OPTS`); bwa remains the default mapper.

Adaptive flank extension (2026-10-30, `ADAPT="600 1000"`, `ADAPT_MAPQ=20`, both flanks, 300 bp, `RESCUE=1
NEST=1`; events right, ADAPT off -> on; the off runs reproduce the previous runs exactly): at ~6 % bwa 430 ->
432 and minimap2 429 -> 432 (repeat flanks 62 -> 63, segmental duplicates 54 -> 55 / 53 -> 55); at ~16 % bwa
302 -> 416 (repeat flanks 26 -> 59, plain 196 -> 261, segmental duplicates 41 -> 52, contig ends 23 -> 25,
right-flank indels 16 -> 19) and minimap2 397 -> 410 (repeat flanks 35 -> 44, segmental duplicates 52 -> 54).
No wrong status or wrong pair. bwa at ~16 % gains 2 paralog rows (segmental duplicates, from the rescue;
the same two events minimap2 pairs with the duplicate with or without ADAPT, which keeps its 2). Every
converted hit is 300 bp.

### minimap2 settings, from the code

Why these flags, from minimap2 2.31 (`file:line`) and bwa 0.7.19, with the index measured on the bat
genomes (rsi 2.12 Gbp, rfe 2.08 Gbp; 2026-10-28 run):

- `-x sr -k 11 -w 3`: the `sr` preset keeps k = 21, w = 11 unless overridden; a 300 bp flank at 84 %
  identity shares about one 21-mer minimizer with its ortholog (280 x 0.84^21 x 2/12) but about 20
  11-mer minimizers at w = 3, and a chain needs 2 anchors and a score of 25 (`-n 2 -m 25`,
  `options.c:151-168`, `lchain.c:52,67`). That is the collapse of plain `-x sr` at ~16 % in the table
  above (101 / 22) and its absence with `-k 11 -w 3`.
- Repeats in the flanks: `sr` fixes the seed-occurrence cutoff at 1000 (`mid_occ`, `options.c:166`; `-f`
  is ignored). On a 2 Gbp genome at k 11 w 3 the index holds 1.8 M distinct minimizers with 600
  occurrences on average (0.4 % singletons), so a flank minimizer inside a repeat, or any 11-mer more
  common than average, is dropped (`seed.c:109-111`), the chain forms from the rarer half, and a flank
  with no chain is re-seeded at 5000 occurrences (`map.c:293-316`). The `rl:i` SAM tag gives the flank
  bases covered by dropped seeds. This is why the toy's `-f 0 -U 1000000` changed nothing and why the
  real pair maps the same number of flanks as bwa (46,785 vs 46,806 of 50,776 rsi left flanks on rfe).
- `-p 0.5 -N 50`: both act on chain scores before alignment (`map.c:206-213`, `hit.c:265-267`); the
  candidates kept are then aligned and the primary is chosen again by alignment score (`map.c:219-223`,
  `align.c:1118`). More candidates through the DP means a locus whose chain is weaker but whose alignment
  is better can win; `--secondary=no` only stops the printing (`map.c:612-613`). Hence the gain on repeat
  flanks (26 -> 35 at ~16 %) at no time cost.
- Rescue and nested-locus test: every placement of a flank is wanted. bwa `-a` prints the secondaries that
  score >= 0.5 x the primary (`-D 0.5`, `bwamem.c:1050`), flag 256, SEQ `*`, with CIGAR and `AS`;
  minimap2 prints at most `-N` secondaries with a chain score >= `-p` x the primary's (`hit.c:265-267`),
  flag 256, SEQ `*`, with CIGAR and `AS` (`format.c:556,587`), once `--secondary=yes` lifts the preset's
  `MM_F_NO_PRINT_2ND`. The rescue uses the 10 best hits of each flank as anchors and re-aligns both flanks
  itself; the nested-locus test keeps a pair only when it leads by 5 % in `AS`. `AS` is on a different
  scale in the two mappers (`-A 2 -B 8` under `sr`, `-A 1 -B 4` in bwa) and is only ever compared within
  one run, never against a threshold.
- Port test on the flank toy (2026-10-30, both flanks, 300 bp, rescue and nested-locus test on, each mapper
  end to end): events right 430 (bwa) / 429 (minimap2) of 450 at ~6 %, 302 / 397 at ~16 %; repeat flanks
  62 / 62 and 26 / 35; segmental duplicates 54 / 53 and 41 / 52. Both mappers produce the same kind of
  wrong rows there (45 / 45 and 27 / 44 events), all of one type: in a duplication that carries the SINE
  in both copies, the rescue pairs the duplicate's copy with the one site of the other genome too, so one
  site gets two partners. This comes from the rescue, which has no "one partner per site" rule, not from
  the mapper; before the rescue existed these events were simply unresolved (9-14 right). The rescue SAM
  of minimap2 held one secondary hit per primary on the toy (62 / 62), bwa 37 / 38. Run times equal
  (116-147 s).
- Index: a `.mmi` stores k, w and HPC and overrides the command line with a warning (`index.c:626-627`),
  so the index is tagged with `MM2_OPTS`. Building it took 40-85 s and 19-27 GB of memory per bat genome
  (the whole flank mapping of 50 k flanks took 53-110 s with 20 threads, index included); `bwa index` of
  the same genome took 34 min (BWT construction 23 min) for 3.7 GB of index files. Supplementary pieces (`-x sr` splits at indels wider than its
  band of 100) are left out by the pipeline in both mappers (flag 2048), as before.
- Real pair (rsi-rfe, 2.1 Gbp each, bank of 21 consensi, 300 bp, both flanks, 40 threads; copy +-10 kb
  check of `tests/compare_real.py`: rows, of which the partner copy agrees / disagrees / is untested):

  | run | rows | agree | disagree | untested | pipeline time from bare genomes |
  |---|---|---|---|---|---|
  | bwa, `ADAPT=""` | 51,075 | 48,870 | 1,179 | 1,026 | 92 min + 2 x 34 min `bwa index` = ~2 h 40 |
  | minimap2, `ADAPT=""` | 51,037 | 48,806 | 1,201 | 1,030 | 90 min + 2 x 38 s index = ~91 min |
  | minimap2, `ADAPT="600 1000"` | 51,090 | 48,944 | 1,023 | 1,123 | ~110 min |
  | minimap2, `ADAPT` + stop at copies | 51,083 | 48,939 | 1,022 | 1,122 | ~114 min |

  The two mappers alone are equal within 0.1 %; the adaptive extension (about 5,000 of 50,000 flanks per
  direction unresolved at 300 bp; 600 bp resolved ~3,100 of them, 1,000 bp ~500 more, 35-86 unmapped flanks
  per direction newly placed; ~6 % of the final hits carry `EXT=`) replaces poorly supported rows by well
  supported ones (rows only in the adaptive run: 542, 310 agree / 117 disagree; only in the plain run: 489,
  191 / 255). The stop at neighbouring copies truncated 109-136 extended flanks per direction and pass
  (20-24 of them to nothing) and cost 20 rows (17 agreeing) for 4 gained: the price of never mapping a
  flank that contains another copy. The rescue's one-partner rule fired 59 times. Each flank mapping with
  the prebuilt `.mmi` took 5 s; the pipeline time is the alignments, not the mapping.

## Multi-copy loci (multimappers)

When a flank matches several places, loci join into clusters of 3–10 (multi) or more
(poly). The multi stage keeps at most one pair per cluster; poly clusters were never
analysed. Since v2.1, every copy of a multi-copy cluster that is not yet in a final
alignment goes through a **two-flank rescue** (`rescue_multi.py`): left and right flanks
are mapped with all their placements (`bwa mem -a`, or minimap2 with its secondary hits), the best
hits of each flank serve as anchors, and the other
flank is aligned in the window next to each anchor where it must lie (same strand, at
the distance of an empty site or of a SINE). The placement whose two flanks score clearly
best (`--margin`, default 20) is accepted and checked with ComPair.sh like a double
(clusters `M<n>R`); ties are reported as ambiguous instead of guessed. A site has one partner:
placements are taken best first, and a copy whose site, or whose target site, already belongs to
an accepted pair is reported as `shared_site:<pair>` instead of paired (a duplication carrying the
SINE in both copies gives one orthologous pair and one `shared_site` copy, not two pairs to one site).

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
