# Pan-SINEome format, version 0.1

A pan-SINEome is the set of SINE copies of several genomes, organised into
**orthologous groups**: one group per ancestral insertion, listing for every genome
the copy, the empty insertion site, or the lack of data. No genome is a reference.
Every genome is described in its own assembly coordinates, and groups are derived
from pairwise orthology evidence that anyone can add to.

`sine_registry.py build` writes these tables from the per-pair `orth_<a>-<b>.tsv`
files of SINE_orth_loc and, optionally, the BED files of all annotated copies.

## Concepts

- **Copy**: one annotated SINE copy in one genome (from the copy search, e.g. ssearch36).
- **Site**: the position of an insertion in one genome: the 5' junction between the
  left flank and the SINE, in SINE orientation. For a copy on `+` it is the copy's
  start, for `-` its end. A genome without the SINE has the same site as an empty
  junction between the orthologous flanks.
- **Evidence**: one validated pairwise alignment (a row of an `orth` table) linking a
  site in genome A to a site in genome B, with presence or absence of the SINE in
  each and the ComPair.sh metrics.
- **Group**: the connected set of sites linked by evidence, after merging the sites of
  one genome that lie within `--tol` bp (default 60) on the same chromosome and strand.

## Conventions

- Coordinates are 0-based, half-open (BED), as in the pipeline outputs.
- A copy is written `chrom:start-end(strand)`; a site without a copy `chrom:pos(strand)`.
- A copy ID is `<species>:<chrom>:<start>-<end>(<strand>)`.
- Species names are the short names used in the pairwise runs (e.g. `dva`, `pmu`).

## States

| State | Meaning |
| :- | :- |
| P | SINE present at the site |
| A | empty site: orthologous flanks, no SINE |
| U | no data: no validated evidence for this genome |
| X | ambiguous: contradicting calls (`inconsistent:<sp>`), or several sites of the genome joined into the group (`multicopy:<sp>`) |

## Tables

All tables are tab-separated with one header line. `PREFIX` is the build's output prefix.

### `PREFIX.groups.tsv`

| Column | Content |
| :- | :- |
| group | group ID (see below) |
| family | most frequent family of the group's copies (`.` without `--copies`) |
| pattern | states in species order, e.g. `PPAU` |
| n_P, n_A, n_U, n_X | number of genomes in each state |
| pairs | number of evidence rows |
| flags | comma-separated: `multicopy:<sp>`, `inconsistent:<sp>`, `unannotated:<sp>` (P site without an annotated copy), `family_mixed:<f1>/<f2>`; `.` if none |
| one column per species | the copy, or the site (`chrom:pos(strand)`) for A/unannotated P; comma-separated when several; `.` for U |

### `PREFIX.copies.tsv` (with `--copies`)

Every annotated copy, whether or not it is in a group.

| Column | Content |
| :- | :- |
| copy | copy ID |
| species, chrom, start, end, strand | location |
| family | from the BED name column |
| group | group ID, or `.` |
| status | `grouped`; `close_copy` (another copy within 300 bp, excluded by the pipeline); `no_validated_pair` (compared, but no evidence passed QC); `not_compared` (genome not in any orth table) |

### `PREFIX.evidence.tsv`

The `orth` rows the groups were built from, with the group they ended in and the
source table: `group`, `source`, then the orth columns (`alignment`, `cluster`, `status`,
`species1`, `locus1`, `sine1`, `species2`, `locus2`, `sine2` and the ComPair.sh metrics).
Rebuilding from the orth tables reproduces the groups.

### `PREFIX.aliases.tsv`

Group IDs of the previous build that no longer exist as such:
`alias`, `group` (the current ID, `.` if none), `event`:

- `merged`: the alias's group was joined into `group`;
- `split`: part of the alias's group now forms the new `group`;
- `retired`: no current group shares a site with it.

### `PREFIX.matrix.tsv`, `PREFIX.patterns.tsv`, `PREFIX.nex`

Group × species states, the number of groups per pattern, and a NEXUS 0/1 matrix of
the variable groups (at least one P and one A; U and X coded `?`).

## Group IDs

IDs are `<prefix><7 digits>` (default prefix `PSG`). A fresh build numbers groups in
coordinate order. A build with `--previous OLD` keeps IDs:

1. each previous group goes to the new group that shares most of its sites (within `--tol`);
2. a new group that receives several previous IDs keeps the lowest; the others become
   `merged` aliases;
3. a new group that receives none gets the next unused number; previous groups it
   shares sites with are recorded as `split`;
4. aliases of earlier builds are carried over.

Cite group IDs together with the build they come from; resolve old IDs through the
aliases table.

## Adding a genome

1. Annotate its copies (BED6, name = family).
2. Run SINE_orth_loc against the genomes already included (all of them, or at least one
   close relative and one from each other clade).
3. `sine_registry.py build <all orth tables> --copies ... --previous <last build> -o <new build>`.

## Planned for later versions

- One comparison per new genome against a shared index of representative flanks of
  all groups, instead of pairwise runs.
- Verification of empty sites (distance between flanks, no residual SINE sequence).
- Finer `copies` statuses from the ComPair.sh stat files (`shortRF`, `badLF`, ...).
- Family assignment confidence (identity and coverage against each consensus).
- Re-using a previous ID when a group that was merged by mistake is split again.
