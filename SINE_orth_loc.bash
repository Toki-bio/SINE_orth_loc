#!/bin/bash

# Finding orthologius loci with SINE ("SINEname") insertions between two genome assemblies ("sp1" and "sp2")
# Arguments: $1 - genome1, $2-genome2, $3 - SINE consensus. Genomes (preferably) should be named with 3 symbols + .bnk extension, like "mmu.bnk", "hsa.bnk"'
# bed file with SINE coordinates must be named "sp1"-"SINEname".bed "sp2"-"SINEname".bed
# SINE consensus file name must exactly match it's FASTA name

# checking requirements

if [[ $# -eq 0 ]]; then
    echo 'No arguments supplied! $1 - genome1, $2-genome2, $3 - SINE consensus. Genomes (preferably) should be named like "mmu.bnk", "hsa.bnk"'
    echo bed file with SINE coordinates should be named "sp1"-"SINEname".bed "sp2"-"SINEname".bed
    exit 1
fi

THREADS=${THREADS:-$(nproc)} # set THREADS to match the scheduler allocation on clusters
RESCUE=${RESCUE:-1}          # 0: skip the two-flank rescue of unresolved multi-copy clusters
RESCUE_PY="$(dirname "$(readlink -f "$0")")/rescue_multi.py"
NEST=${NEST:-1}              # 0: skip nested/close/satellite analysis (sine_nest.py; needs nhmmer)
NEST_PY="$(dirname "$(readlink -f "$0")")/sine_nest.py"
PYTHON=${PYTHON:-python3}    # Python >= 3.7 (e.g. PYTHON=/usr/local/bin/python3.12 where python3 is older)
FLANK=${FLANK:-300}          # flank length (bp) mapped to the other genome; also the distance below which copies form one locus
MAPPER=${MAPPER:-minimap2}   # mapper for the flanks, the rescue and the nested-locus test: minimap2 (minimap2 -a
                             # $MM2_OPTS; a .mmi index is built once per genome, ~40 s) or bwa (bwa mem, needs a
                             # bwa index: ~34 min per 2 Gbp genome; same results, less memory, README: minimap2 settings)
MM2_OPTS=${MM2_OPTS:--x sr -k 11 -w 3 -p 0.5 -N 50}  # minimap2 options (README: minimap2 settings)
BOTHFLANKS=${BOTHFLANKS:-1}  # 1 (default): map the right flanks too and keep a left-flank hit only where the right
                             # flank of the same copy maps on the same target sequence and strand, downstream of it
                             # (-30 to SINElength + 100 bp away), or where it cannot be tested (contig end, see
                             # both_flanks); 0: left flanks only
ADAPT="${ADAPT-600 1000}"    # adaptive flank extension: flanks left unresolved (no hit, or MAPQ < ADAPT_MAPQ) are cut
                             # again at these lengths (bp, in order) and mapped again; a better hit replaces the old
                             # one, cut back to FLANK bp at its SINE side (adapt_flanks); ADAPT="" turns it off
ADAPT_MAPQ=${ADAPT_MAPQ:-20} # MAPQ below which a flank hit counts as unresolved
# A pair can be run in pieces (e.g. on CI runners or array jobs with a time limit), all in one directory
# or with the directory copied between steps:
#   STAGE=prep                 flank mapping and clustering; writes the batch files (PART*, MULTIPART*)
#   STAGE=align SHARD=k/N      aligns and checks batches k, k+N, k+2N, ... into shard_<k>/ (needs only
#                              the .fai files, the BED files and the consensus; no genome sequence)
#   STAGE=finish               merges shard_*/, then rescue, tables and the nested-locus analysis
# STAGE=all (default) runs everything in one go, as before.
STAGE=${STAGE:-all}
SHARD=${SHARD:-1/1}
case "$STAGE" in all|prep|align|finish) ;; *) echo "STAGE must be all, prep, align or finish"; exit 1;; esac
if ! [[ "$SHARD" =~ ^[0-9]+/[0-9]+$ ]] || (( ${SHARD%/*} < 1 || ${SHARD%/*} > ${SHARD#*/} )); then
    echo "SHARD must be k/N with 1 <= k <= N"; exit 1
fi

case "$MAPPER" in bwa|minimap2) ;; *) echo "MAPPER must be bwa or minimap2"; exit 1;; esac
for tool in mafft esl-alipid seqkit bedtools samtools sam2bed "$MAPPER" ComPair.sh; do
    if ! command -v "$tool" > /dev/null 2>&1; then
    echo "$tool: not found"
    exit 1
    fi
done
if ! { [ -f "$RESCUE_PY" ] && "$PYTHON" -c 'import sys; sys.exit(sys.version_info < (3, 7))' 2> /dev/null; }; then
    echo "$PYTHON is not Python >= 3.7 or $RESCUE_PY is missing: set PYTHON=/path/to/python3.7+"; exit 1
fi
[[ "$FLANK" =~ ^[0-9]+$ ]] && (( FLANK >= 100 )) || { echo "FLANK must be a number >= 100"; exit 1; }
for L in $ADAPT; do [[ "$L" =~ ^[0-9]+$ ]] && (( L > FLANK )) || { echo "ADAPT lengths must be numbers > FLANK"; exit 1; }; done
# the cluster grouping step uses gawk-only features; with mawk it silently finds no loci
if ! awk --version 2>/dev/null | grep -q "GNU Awk"; then echo "awk is not GNU awk (gawk): install gawk"; exit 1; fi

seqkit seq -w 0 $3 > SINE.q # linearising SINE sequence
SINEname="${3%.*}"
# Bank mode: $3 holds several consensi (e.g. a SINEderella bank). Column 4 of both BED files then names the
# consensus of each copy ("subfam", or "subfam:a-b" for a partial copy covering positions a-b), and every
# alignment gets the consensus (or span) of its own copy as reference row, named >$SINEname.
BANK=$([[ $(grep -c '^>' "$3") -gt 1 ]] && echo 1 || echo 0)
if [[ "$BANK" == 1 ]]; then
    SINElength=$(seqkit fx2tab -n -l "$3" | awk -F'\t' '$2>m {m=$2} END {print m}') # longest consensus: window size
    export CP_SINE_FRAC=${CP_SINE_FRAC:-50}
    echo "bank mode: $(grep -c '^>' "$3") consensi, windows for $SINElength bp, CP_SINE_FRAC=$CP_SINE_FRAC"
else
    SINElength=$(seqkit seq -w 0 $3 | awk 'NR==2 {print length}') # calculating SINE length
    # the multi stage finds the consensus row in alignments by this name
    if [[ "$(head -n1 "$3" | sed 's/^>//' | awk '{print $1}')" != "$SINEname" ]]; then
        echo "SINE consensus header must equal the file name without extension (>$SINEname in $3)"; exit 1
    fi
fi

sp1="${1%.*}"
sp2="${2%.*}"

# prepare genomes - renaming to simple 3-letters FASTA names

# indexing genome files
echo "Checking FASTA indexes"

FastaIndexsp1=$(echo $1.fai)
FastaIndexsp2=$(echo $2.fai)

if [ -f "$FastaIndexsp1" ]; then echo $1 index found
    else echo indexing $1
    samtools faidx $1
fi

if [ -f "$FastaIndexsp2" ]; then echo $2 index found
    else echo indexing $2
    samtools faidx $2
fi

# SINE search

SINEs_in_genome1=$(echo $sp1-$SINEname.bed)
SINEs_in_genome2=$(echo $sp2-$SINEname.bed)

if [ -f "$SINEs_in_genome1" ]; then echo Bank of $(wc -l $SINEs_in_genome1) SINEs in $1 found
    else echo no SINE bed file found - please provide $SINEs_in_genome1 file or use sear2kloop script
    exit 1
fi

if [ -f "$SINEs_in_genome2" ]; then echo Bank of $(wc -l $SINEs_in_genome2) SINEs in $2 found
    else echo no SINE bed file found - please provide $SINEs_in_genome2 file or use sear2kloop script
    exit 1
fi

if [[ "$STAGE" == all || "$STAGE" == prep ]]; then
# preparing left flanks

echo clustering coordinates
sort -k1,1 -k2,2n $sp1-$SINEname.bed | bedtools cluster -d $FLANK -i - > $sp1-"$SINEname"_clust.bed
uniq -f6 -u $sp1-"$SINEname"_clust.bed > $sp1-"$SINEname"_uniq.bed # SINE loci at distance more than FLANK bp from each other
uniq -f6 -D $sp1-"$SINEname"_clust.bed > $sp1-"$SINEname"_clusters.bed # SINE loci that are within FLANK bp from each other
sort -k1,1 -k2,2n $sp2-$SINEname.bed | bedtools cluster -d $FLANK -i - > $sp2-"$SINEname"_clust.bed
uniq -f6 -u $sp2-"$SINEname"_clust.bed > $sp2-"$SINEname"_uniq.bed
uniq -f6 -D $sp2-"$SINEname"_clust.bed > $sp2-"$SINEname"_clusters.bed

echo "creating ${FLANK}bp flanks ONLY for SINE copies more than $FLANK bp apart"
for sp in $sp1 $sp2; do
    g=$1; [[ "$sp" == "$sp2" ]] && g=$2
    bedtools flank -s -i $sp-"$SINEname"_uniq.bed -g $g.fai -l $FLANK -r 0 > "$sp"_flanks_$FLANK.bed # left flanks of all SINE loci
    bedtools getfasta -s -bed "$sp"_flanks_$FLANK.bed -fi $g > "$sp"_flanks_$FLANK.bnk
    if [[ "$BOTHFLANKS" == 1 ]]; then
        # right flanks, named by the left flank of the same copy (the name the left-flank hits carry)
        awk -v F=$FLANK 'BEGIN {OFS="\t"} NR==FNR {len[$1]=$2; next}
         {if ($6=="+") {ls=($2>F)?$2-F:0; le=$2; rs=$3; re=($3+F<len[$1])?$3+F:len[$1]}
          else {ls=$3; le=($3+F<len[$1])?$3+F:len[$1]; rs=($2>F)?$2-F:0; re=$2}
          if (le>ls && re>rs) print $1, rs, re, $1":"ls"-"le"("$6")", 0, $6}' $g.fai $sp-"$SINEname"_uniq.bed > "$sp"_rflanks_$FLANK.bed
        bedtools getfasta -s -name -bed "$sp"_rflanks_$FLANK.bed -fi $g | sed '/^>/s/::.*//' > "$sp"_rflanks_$FLANK.bnk
    fi
done

# minimap2 index of a genome, built once per run directory with the k-mer and window settings of MM2_OPTS (a .mmi
# fixes k, w and HPC and overrides the command line; about 9 GB and 1-2 min for a 2 Gbp genome at -k 11 -w 3),
# rebuilt when MM2_OPTS changed
mm2_index() {  # $1: genome; writes $1.mmi and $1.mmi.opts
    if [ -f "$1.mmi" ] && [ "$(cat "$1.mmi.opts" 2> /dev/null)" == "$MM2_OPTS" ]; then echo "$1.mmi found"; else
        echo "indexing $1 (minimap2 $MM2_OPTS)"
        minimap2 $MM2_OPTS -t $THREADS -d "$1.mmi" "$1" && echo "$MM2_OPTS" > "$1.mmi.opts"
    fi
}

# mapping flank sequences, one placement per flank (minimap2, or bwa mem with MAPPER=bwa); hits to BED:
# target, start, end, flank name, MAPQ/LENGTH/AS, strand (SAM records with fewer than 100 bases of SEQ are left out)
map_flanks() {  # $1: genome, $2: flank FASTA, $3: output prefix (writes $3.sam and $3.bed)
    if [ -f "$3.sam" ]; then echo "$3.sam found"; else
        echo "mapping $2 on $1 ($MAPPER)"
        if [[ "$MAPPER" == minimap2 ]]; then minimap2 -a --secondary=no -t $THREADS $MM2_OPTS "$1.mmi" "$2" > "$3.sam"
        else bwa mem -t $THREADS "$1" "$2" > "$3.sam"; fi
    fi
    sam2bed < "$3.sam" | awk 'BEGIN {OFS="\t"} length($12) > 99 {as="AS:i:NA"
     for (i=14; i<=NF; i++) if ($i ~ /^AS:i:/) {as=$i; break}
     sub(/:i:/, "=", as); print $1, $2, $3, $4, "MAPQ="$5",LENGTH="length($12)","as, $6}' | bedtools sort -i > "$3.bed"
}

# all placements of each flank (multi-copy clusters, compound loci): bwa mem -a, or minimap2 with its secondary
# hits printed (the sr preset hides them), up to 50 per primary with a chain score >= 0.5 x the primary's
# (-N 50; -p 0.5 of the sr preset, the counterpart of bwa's -D 0.5); same flags as rescue_multi.py mapper_command()
map_all() {  # $1: genome, $2: FASTA, $3: output SAM
    if [[ "$MAPPER" == minimap2 ]]; then minimap2 -a -N 50 -t $THREADS $MM2_OPTS --secondary=yes "$1.mmi" "$2" > "$3"
    else bwa mem -a -t $THREADS "$1" "$2" > "$3"; fi
}

# BOTHFLANKS=1: keep a left-flank hit only where the right flank of the same copy has a hit on the same target
# sequence and strand, -30 to SINElength + 100 bp downstream of the left-flank hit (in the flank orientation).
# Contig-end rule: a hit whose right flank cannot be tested is kept and tagged in the MAPQ field --
# RF=short: the copy's own right flank is shorter than MINRF bp (copy at a contig end of its genome; flanks
# under 100 bases are not mapped at all); RF=edge: the hit lies closer than SINElength + 100 + MINRF bp to
# the target contig end on the right-flank side, so the right flank has no room to land.
MINRF=${MINRF:-100}
both_flanks() {  # $1: left-flank hits BED (filtered in place), $2: right-flank hits BED, $3: right-flank BED of
                 # the source genome (lengths), $4: target genome .fai
    awk -v M=$((SINElength+100)) -v MINRF=$MINRF 'BEGIN {OFS="\t"}
     FILENAME==ARGV[1] {r[$4]=r[$4] " " $1"|"$2"|"$3"|"$6; next}
     FILENAME==ARGV[2] {rl[$4]=$3-$2; next}
     FILENAME==ARGV[3] {len[$1]=$2; next}
     {if (rl[$4]+0 < MINRF) {$5=$5",RF=short"; short++; print; next}
      room=($6=="+") ? len[$1]-$3 : $2
      if (room < M+MINRF) {$5=$5",RF=edge"; edge++; print; next}
      n=split(r[$4], h, " "); ok=0
      for (i=1; i<=n; i++) {split(h[i], x, "|"); if (x[1]!=$1 || x[4]!=$6) continue
       d=($6=="+") ? x[2]-$3 : $2-x[3]; if (d>=-30 && d<=M) ok=1}
      if (ok) print; else bad++}
     END {print FILENAME ": " bad+0 " left-flank hits without a matching right-flank hit removed; kept untested: " \
          short+0 " RF=short, " edge+0 " RF=edge" > "/dev/stderr"}' "$2" "$3" "$4" "$1" > "$1.tmp" &&
     mv "$1.tmp" "$1"
}

# ADAPT: flanks still unresolved (no hit, or best MAPQ < ADAPT_MAPQ) are cut again at each length L of ADAPT
# from the source genome, under their FLANK-bp name, and mapped again; a hit with MAPQ >= ADAPT_MAPQ (and above
# the old one) replaces the old hit, cut back to FLANK bp at the far side so that its SINE-side end is kept
# (left flank: the hit end on +, the start on -; right flank: the start on +, the end on -) and tagged ",EXT=<n>"
# in the MAPQ field. The hits BED is rewritten in place, so both_flanks and all later steps see FLANK-bp hits.
# The extended flank stops at the nearest other annotated copy of the species ($4-$SINEname.bed, clustered copies
# included) on its far side, so it never contains another SINE copy; a flank left at FLANK + 50 bp or shorter is
# not mapped again; EXT=<n> gives the length n actually reached (from $p.in.bed).
adapt_flanks() {  # $1: hits BED (updated in place), $2: flank FASTA, $3: source genome, $4: source species,
                  # $5: target genome, $6: l (left flanks) or r (right flanks)
    local L p
    for L in $ADAPT; do
        p="${1%.bed}.ext$L"
        awk -v Q=$ADAPT_MAPQ 'FNR==NR {split($5, m, /[=,]/); if (m[2]+0 >= Q) ok[$4]=1; next}
         /^>/ {n=substr($1, 2); if (!(n in ok)) print n}' "$1" "$2" > "$p.names"
        # L-bp flanks of the unresolved copies, named by the FLANK-bp left flank (as bedtools flank + getfasta -s),
        # cut back at the far side to the nearest other copy of $4-$SINEname.bed (copies indexed in 10 kb bins);
        # kept only if longer than FLANK + 50 bp. Counts (truncated, dropped after truncation, dropped short at a
        # contig end) go to $p.cnt.
        awk -v F=$FLANK -v L=$L -v side=$6 -v C="$p.cnt" 'BEGIN {OFS="\t"; W=10000} FILENAME==ARGV[1] {u[$1]=1; next}
         FILENAME==ARGV[2] {len[$1]=$2; next}
         FILENAME==ARGV[3] {for (b=int($2/W); b<=int(($3-1)/W); b++) bin[$1,b]=bin[$1,b] " " $2 " " $3; next}
         {lo=($2>L)?$2-L:0; hi=($3+L<len[$1])?$3+L:len[$1]
          if ($6=="+") {n=$1":"(($2>F)?$2-F:0)"-"$2"(+)"; s=(side=="l")?lo:$3; e=(side=="l")?$2:hi}
          else {n=$1":"$3"-"(($3+F<len[$1])?$3+F:len[$1])"(-)"; s=(side=="l")?$3:lo; e=(side=="l")?hi:$2}
          if (!(n in u) || e<=s) next
          low=(($6=="+") == (side=="l")); s0=s; e0=e  # far side: low coordinates (left +, right -) or high
          for (b=int(s0/W); b<=int((e0-1)/W); b++) {k=split(bin[$1,b], c, " ")
           for (i=1; i<k; i+=2) if (c[i+1]+0 > s0 && c[i]+0 < e0) {if (low) {if (c[i+1]+0 > s) s=c[i+1]+0}
                                                                 else {if (c[i]+0 < e) e=c[i]+0}}}
          t=(s!=s0 || e!=e0); nt+=t
          if (e-s <= F+50) {if (t) nd++; else ns++; next}
          print $1, s, e, n, 0, $6}
         END {print nt+0, nd+0, ns+0 > C}' "$p.names" $3.fai $4-"$SINEname".bed $4-"$SINEname"_uniq.bed > "$p.in.bed"
        local nt nd ns; read -r nt nd ns < "$p.cnt"; rm -f "$p.cnt"
        local nu=$(wc -l < "$p.names") nc=$(cut -f4 "$p.in.bed" | sort -u | wc -l); rm -f "$p.names"
        local tr="$nt truncated at a neighbouring copy ($nd dropped), $ns too short at a contig end"
        (( nu == nc + nd + ns )) || echo "$1: L=$L: $nu unresolved flanks but $((nc + nd + ns)) found in $4-${SINEname}_uniq.bed" >&2
        (( nc > 0 )) || { echo "$1: L=$L: $nu unresolved, $tr, 0 improved (MAPQ >= $ADAPT_MAPQ), 0 newly placed" >&2; break; }
        bedtools getfasta -s -name -bed "$p.in.bed" -fi $3 | sed '/^>/s/::.*//' > "$p.bnk"
        map_flanks $5 "$p.bnk" "$p"
        awk -v Q=$ADAPT_MAPQ -v F=$FLANK -v L=$L -v side=$6 -v nu=$nu -v tr="$tr" -v B="$1" 'BEGIN {OFS="\t"}
         function mq(x,  m) {split(x, m, /[=,]/); return m[2]+0}
         FILENAME==ARGV[1] && FNR==NR {if (!($4 in old) || mq($5) > old[$4]) old[$4]=mq($5); next}
         FILENAME==ARGV[2] {xl[$4]=$3-$2; next}  # length of the extended flank actually cut
         FILENAME==ARGV[3] {q=mq($5); if (q >= Q && (!($4 in old) || q > old[$4]) && (!($4 in nq) || q > nq[$4])) {
           if (($6=="+") == (side=="l")) {if ($2 < $3-F) $2=$3-F} else {if ($3 > $2+F) $3=$2+F}  # SINE side kept
           $5=$5",EXT="xl[$4]; nq[$4]=q; hit[$4]=$0}; next}
         !($4 in hit) {print}
         END {for (n in hit) {print hit[n]; m++; k+=!(n in old)}
              print B ": L=" L ": " nu " unresolved, " tr ", " m+0 " improved (MAPQ >= " Q "), " k+0 " newly placed" > "/dev/stderr"}' \
         "$1" "$p.in.bed" "$p.bed" "$1" | bedtools sort -i > "$1.tmp" && mv "$1.tmp" "$1"
    done
}

for g in $1 $2; do
    if [[ "$MAPPER" == minimap2 ]]; then mm2_index $g
    elif [ -f "$g.sa" ]; then echo $g bwa index found; else echo indexing $g; bwa index $g; fi
done
sfx=""; [[ "$MAPPER" == minimap2 ]] && sfx=".mm2"
map_flanks $1 "$sp2"_flanks_$FLANK.bnk "$sp2"_flanks-$sp1$sfx
map_flanks $2 "$sp1"_flanks_$FLANK.bnk "$sp1"_flanks-$sp2$sfx
adapt_flanks "$sp2"_flanks-$sp1$sfx.bed "$sp2"_flanks_$FLANK.bnk $2 $sp2 $1 l
adapt_flanks "$sp1"_flanks-$sp2$sfx.bed "$sp1"_flanks_$FLANK.bnk $1 $sp1 $2 l
if [[ "$BOTHFLANKS" == 1 ]]; then
    map_flanks $1 "$sp2"_rflanks_$FLANK.bnk "$sp2"_rflanks-$sp1$sfx
    map_flanks $2 "$sp1"_rflanks_$FLANK.bnk "$sp1"_rflanks-$sp2$sfx
    adapt_flanks "$sp2"_rflanks-$sp1$sfx.bed "$sp2"_rflanks_$FLANK.bnk $2 $sp2 $1 r
    adapt_flanks "$sp1"_rflanks-$sp2$sfx.bed "$sp1"_rflanks_$FLANK.bnk $1 $sp1 $2 r
    both_flanks "$sp2"_flanks-$sp1$sfx.bed "$sp2"_rflanks-$sp1$sfx.bed "$sp2"_rflanks_$FLANK.bed $1.fai
    both_flanks "$sp1"_flanks-$sp2$sfx.bed "$sp1"_rflanks-$sp2$sfx.bed "$sp1"_rflanks_$FLANK.bed $2.fai
fi
if [[ -n "$sfx" ]]; then cp "$sp2"_flanks-$sp1$sfx.bed "$sp2"_flanks-$sp1.bed; cp "$sp1"_flanks-$sp2$sfx.bed "$sp1"_flanks-$sp2.bed; fi
echo mapping complete

slop=$((SINElength+FLANK))

echo "preparing coordinates"

bedtools slop -s -i "$sp1"_flanks-$sp2.bed -g $2.fai -l 0 -r $slop > "$sp1"_flanks-"$sp2"_slop.bed
awk '{sub(/-/,"\t",$4); sub(/:/,"\t",$4); {gsub(/\)/,"")}; {gsub(/\(/,"\t0\t0\t")}; {print $4,$5,$6,$1":"$2"-"$3"("$11")",$10,$9}}' OFS="\t" "$sp1"_flanks-"$sp2"_slop.bed |
 bedtools sort | bedtools slop -s -i - -g $1.fai -l 0 -r $slop |
 awk -v sp1=$sp1 -v sp2=$sp2 '{print $0"\t"sp1"_mapped_on_"sp2}' > inv_"$sp1"_flanks-"$sp2"_slop.bed

bedtools slop -s -i "$sp2"_flanks-$sp1.bed -g $1.fai -l 0 -r $slop > "$sp2"_flanks-"$sp1"_slop.bed
awk '{sub(/-/,"\t",$4); sub(/:/,"\t",$4); {gsub(/\)/,"")}; {gsub(/\(/,"\t0\t0\t")}; {print $4,$5,$6,$1":"$2"-"$3"("$11")",$10,$9}}' OFS="\t" "$sp2"_flanks-"$sp1"_slop.bed |
 bedtools sort | bedtools slop -s -i - -g $2.fai -l 0 -r $slop |
 awk '{sub(/-/,"\t",$4); sub(/:/,"\t",$4); {gsub(/\)/,"")}; {gsub(/\(/,"\t0\t0\t")}; {print $4,$5,$6,$1":"$2"-"$3"("$11")",$10,$9}}' OFS="\t" |
 bedtools sort |
 awk -v sp1=$sp1 -v sp2=$sp2 '{print $0"\t"sp2"_mapped_on_"sp1}'> "$sp2"_flanks-"$sp1"_slop_2i.bed

cat "$sp2"_flanks-"$sp1"_slop_2i.bed inv_"$sp1"_flanks-"$sp2"_slop.bed | sort -k1,1 -k2,2n > $sp1-"$sp2"_m.bed

echo "clustering"
bedtools cluster -s -i "$sp1"-"$sp2"_m.bed | awk -v sp1=$sp1 'OFS="\t" {print $1,$2,$3,$4,$5,$6,sp1$8"c",$7}'> clu_$sp1.bed

awk '{sub(/-/,"\t",$4); sub(/:/,"\t",$4); {gsub(/\)/,"")}; {gsub(/\(/,"\t0\t0\t")}; {print $4,$5,$6,$1":"$2"-"$3"("$11")",$10,$9,$12}}' OFS="\t" $sp1-"$sp2"_m.bed |
 bedtools sort | bedtools cluster -s -i - | awk -v sp2=$sp2 'OFS="\t" {print $1,$2,$3,$4,$5,$6,sp2$8"c",$7}'> clu_$sp2.bed

awk '{sub(/-/,"\t",$4); sub(/:/,"\t",$4); {gsub(/\)/,"")}; {gsub(/\(/,"\t0\t0\t")}; {print $4,$5,$6,$1":"$2"-"$3"("$11")",$10,$9,$12,$13}}' OFS="\t" clu_$sp2.bed > "$sp1"_from_clu_$sp2.bed
awk '{sub(/-/,"\t",$4); sub(/:/,"\t",$4); {gsub(/\)/,"")}; {gsub(/\(/,"\t0\t0\t")}; {print $4,$5,$6,$1":"$2"-"$3"("$11")",$10,$9,$12,$13}}' OFS="\t" clu_$sp1.bed > "$sp2"_from_clu_$sp1.bed

cat "$sp1"_from_clu_"$sp2".bed clu_"$sp1".bed | sort -k1,1 -k2,2n > "$sp1"_comb_clu.bed
cat "$sp2"_from_clu_"$sp1".bed clu_"$sp2".bed | sort -k1,1 -k2,2n > "$sp2"_comb_clu.bed
cat "$sp1"_comb_clu.bed "$sp2"_comb_clu.bed > $sp1$sp2

bedtools merge -i "$sp1"_comb_clu.bed -c 4,5,6,7,8 -o collapse,collapse,collapse,collapse,collapse > "$sp1"_comb_clu_mer.bed
bedtools merge -i "$sp2"_comb_clu.bed -c 4,5,6,7,8 -o collapse,collapse,collapse,collapse,collapse > "$sp2"_comb_clu_mer.bed

echo "working with clusters"

awk '{print $7}' "$sp2"_comb_clu_mer.bed "$sp1"_comb_clu_mer.bed > $sp1-"$sp2"_clusters

awk '{gsub(/,/,"\t")}1' $sp1-"$sp2"_clusters |
awk '{delete seen; for(i=1;i<=NF;i++) if (!($i in seen) ) {printf "%s ",$i; seen[$i]=i;}; printf "\n"}' |
awk '!seen[$0]++' |
awk '{
    for ( fldNrA=1; fldNrA<NF; fldNrA++ ) {
        flValA = $fldNrA
        for ( fldNrB=fldNrA+1; fldNrB<=NF; fldNrB++ ) {
            flValB = $fldNrB
            val_pairs[flValA][flValB]
            val_pairs[flValB][flValA]
        }
    }
}

function descend(flValA,       flValB) {
    if ( !seen[flValA]++ ) {
        all_vals[flValA]
        for ( flValB in val_pairs[flValA] ) {
            descend(flValB)
        }
    }
}

END {
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for ( flValA in val_pairs ) {
        delete all_vals
        descend(flValA)
        if ( flValA in all_vals ) {
            sep = ""
            for ( flValB in all_vals ) {
                printf "%s%s", sep, flValB
                sep = OFS
            }
            print ""
        }
    }
}' > $sp1-"$sp2"_clust.uniq.output.clusters

echo "stage1 completed"

echo "preparing coordinates"

awk '{print "C"NR"R",$0}' $sp1-"$sp2"_clust.uniq.output.clusters | awk '{i = 2}{ while ( i <= NF ) { print $1,$i ; i++ } }' > list

awk 'NR==FNR{a[$2]=$1;next} $7 in a {$7=a[$7]}1' OFS="\t" list $sp1$sp2 | awk '{print $1,$2,$3,$7,$5,$6,$4,$8}' OFS="\t" | sort | uniq > $sp1$sp2.bed

bedtools sort -i $sp1$sp2.bed | bedtools merge -s -c 4,5,6 -o distinct,collapse,distinct -i - | awk '{print $1,$2,$3,$4":"$5,"0",$6}' OFS="\t" > $sp1"$sp2"_m.bed

echo "merging banks"

cat $1 $2 > $sp1-$sp2.bnk
samtools faidx $sp1-$sp2.bnk

echo "extracting"

bedtools getfasta -s -name -fi $sp1-"$sp2".bnk -bed $sp1"$sp2"_m.bed |
 awk -F: '{if ($0 ~ ">") {print substr($1,2),">"$4":"$5"::"$2} else print}' |
 awk 'NR%2{printf "%sSeqStart",$0;next;}1' |
 awk '{if ($1 in clusters) clusters[$1] = clusters[$1] "," $2; else clusters[$1] = $2} END {for (cluster in clusters) {print cluster, clusters[cluster]}}' > $sp1-"$sp2"_oneline

echo "separating"

cat $sp1-"$sp2"_oneline | awk -v sp1=$sp1 -v sp2=$sp2 -F">" '{if (NF<=3) print $0 >> sp1"-"sp2"_double"; else if (NF<11) print $0 >> sp1"-"sp2"_multi"; else print $0 >> sp1"-"sp2"_poly"}'

echo "stage2 completed"

# doubles and multies are processed in batches of 100 clusters (the units of STAGE=align shards)
rm -f PART* MULTIPART*
awk 'NR%100==1{x="PART"++i;}{print > x}' $sp1-"$sp2"_double
awk 'NR%100==1{x="MULTIPART"++i;}{print > x}' $sp1-"$sp2"_multi
fi
if [[ "$STAGE" == prep ]]; then
    echo "STAGE=prep done: $(ls PART* 2>/dev/null | wc -l) double and $(ls MULTIPART* 2>/dev/null | wc -l) multi batches"
    exit 0
fi

# Doubles and multies are processed in batches of 100 clusters. Each batch runs
# in its own scratch directory (node-local if SCRATCH_DIR, SLURM_TMPDIR or TMPDIR
# is set), and its alignments are appended to one bundle per category, so the
# shared work directory never holds one file per locus:
#   aln_<sp1>-<sp2>_PM.aln.gz  _MP.aln.gz  _SINE.aln.gz  _rejected.aln.gz
# A bundle is plain text: "##FILE <name>" followed by that alignment (FASTA).
# Use aln_bundle.sh to list/extract, or open them in sine_loci_browser.html.

WD=$(pwd)
SINEfa=$(readlink -f "$3")
REFDIR="$WD/refs" COPIES="$WD/bank_copies.bed"
if [[ "$BANK" == 1 ]]; then
    # one reference file per bank consensus and per BED label (whole consensus, or the span a-b of it), header
    # >$SINEname; containers.tsv: consensus X, longer consensus Y holding >= 80 % of the 16-mers of X (a unit in
    # a composite), used when a copy is longer than the consensus it is labelled with (see pick_refs)
    mkdir -p "$REFDIR"
    cut -f1-6 $sp1-$SINEname.bed $sp2-$SINEname.bed > "$COPIES"
    seqkit seq -w 0 "$3" | awk -v d="$REFDIR" -v n="$SINEname" -F'\t' '
     FNR==NR {if (/^>/) {h=substr($1,2); sub(/[ \t].*/,"",h)} else {seq[h]=toupper($0); f=h; gsub(/[:\/]/,"_",f)
              print ">" n " " h "\n" seq[h] > (d "/" f ".fa"); print h "\t" length(seq[h]) > (d "/lengths.tsv")}; next}
     FNR==1 {for (x in seq) for (y in seq) if (x != y && length(seq[y]) > length(seq[x])) {
               t=0; k=0; for (i=1; i+15<=length(seq[x]); i++) {t++; if (index(seq[y], substr(seq[x], i, 16))) k++}
               if (t && k >= 0.8*t) print x "\t" y > (d "/containers.tsv")}}
     !($4 in done) && !($4 in seq) {done[$4]; l=$4; c=l; a=1; b=0
      if (match(l, /^(.*):([0-9]+)-([0-9]+)$/, m)) {c=m[1]; a=m[2]; b=m[3]}
      if (!(c in seq)) {print "bank mode: consensus " c " (BED label " l ") is not in the bank" > "/dev/stderr"; bad=1; next}
      s=(b ? substr(seq[c], a, b-a+1) : seq[c]); f=l; gsub(/[:\/]/,"_",f)
      print ">" n " " l "\n" s > (d "/" f ".fa"); print l "\t" length(s) > (d "/lengths.tsv")}
     END {exit bad}' - "$COPIES" || exit 1
    sort -u -o "$REFDIR/lengths.tsv" "$REFDIR/lengths.tsv"; touch "$REFDIR/containers.tsv"
    echo "bank mode: $(wc -l < "$REFDIR/lengths.tsv") reference rows in $REFDIR"
fi
cut -f1 "$1.fai" > "$WD/$sp1.chroms"
# sequence sizes of both genomes: ComPair.sh reports loci whose window is clipped by a contig end as contig_end
cut -f1,2 "$1.fai" "$2.fai" > "$WD/$sp1-$sp2.sizes"
export CP_SIZES="$WD/$sp1-$sp2.sizes" PYTHON
if [[ "$STAGE" == align ]]; then
    # this shard's batches go to shard_<k>/, which then serves as the work directory
    k=${SHARD%/*}; N=${SHARD#*/}
    rm -rf "shard_$k"; mkdir -p "shard_$k"
    for i in PART* MULTIPART*; do [ -f "$i" ] || continue
        n=${i##*PART}; (( (n - 1) % N == k - 1 )) && cp "$i" "shard_$k/"
    done
    cp "$sp1.chroms" "shard_$k/"
    cd "shard_$k" && WD=$(pwd)
    echo "STAGE=align shard $SHARD: $(ls PART* 2>/dev/null | wc -l) double and $(ls MULTIPART* 2>/dev/null | wc -l) multi batches"
fi
SCRATCH=$(mktemp -d "${SCRATCH_DIR:-${SLURM_TMPDIR:-${TMPDIR:-$WD}}}/SINE_orth_loc_${sp1}-${sp2}.XXXXXX") || exit 1
trap 'rm -rf "$SCRATCH"' EXIT
BUNDLE="$WD/aln_${sp1}-${sp2}"
for t in PM MP SINE rejected; do rm -f "${BUNDLE}_$t.aln.gz"; echo "##SINE_orth_loc bundle v1" > "${BUNDLE}_$t.aln"; done
rm -f statcoords statpairs refs_chosen

# each parallel ComPair.sh job writes its own *.stat file; merge them into "stat"
# serially so concurrent appends cannot interleave (unsafe on NFS/Lustre)
collect_stats() {
    find . -maxdepth 1 -type f -name "*.stat" -exec cat {} + >> "$WD/stat"
    find . -maxdepth 1 -type f -name "*.stat" -delete
}

# coordinates of both loci in every validated alignment of the current batch
# (statcoords: "cluster locus" per sequence; statpairs: alignment, seq1 locus, seq2 locus)
collect_coords() {
    for i in *.cl.*.PM *.cl.*.MP *.cl.*.SINE *.cl.*.contig_end; do [ -f "$i" ] || continue
     [[ "$i" == *.contig_end ]] || awk -F:: 'NR==1 || NR==3 {c=FILENAME; sub(/\..*/,"",c); h=$1; sub(/^>/,"",h); print c, h}' "$i" >> "$WD/statcoords"
     awk -F:: 'NR==1 {h1=$1} NR==3 {h2=$1} END {sub(/^>/,"",h1); sub(/^>/,"",h2); print FILENAME"\t"h1"\t"h2}' "$i" >> "$WD/statpairs"
    done
}

# append every alignment left in the current batch directory to the bundles
bundle_batch() {
    find . -maxdepth 1 -type f ! -name "*.fai" ! -name "*.stats" ! -name "*.off" -print0 |
     xargs -0 -r awk -v b="$BUNDLE" 'FNR==1 {n=FILENAME; sub(/^\.\//,"",n); t=n; sub(/.*\./,"",t)
      if (t!="PM" && t!="MP" && t!="SINE") t="rejected"; out=b"_"t".aln"; print "##FILE " n >> out}
      {print >> out}'
}

# write one <cluster>.cl FASTA (cluster sequences + SINE consensus) per line of a PART file
make_clusters() {
    clstrs=$1
    Lines=$(awk 'END {print NR}' $clstrs)
    while [ $Lines -ge 1 ]
     do
      echo "Processing $Lines line"
      CurrCluster=$(awk -v line=$Lines 'NR == line' $clstrs)
      cName=$(echo $CurrCluster | awk '{print $1}')
      echo $CurrCluster |
       awk '{sub(/ /,"\n")}{gsub(/,>/,"\n>")}{gsub(/SeqStart/,"\n")} {print}' |
       awk 'NR==1 {a=$0} NR!=1 {if ($0 ~ ">") {print $0"::"a} else print}' > "$cName".cl
      [[ "$BANK" == 1 ]] || cat "$SINEfa" >> "$cName".cl
      Lines=$(( $Lines - 1 ))
    done
    [[ "$BANK" == 1 ]] && pick_refs
}

# bank mode: the reference row of a cluster is the consensus of the annotated copy at its windows' insertion
# junction (window start + FLANK on +, end - FLANK on -; nearest copy within FLANK bp); when the loci carry different
# labels, the longest reference is used. A copy more than 10 % longer than its labelled consensus (e.g. a
# composite copy whose votes went to a unit it contains) gets the container consensus closest to its length.
# Without a copy at any junction: the longest consensus of the bank.
# "$WD/refs_chosen": cluster, reference label, labels found (locus=label, or label>container when widened).
# <cluster>.off: locus, start of its annotated copy within the window (ComPair.sh takes the insertion junction
# from it; see compair_flank.py junction).
pick_refs() {
    grep -H '^>' *.cl | awk '{f=$0; sub(/:>.*/,"",f); h=$0; sub(/^[^>]*>/,"",h); sub(/::.*/,"",h)
      if (match(h, /^(.*):([0-9]+)-([0-9]+)\(([+-])\)$/, m)) print m[1]"\t"m[2]"\t"m[3]"\t"f"|"h"\t0\t"m[4]}' > win.bed
    bedtools intersect -a win.bed -b "$COPIES" -wa -wb 2> /dev/null |
     awk -F'\t' -v F=$FLANK -v L="$REFDIR/lengths.tsv" -v CT="$REFDIR/containers.tsv" \
      -v fb="$(sort -k2,2nr "$REFDIR/lengths.tsv" | head -n1 | cut -f1)" '
      BEGIN {while ((getline l < L) > 0) {split(l, x, "\t"); len[x[1]]=x[2]}
             while ((getline l < CT) > 0) {split(l, x, "\t"); cont[x[1]]=cont[x[1]] " " x[2]}}
      {d=($6=="+") ? $8-($2+F) : ($3-F)-$9; if (d<0) d=-d; if ($6!=$12) d+=1000
       if (d<=F && (!($4 in best) || d<best[$4])) {best[$4]=d; lab[$4]=$10; cl[$4]=$9-$8; off[$4]=($6=="+") ? $8-$2 : $3-$9}}
      END {for (w in off) {split(w, p, "|"); c=p[1]; sub(/\.cl$/,"",c); print p[2] "\t" off[w] > (c ".off")}
           for (w in lab) {r=lab[w]; shown=r
            if (cl[w] > 1.1*len[r] && (r in cont)) {n=split(cont[r], y, " "); b=r
             for (i=1; i<=n; i++) {e1=len[y[i]]-cl[w]; e0=len[b]-cl[w]; if ((e1<0?-e1:e1) < (e0<0?-e0:e0)) b=y[i]}
             if (b != r) {lab[w]=b; shown=r ">" b}}
            split(w, p, "|"); c=p[1]; sub(/\.cl$/,"",c); seen[c]=seen[c] (seen[c]=="" ? "" : ",") p[2]"="shown
            if (!(c in ref) || len[lab[w]]>len[ref[c]]) ref[c]=lab[w]}
           for (c in seen) print c"\t"ref[c]"\t"seen[c]
           print "#fallback\t" fb}' > win.refs
    fb=$(awk -F'\t' '$1=="#fallback" {print $2}' win.refs)
    for f in *.cl; do
        c=${f%.cl}; r=$(awk -F'\t' -v c="$c" '$1==c {print $2}' win.refs)
        if [[ -z "$r" ]]; then r=$fb; printf '%s\t%s\t-\n' "$c" "$r" >> "$WD/refs_chosen"
        else awk -F'\t' -v c="$c" '$1==c' win.refs >> "$WD/refs_chosen"; fi
        x=$r; x=${x//:/_}; x=${x//\//_}; cat "$REFDIR/$x.fa" >> "$f"
    done
    rm -f win.bed win.refs
}

# run a per-batch function in a fresh scratch directory, then bundle and clean up
run_batch() {
    local part=$1 func=$2 B="$SCRATCH/$1"
    mkdir -p "$B" && cd "$B" || exit 1
    make_clusters "$WD/$part"
    $func
    collect_stats
    collect_coords
    bundle_batch
    cd "$WD" && rm -rf "$B" && rm "$part"
}

align_doubles() {
    find . -type f -name "*.cl" -print0 |
     xargs -0 -t -I % -P $THREADS sh -c "mafft --op 5 --quiet '%' |
     seqkit seq -w 0 > '%.dbl'; COMPAIR_STAT='%.dbl.stat' ComPair.sh '%.dbl'; rm '%' '%.dbl'"
}

align_multies() {
    find . -type f -name "*.cl" -print0 |
     xargs -0 -t -I % -P $THREADS sh -c "mafft --op 5 --quiet '%' |
     seqkit seq -w 0 > '%.mul'; rm '%'"

    echo finding best pair of sequences in each cluster
    # species of a sequence = genome its chromosome belongs to (chromosome names of the two
    # genomes need not differ in their first characters, e.g. GenBank accessions)

    find . -type f -name "*.mul" -exec sh -c "esl-alipid {} |
     awk -v g1=\"$WD/$sp1.chroms\" -v SINEname=$SINEname 'BEGIN {while ((getline c < g1) > 0) in1[c]}
      function sp(n) {sub(/::.*/, \"\", n); sub(/:[^:]*\$/, \"\", n); return (n in in1) ? 1 : 2}
      NR>1 && \$0 !~ SINEname && sp(\$1) != sp(\$2)' |
     LC_ALL=C sort -r -k 3 -g |
     awk  -v SINEname=$SINEname 'NR==1 {  print \$1,\$2,SINEname,"\n" }' | sed 's/\s\+/\n/g' >  '{}.list'
     seqkit grep -f '{}.list' '{}' > '{}.doub'; rm '{}.list'" \;

    find . -type f -name "*.doub" -print0 | xargs -0 -t -I % -P $THREADS sh -c "mafft --op 5 --quiet '%' | seqkit seq -w 0 > '%.doubles'; COMPAIR_STAT='%.doubles.stat' ComPair.sh '%.doubles'; rm '%' '%.doubles.stat'"

    for i in *.mul.*bad*F *.mul.*shortRF *.mul.*.lfSINE *.mul.*.rfSINE; do [ -f "$i" ] || continue; rm $i "${i//.doub.doubles*/}" ; done

    echo "finding all good sequences in each cluster by similarity in right flanks"

    for i in *.mul.*.MP *.mul.*.PM *.mul.*.SINE; do [ -f "$i" ] || continue;
     name1=$(awk 'NR==1{print substr($1,2)}' "$i");
     name2=$(awk 'NR==3{print substr($1,2)}' "$i");
    echo name1=$name1 name2=$name2;
     SINEend=$(seqkit range -r -1:-1 -w 0 "${i//.doub.doubles*/}" | awk -F '[^-]+' 'END {print length($NF)}');
     ALIGNMENTend=$(seqkit range -r -1:-1 -w 0 "${i//.doub.doubles*/}" | awk 'END {print length}');
     CUT=$((ALIGNMENTend-SINEend+1));
    awk -v CUT=$CUT -v ALIGNMENTend=$ALIGNMENTend 'BEGIN{RS=">";FS="\n"}NR>1{seq="";
    for (i=2;i<=NF;i++) seq=seq""$i; print ">"$1"\n"substr(seq,CUT,ALIGNMENTend-CUT)}' "${i//.doub.doubles*/}" |
     esl-alipid - | awk -v name1="$name1" -v name2="$name2" -v SINEname=$SINEname '$0!~SINEname {if ($1==name1 || $1==name2) print}' |
     awk '{if ($3>40.00) print $1"\n"$2}' |
     awk -v name1="$name1" -v name2="$name2" '{print}END {print name1"\n"name2}' | awk '!seen[$0]++' > $i.list; ListSize=$(awk 'END{print NR}' $i.list)
    if [[ $ListSize -eq 0 ]]; then mv $i $i."sel.aligned"; rm $i.list
     else
     seqkit grep -f $i.list -w 0 "${i//.doub.doubles*/}" > $i.SelSeq; rm $i.list
     seqkit range -r -1:-1 -w 0 $i > $i.SelSINE; cat $i.SelSeq $i.SelSINE > $i.selected;
     rm $i "${i//.doub.doubles*/}" $i.SelSeq $i.SelSINE;
    fi
    done

    for k in *.selected; do [ -f "$k" ] || continue; Count=$(grep -c ">" $k); if [[ $Count -ne "3" ]]; then rm $k; fi; done

    find . -type f -name "*.selected" -print0 |
     xargs -0 -t -I % -P $THREADS sh -c "mafft --op 5 --quiet '%' |
     seqkit seq -w 0 > '%.sel.aligned'; COMPAIR_STAT='%.sel.aligned.stat' ComPair.sh '%.sel.aligned'; rm '%'"

    for i in *.mul.*.MP *.mul.*.PM *.mul.*.SINE; do [ -f "$i" ] || continue; Len=$(awk '{if (NR==2) print length}' $i);
     if [[ $Len -ge 1600 ]]
         then rm $i
     fi
    done
}

if [[ "$STAGE" == all || "$STAGE" == align ]]; then
echo working with doubles

: > stat

for i in PART*; do [ -f "$i" ] || continue
    run_batch "$i" align_doubles
done

mv stat stat_doubles_"$sp1"-"$sp2"

echo "stage3 - aligning doubles - completed"

echo "working with multies"

: > stat
echo aligning clusters

for i in MULTIPART*; do [ -f "$i" ] || continue
    run_batch "$i" align_multies
done

mv stat stat_multi_"$sp1"-"$sp2"
fi
if [[ "$STAGE" == align ]]; then
    for t in PM MP SINE rejected; do gzip -f "${BUNDLE}_$t.aln"; done
    touch statcoords statpairs refs_chosen
    echo "STAGE=align shard $SHARD done"
    exit 0
fi
if [[ "$STAGE" == finish ]]; then
    # merge the shards: stat files, coordinates and alignment bundles (one header per bundle)
    ls -d shard_*/ > /dev/null 2>&1 || { echo "STAGE=finish: no shard_*/ directories"; exit 1; }
    for f in stat_doubles_"$sp1"-"$sp2" stat_multi_"$sp1"-"$sp2" statcoords statpairs refs_chosen; do
        cat shard_*/"$f" > "$f" 2> /dev/null
    done
    for t in PM MP SINE rejected; do
        for z in shard_*/aln_"$sp1"-"$sp2"_$t.aln.gz; do zcat "$z" | tail -n +2; done >> "${BUNDLE}_$t.aln"
    done
    echo "merged $(ls -d shard_*/ | wc -l) shards"
fi

# Multi-copy clusters (3-10 loci, of which the multi stage keeps at most one pair, and
# clusters of more than 10 loci) are usually due to a repetitive left flank. For every SINE
# copy in them that is not yet in a final alignment, map left and right flanks with all their
# placements (map_all) and accept a target only where both flanks agree on one placement
# (rescue_multi.py); accepted pairs are aligned and checked like doubles (clusters M<n>R).
: > stat
if [[ "$RESCUE" == 1 ]]; then
    echo "rescuing unresolved multi-copy clusters with both flanks"
    touch "$sp1-$sp2"_multi "$sp1-$sp2"_poly
    "$PYTHON" "$RESCUE_PY" prepare --clusters "$sp1-$sp2"_multi "$sp1-$sp2"_poly --resolved statpairs \
     --copies "$sp1=$sp1-${SINEname}_uniq.bed" "$sp2=$sp2-${SINEname}_uniq.bed" --genome "$sp1=$1" "$sp2=$2" --sine-length $SINElength --flank $FLANK --out rescue
    [ -s rescue_"$sp1".fa ] && map_all $2 rescue_"$sp1".fa rescue_"$sp1".sam
    [ -s rescue_"$sp2".fa ] && map_all $1 rescue_"$sp2".fa rescue_"$sp2".sam
    "$PYTHON" "$RESCUE_PY" resolve --prep rescue --sam "$sp1=rescue_$sp1.sam" "$sp2=rescue_$sp2.sam" \
     --genome "$sp1=$1" "$sp2=$2" --sine-length $SINElength --flank $FLANK --out "$sp1-$sp2"_rescue --report rescue_"$sp1-$sp2".tsv
    awk 'NR%100==1{x="RPART"++i;}{print > x}' "$sp1-$sp2"_rescue
    for i in RPART*; do [ -f "$i" ] || continue
        run_batch "$i" align_doubles
    done
fi
mv stat stat_rescue_"$sp1"-"$sp2"

echo compressing alignment bundles
for t in PM MP SINE rejected; do gzip -f "${BUNDLE}_$t.aln"; done

echo "gathering statistics"

awk -F"./" '{print $2}' stat_multi_"$sp1"-"$sp2" stat_doubles_"$sp1"-"$sp2" stat_rescue_"$sp1"-"$sp2" | awk '{if ($2=="PM" || $2=="MP" || $2=="SINE") print}' | awk -F. '{print $1,$0}' > statdata

touch statcoords statpairs

cat statcoords statdata | sort | awk '$1 in a {print a[$1]"\n"$0; next} {a[$1]=$0}' | sort | uniq > statcomb

awk '$1 in a {print a[$1]"\n"$2; next} {a[$1]=$0}' statcomb | awk '!seen[$0]++' | tac | awk '{if (NF==1) printf $0" "; else print }' > statbed_"$sp1"-"$sp2"

while read col1 col2 rest; do echo $col2 $col1 $rest; done < statbed_"$sp1"-"$sp2" > statbed_"$sp2"-"$sp1"_reversed

awk '{sub(/\(/,"\t0\t0\t")}{sub(/:/,"\t")}{sub(/)/,"\t")}{sub(/-/,"\t")}'1 statbed_"$sp1"-"$sp2" > statbed_"$sp1"-"$sp2".bed
awk '{sub(/\(/,"\t0\t0\t")}{sub(/:/,"\t")}{sub(/)/,"\t")}{sub(/-/,"\t")}'1 statbed_"$sp2"-"$sp1"_reversed > statbed_"$sp2"-"$sp1".bed

awk '{if (NF>3) print $3}' statcomb | sort | uniq -c > MP_PM_SINE_"$sp1"-"$sp2".txt

# one row per validated alignment, loci in alignment order (seq1 = the sequence SO refers to);
# sine1/sine2: 1 if that locus carries the SINE. Input for sine_registry.py build.
# reference: the consensus (bank mode: subfamily, or span a-b of it) that is the alignment's reference row
echo "writing orth_$sp1-$sp2.tsv"
touch refs_chosen
awk -v sp1=$sp1 -v sp2=$sp2 -v SINEname=$SINEname -v RC=refs_chosen -F'\t' '
 BEGIN {OFS="\t"; n=split("LF FL LFcp RFlength FR RFcp OneTwo OneSINE TwoSINE SL SO ST FI",K," ")
     while ((getline l < RC) > 0) {split(l, x, "\t"); ref[x[1]]=x[2]}
     printf "alignment\tcluster\tstatus\tspecies1\tlocus1\tsine1\tspecies2\tlocus2\tsine2"; for (i=1;i<=n;i++) printf "\t%s", K[i]; print "\tanchor1\tanchor2\treference"}
 FILENAME==ARGV[1] {g1[$1]; next}
 FILENAME ~ /(^|\/)stat_/ {split($0,f," "); b=f[1]; sub(/^\.\//,"",b); delete kv
     for (i=3;i<=length(f);i++) {p=index(f[i],"="); kv[substr(f[i],1,p-1)]=substr(f[i],p+1)}
     m=""; for (i=1;i<=n;i++) m=m"\t"((K[i] in kv) ? kv[K[i]] : "")
     met[b"."f[2]]=m; edge[b"."f[2]]=("EDGE" in kv) ? kv["EDGE"] : -1
     ls1[b"."f[2]]=kv["LS1"]; ls2[b"."f[2]]=kv["LS2"]; next}
 ($1 in met) {st=$1; sub(/.*\./,"",st); cl=$1; sub(/\..*/,"",cl)
     c1=$2; sub(/:[0-9]+-[0-9]+\([+-]\)$/,"",c1); c2=$3; sub(/:[0-9]+-[0-9]+\([+-]\)$/,"",c2)
     if (st=="contig_end") {s1=(edge[$1]==0)?"M":"?"; s2=(edge[$1]==1)?"M":"?"; st="MISSING"}
     else {s1=(st=="PM"||st=="SINE")?1:0; s2=(st=="MP"||st=="SINE")?1:0}
     print $1,cl,st,((c1 in g1)?sp1:sp2),$2,s1,((c2 in g1)?sp1:sp2),$3,s2 met[$1], anchor($2,ls1[$1]), anchor($3,ls2[$1]), ((cl in ref)?ref[cl]:SINEname)}
 # insertion junction in genome coordinates: window start + residues before the SINE (+), window end - them (-)
 function anchor(locus, ls,   m) {if (ls=="") return ""; match(locus, /:([0-9]+)-([0-9]+)\(([+-])\)$/, m)
     return (m[3]=="+") ? m[1]+ls : m[2]-ls}' \
 "$1.fai" stat_doubles_"$sp1"-"$sp2" stat_multi_"$sp1"-"$sp2" stat_rescue_"$sp1"-"$sp2" statpairs > orth_"$sp1"-"$sp2".tsv

# bank mode: the reference row used for every cluster (labels: locus=BED label of the copy at each junction)
if [[ "$BANK" == 1 ]]; then
    { printf 'cluster	reference	labels
'; sort -u refs_chosen; } > refs_"$sp1"-"$sp2".tsv
fi

# what became of every cluster: resolved (and how), rejected by ComPair.sh, or left unresolved
echo "writing clusters_$sp1-$sp2.tsv"
"$PYTHON" "$RESCUE_PY" account --double "$sp1-$sp2"_double --multi "$sp1-$sp2"_multi --poly "$sp1-$sp2"_poly \
 --stat stat_doubles_"$sp1"-"$sp2" stat_multi_"$sp1"-"$sp2" stat_rescue_"$sp1"-"$sp2" --statpairs statpairs \
 --rescue-report rescue_"$sp1-$sp2".tsv --out clusters_"$sp1"-"$sp2".tsv

# nested (SINE-in-SINE), split, dimeric, close and satellite copies: sine_nest.py scan finds them with
# nhmmer (incl. fragments below the annotation thresholds); compound loci are then compared through the
# flanks outside the whole locus, and each element is tested by its junctions (orth rows N<k>R).
# nest_<sp>.caution.bed lists copies whose flank-based calls need caution (sine_registry.py --caution).
if [[ "$NEST" == 1 ]]; then
    if command -v nhmmer > /dev/null; then
        for sp in "$sp1" "$sp2"; do
            g=$1; [[ "$sp" == "$sp2" ]] && g=$2
            [[ -s nest_"$sp".loci.tsv ]] || "$PYTHON" "$NEST_PY" scan --genome "$g" --copies "$sp-$SINEname".bed \
                --library "$SINEfa" -o nest_"$sp" --close $FLANK --cpu "$THREADS"
        done
        idx=""; [[ "$MAPPER" == minimap2 ]] && idx=.mmi
        "$PYTHON" "$NEST_PY" orth --a "$sp1=$1" --b "$sp2=$2" --scan nest_"$sp1" -o nested_"$sp1"-in-"$sp2".tsv \
            --sine-length "$SINElength" --flank $FLANK --threads "$THREADS" \
            --mapper "$MAPPER" --mapper-opts "$MM2_OPTS" --index "$2$idx"
        "$PYTHON" "$NEST_PY" orth --a "$sp2=$2" --b "$sp1=$1" --scan nest_"$sp2" -o nested_"$sp2"-in-"$sp1".tsv \
            --sine-length "$SINElength" --flank $FLANK --threads "$THREADS" \
            --mapper "$MAPPER" --mapper-opts "$MM2_OPTS" --index "$1$idx"
        "$PYTHON" "$NEST_PY" supersede orth_"$sp1"-"$sp2".tsv --scan "$sp1=nest_$sp1" "$sp2=nest_$sp2" \
            --nested nested_"$sp1"-in-"$sp2".tsv nested_"$sp2"-in-"$sp1".tsv --flank $FLANK -o orth.tmp \
            && mv orth.tmp orth_"$sp1"-"$sp2".tsv && mv orth.tmp.superseded orth_"$sp1"-"$sp2".superseded.tsv
        tail -n +2 -q nested_"$sp1"-in-"$sp2".tsv nested_"$sp2"-in-"$sp1".tsv >> orth_"$sp1"-"$sp2".tsv
    else
        echo "nhmmer (HMMER) not found: skipping the nested/close/satellite analysis (NEST=0 to silence)"
    fi
fi

# per-reference summary of the pair (loci, statuses, identities): summary_<sp1>-<sp2>.tsv
"$PYTHON" "$(dirname "$(readlink -f "$0")")/orth_summary.py" orth_"$sp1"-"$sp2".tsv \
    --copies "$sp1=$sp1-$SINEname.bed" "$sp2=$sp2-$SINEname.bed"
# the SINE/PM/MP bundles split per reference, for viewing one subfamily's alignments: by_reference/<reference>/
"$PYTHON" "$(dirname "$(readlink -f "$0")")/orth_split.py" orth_"$sp1"-"$sp2".tsv
