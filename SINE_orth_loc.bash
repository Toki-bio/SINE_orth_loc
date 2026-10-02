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
PYTHON=${PYTHON:-python3}    # Python >= 3.7 (e.g. PYTHON=/usr/local/bin/python3.12 where python3 is older)

for tool in mafft esl-alipid seqkit bedtools samtools sam2bed bwa ComPair.sh; do
    if ! command -v "$tool" > /dev/null 2>&1; then
    echo "$tool: not found"
    exit 1
    fi
done
if ! { [ -f "$RESCUE_PY" ] && "$PYTHON" -c 'import sys; sys.exit(sys.version_info < (3, 7))' 2> /dev/null; }; then
    echo "$PYTHON is not Python >= 3.7 or $RESCUE_PY is missing: set PYTHON=/path/to/python3.7+"; exit 1
fi
# the cluster grouping step uses gawk-only features; with mawk it silently finds no loci
if ! awk --version 2>/dev/null | grep -q "GNU Awk"; then echo "awk is not GNU awk (gawk): install gawk"; exit 1; fi

seqkit seq -w 0 $3 > SINE.q # linearising SINE sequence
SINElength=$(seqkit seq -w 0 $3 | awk 'NR==2 {print length}') # calculating SINE length
SINEname="${3%.*}"
# the multi stage finds the consensus row in alignments by this name
if [[ "$(head -n1 "$3" | sed 's/^>//' | awk '{print $1}')" != "$SINEname" ]]; then
    echo "SINE consensus header must equal the file name without extension (>$SINEname in $3)"; exit 1
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

# preparing left flanks

echo clustering coordinates
sort -k1,1 -k2,2n $sp1-$SINEname.bed | bedtools cluster -d 300 -i - > $sp1-"$SINEname"_clust.bed
uniq -f6 -u $sp1-"$SINEname"_clust.bed > $sp1-"$SINEname"_uniq.bed # SINE loci at distance more than 300bp from each other 
uniq -f6 -D $sp1-"$SINEname"_clust.bed > $sp1-"$SINEname"_clusters.bed # SINE loci that are within 300 bp from each other
sort -k1,1 -k2,2n $sp2-$SINEname.bed | bedtools cluster -d 300 -i - > $sp2-"$SINEname"_clust.bed
uniq -f6 -u $sp2-"$SINEname"_clust.bed > $sp2-"$SINEname"_uniq.bed
uniq -f6 -D $sp2-"$SINEname"_clust.bed > $sp2-"$SINEname"_clusters.bed

echo "creating 300bp flanks ONLY for SINE copies more than 300 bp apart" 
bedtools flank -s -i $sp1-"$SINEname"_uniq.bed -g $1.fai -l 300 -r 0 > "$sp1"_flanks_300.bed # creating 300bp left flanks for all SINE loci
bedtools getfasta -s -bed "$sp1"_flanks_300.bed -fi $1 > "$sp1"_flanks_300.bnk
bedtools flank -s -i $sp2-"$SINEname"_uniq.bed -g $2.fai -l 300 -r 0 > "$sp2"_flanks_300.bed
bedtools getfasta -s -bed "$sp2"_flanks_300.bed -fi $2 > "$sp2"_flanks_300.bnk

# mapping left flanks sequences with bwa mem 

echo "checking if bwa index exists"
bwaIndexsp1=$(echo $1.sa)
bwaIndexsp2=$(echo $2.sa)

if [ -f "$bwaIndexsp1" ]; then echo $1 bwa index found
    else echo indexing $1
    bwa index $1
fi

if [ -f "$bwaIndexsp2" ]; then echo $2 bwa index found
    else echo indexing $2
    bwa index $2
fi

SAMsp1=$(echo "$sp2"_flanks-$sp1.sam)
SAMsp2=$(echo "$sp1"_flanks-$sp2.sam)

if [ -f "$SAMsp1" ]; then echo $1 sam file found
    else echo mapping $sp2 flanks on $sp1 genome
    bwa mem -t $THREADS $1 "$sp2"_flanks_300.bnk > "$sp2"_flanks-$sp1.sam
    sam2bed < "$sp2"_flanks-$sp1.sam | awk '{if (length($12)>99) print $1,$2,$3,$4,"MAPQ="$5",LENGTH="length($12)","$16,$6}' |
     awk '{sub(/:i:/,"=",$5)}1' OFS="\t" | bedtools sort -i > "$sp2"_flanks-$sp1.bed
fi

if [ -f "$SAMsp2" ]; then echo $2 sam file found
    else echo mapping $sp1 flanks on $sp2 genome
    bwa mem -t $THREADS $2 "$sp1"_flanks_300.bnk > "$sp1"_flanks-$sp2.sam
    sam2bed < "$sp1"_flanks-$sp2.sam | awk '{if (length($12)>99) print $1,$2,$3,$4,"MAPQ="$5",LENGTH="length($12)","$16,$6}' |
     awk '{sub(/:i:/,"=",$5)}1' OFS="\t" | bedtools sort -i > "$sp1"_flanks-$sp2.bed
    echo mapping complete
fi

slop=$((SINElength+300))

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

# Doubles and multies are processed in batches of 100 clusters. Each batch runs
# in its own scratch directory (node-local if SCRATCH_DIR, SLURM_TMPDIR or TMPDIR
# is set), and its alignments are appended to one bundle per category, so the
# shared work directory never holds one file per locus:
#   aln_<sp1>-<sp2>_PM.aln.gz  _MP.aln.gz  _SINE.aln.gz  _rejected.aln.gz
# A bundle is plain text: "##FILE <name>" followed by that alignment (FASTA).
# Use aln_bundle.sh to list/extract, or open them in sine_loci_browser.html.

WD=$(pwd)
SINEfa=$(readlink -f "$3")
SCRATCH=$(mktemp -d "${SCRATCH_DIR:-${SLURM_TMPDIR:-${TMPDIR:-$WD}}}/SINE_orth_loc_${sp1}-${sp2}.XXXXXX") || exit 1
trap 'rm -rf "$SCRATCH"' EXIT
BUNDLE="$WD/aln_${sp1}-${sp2}"
for t in PM MP SINE rejected; do rm -f "${BUNDLE}_$t.aln.gz"; echo "##SINE_orth_loc bundle v1" > "${BUNDLE}_$t.aln"; done
rm -f statcoords statpairs
cut -f1 "$1.fai" > "$WD/$sp1.chroms"

# each parallel ComPair.sh job writes its own *.stat file; merge them into "stat"
# serially so concurrent appends cannot interleave (unsafe on NFS/Lustre)
collect_stats() {
    find . -maxdepth 1 -type f -name "*.stat" -exec cat {} + >> "$WD/stat"
    find . -maxdepth 1 -type f -name "*.stat" -delete
}

# coordinates of both loci in every validated alignment of the current batch
# (statcoords: "cluster locus" per sequence; statpairs: alignment, seq1 locus, seq2 locus)
collect_coords() {
    for i in *.cl.*.PM *.cl.*.MP *.cl.*.SINE; do [ -f "$i" ] || continue
     awk -F:: 'NR==1 || NR==3 {c=FILENAME; sub(/\..*/,"",c); h=$1; sub(/^>/,"",h); print c, h}' "$i" >> "$WD/statcoords"
     awk -F:: 'NR==1 {h1=$1} NR==3 {h2=$1} END {sub(/^>/,"",h1); sub(/^>/,"",h2); print FILENAME"\t"h1"\t"h2}' "$i" >> "$WD/statpairs"
    done
}

# append every alignment left in the current batch directory to the bundles
bundle_batch() {
    find . -maxdepth 1 -type f ! -name "*.fai" ! -name "*.stats" -print0 |
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
       awk 'NR==1 {a=$0} NR!=1 {if ($0 ~ ">") {print $0"::"a} else print}' | cat - "$SINEfa" > "$cName".cl
      Lines=$(( $Lines - 1 ))
    done
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

echo working with doubles
echo splitting doubles into parts of 100

: > stat
awk 'NR%100==1{x="PART"++i;}{print > x}' $sp1-"$sp2"_double

echo analyzing parts

for i in PART*; do [ -f "$i" ] || continue
    run_batch "$i" align_doubles
done

mv stat stat_doubles_"$sp1"-"$sp2"

echo "stage3 - aligning doubles - completed"

echo "working with multies"

: > stat
awk 'NR%100==1{x="MULTIPART"++i;}{print > x}' $sp1-"$sp2"_multi

echo aligning clusters

for i in MULTIPART*; do [ -f "$i" ] || continue
    run_batch "$i" align_multies
done

mv stat stat_multi_"$sp1"-"$sp2"

# Multi-copy clusters (3-10 loci, of which the multi stage keeps at most one pair, and
# clusters of more than 10 loci) are usually due to a repetitive left flank. For every SINE
# copy in them that is not yet in a final alignment, map left and right flanks with
# bwa mem -a and accept a target only where both flanks agree on one placement
# (rescue_multi.py); accepted pairs are aligned and checked like doubles (clusters M<n>R).
: > stat
if [[ "$RESCUE" == 1 ]]; then
    echo "rescuing unresolved multi-copy clusters with both flanks"
    touch "$sp1-$sp2"_multi "$sp1-$sp2"_poly
    "$PYTHON" "$RESCUE_PY" prepare --clusters "$sp1-$sp2"_multi "$sp1-$sp2"_poly --resolved statpairs \
     --copies "$sp1=$sp1-${SINEname}_uniq.bed" "$sp2=$sp2-${SINEname}_uniq.bed" --genome "$sp1=$1" "$sp2=$2" --sine-length $SINElength --out rescue
    [ -s rescue_"$sp1".fa ] && bwa mem -a -t $THREADS $2 rescue_"$sp1".fa > rescue_"$sp1".sam
    [ -s rescue_"$sp2".fa ] && bwa mem -a -t $THREADS $1 rescue_"$sp2".fa > rescue_"$sp2".sam
    "$PYTHON" "$RESCUE_PY" resolve --prep rescue --sam "$sp1=rescue_$sp1.sam" "$sp2=rescue_$sp2.sam" \
     --genome "$sp1=$1" "$sp2=$2" --sine-length $SINElength --out "$sp1-$sp2"_rescue --report rescue_"$sp1-$sp2".tsv
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
echo "writing orth_$sp1-$sp2.tsv"
awk -v sp1=$sp1 -v sp2=$sp2 -F'\t' '
 FILENAME==ARGV[1] {g1[$1]; next}
 FILENAME ~ /(^|\/)stat_/ {split($0,f," "); b=f[1]; sub(/^\.\//,"",b); m=""
     for (i=3;i<=length(f);i++) m=m"\t"substr(f[i],index(f[i],"=")+1); met[b"."f[2]]=m; next}
 ($1 in met) {st=$1; sub(/.*\./,"",st); cl=$1; sub(/\..*/,"",cl)
     c1=$2; sub(/:[0-9]+-[0-9]+\([+-]\)$/,"",c1); c2=$3; sub(/:[0-9]+-[0-9]+\([+-]\)$/,"",c2)
     print $1,cl,st,((c1 in g1)?sp1:sp2),$2,((st=="PM"||st=="SINE")?1:0),((c2 in g1)?sp1:sp2),$3,((st=="MP"||st=="SINE")?1:0) met[$1]}
 BEGIN {OFS="\t"; print "alignment","cluster","status","species1","locus1","sine1","species2","locus2","sine2","LF","FL","LFcp","RFlength","FR","RFcp","OneTwo","OneSINE","TwoSINE","SL","SO","ST"}' \
 "$1.fai" stat_doubles_"$sp1"-"$sp2" stat_multi_"$sp1"-"$sp2" stat_rescue_"$sp1"-"$sp2" statpairs > orth_"$sp1"-"$sp2".tsv

# what became of every cluster: resolved (and how), rejected by ComPair.sh, or left unresolved
echo "writing clusters_$sp1-$sp2.tsv"
"$PYTHON" "$RESCUE_PY" account --double "$sp1-$sp2"_double --multi "$sp1-$sp2"_multi --poly "$sp1-$sp2"_poly \
 --stat stat_doubles_"$sp1"-"$sp2" stat_multi_"$sp1"-"$sp2" stat_rescue_"$sp1"-"$sp2" --statpairs statpairs \
 --rescue-report rescue_"$sp1-$sp2".tsv --out clusters_"$sp1"-"$sp2".tsv
