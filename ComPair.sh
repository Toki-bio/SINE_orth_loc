#!/bin/bash
# Analysis of 3-way alignments (plus-minus SINE loci).
# Each aligment of SINE-containing(?) loci is tested if it has 
# 1) SINE in one genome and gap in otherm or SINE in both genomes 
# 2) good left and right flanks.

# searching for consecutive gaps at both sides of SINE consensus in alignment
# find coordinates of SINE consensus relative to multiple alignment

echo $1
bank=$1
STAT=${COMPAIR_STAT:-stat} # per-job stat file when run in parallel (set by SINE_orth_loc.bash)

# thresholds (defaults = published criteria); can be tuned in sine_loci_browser.html and exported
CP_FLANK_LEN=${CP_FLANK_LEN:-150}   # min length of left/right flank (alignment columns)
CP_FLANK_NID=${CP_FLANK_NID:-65}    # min identical nt in left flank
CP_FLANK_PID=${CP_FLANK_PID:-65}    # min % identity of left and right flanks
CP_SINE_PID=${CP_SINE_PID:-65}      # SINE in both: min % identity between the two SINEs
CP_SINE_NID=${CP_SINE_NID:-100}     # SINE in both: min identical nt between SINEs and with consensus
CP_CUT_PID=${CP_CUT_PID:-70}        # SINE in one: min % identity of the trimmed right flank
CP_CUT_NID=${CP_CUT_NID:-120}       # SINE in one: min identical nt in the trimmed right flank
CP_INDEL_MIN=${CP_INDEL_MIN:-50}    # right flank: one-sided gap runs >= this are masked before re-testing (0: off)
# CP_SIZES: chrom<TAB>size of both genomes; rejected loci whose window touches a sequence end -> contig_end
HELPER="$(dirname "$(readlink -f "$0")")/compair_flank.py"
PY=${PYTHON:-python3}
FI=0                                 # 1: right flank passed only after masking a large indel
EDGE=-1                              # which locus (0/1) has its window clipped by a sequence end
at_edge() {  # $1: alignment file (default $bank)
    [[ -n "$CP_SIZES" && -s "$CP_SIZES" ]] || return 1
    EDGE=$("$PY" "$HELPER" edge "${1:-$bank}" "$CP_SIZES"); [[ "$EDGE" =~ ^[0-9]+$ ]] || { EDGE=-1; return 1; }
}
# contig_end if the left flank (orthology evidence) is good, contig_end_lf otherwise
edge_status() { if [[ ${LFcp:-0} -gt $CP_FLANK_NID && $LF -gt $CP_FLANK_LEN && ${FL:-0} -gt $CP_FLANK_PID ]]; then echo contig_end; else echo contig_end_lf; fi; }
# masked re-test of right-flank columns $1..$2; needs > $3 % identity and > $4 identical nt over >= CP_FLANK_LEN kept columns
masked_rf_ok() {
    [[ "$CP_INDEL_MIN" -gt 0 ]] || return 1
    read -r mpid mnid mcols mmasked <<< "$("$PY" "$HELPER" mask "$bank" "$1" "$2" "$CP_INDEL_MIN")"
    echo "masked right flank: pid=$mpid nid=$mnid cols=$mcols masked=$mmasked"
    [[ "$mmasked" -gt 0 && "$mcols" -ge "$CP_FLANK_LEN" && "$mpid" -gt "$3" && "$mnid" -gt "$4" ]]
}
seqkit range -r -1:-1 -w 0 $bank | awk -F '[^-]+' 'NR!=1 {for (i=1; i<=NF; i++) if ($i != "") print length($i)}' > $bank.gaps # first and last values are lengths of 3' and 5' rows of consecutive gaps in SINE sequence
LF=$(awk 'NR==1' $bank.gaps) # right coordinate of the left flank
RFlength=$(awk 'END {print}' $bank.gaps) # length of the right flank
len=$(seqkit range -r -1:-1 -w 0 $bank | awk 'NR==2 { print length($0) }') # alignment length
RF=$((len-RFlength)) # left coordinate of the right flank
echo LF=$LF len=$len RF=$RF RFlength=$RFlength
rm $bank.gaps

#report problematic flanks length
if [[ "$LF" -le "$CP_FLANK_LEN" ]]; then echo "Short FLANK! Trying to fix it"
    part=$(seqkit range -r -1:-1 -w 0 $bank | awk 'NR==2' | perl -pe 's/\p{L}([a-zA-Z]*)/"-" x (length($1)+1)/ei')
    awk -v part=$part 'NR<6 {print $0} END {print part}' $bank > $bank.temp
    mv $bank.temp $bank
    LF=$(seqkit range -r -1:-1 -w 0 $bank | awk -F '[^-]+' 'NR!=1 {print length($1)}') # recalculating new left flank coordinate
    echo new left flank $LF
fi
if [[ "$RFlength" -lt "$CP_FLANK_LEN" ]]; then echo "Short FRANK!"
    if at_edge; then
        seqkit subseq -r 1:$LF -w 0 $bank | esl-alipid - | awk 'NR==2 {print $3,$4}' > $bank.lfstat
        FL=$(awk -F. '{print $1}' $bank.lfstat); LFcp=$(awk '{print $2}' $bank.lfstat); rm -f $bank.lfstat $bank.seqkit.fai
        status=$(edge_status); mv $bank $bank.$status
        echo $bank $status LF=$LF FL=$FL LFcp=$LFcp RFlength=$RFlength EDGE=$EDGE >> "$STAT"; exit
    fi
    mv $bank $bank.shortRF; status=$(echo "shortRF"); echo $bank $status >> "$STAT"; exit; fi

leftSINE=$(seqkit range -r -1:-1 -w 0 $bank | awk '/^[^-]/{if (NR==2) print "0"}')
if [[ "$leftSINE" == 0 ]]
       then
        echo "Skipping (leftSINE)"
        mv $bank $bank.lfSINE; status=$(echo "lfSINE")
echo $bank $status LF=$LF FL=$FL LFcp=$LFcp RFlength=$RFlength FR=$FR RFcp=$RFcp OneTwo=$OneTwo OneSINE=$OneSINE TwoSINE=$TwoSINE SL=$SL SO=$SO ST=$ST FI=$FI EDGE=$EDGE >> "$STAT"
        exit 111 # aborting search in this locus
fi

rightSINE=$(seqkit range -r -1:-1 -w 0 $bank | awk '{if (NR==2) {print substr($0,length($0),1)}}')
if [[ "$rightSINE" == "-" ]]
       then
        echo "right"
        else
        mv $bank $bank.rfSINE; status=$(echo "rfSINE")
        echo "Skipping (rightSINE)"
echo $bank $status LF=$LF FL=$FL LFcp=$LFcp RFlength=$RFlength FR=$FR RFcp=$RFcp OneTwo=$OneTwo OneSINE=$OneSINE TwoSINE=$TwoSINE SL=$SL SO=$SO ST=$ST FI=$FI EDGE=$EDGE >> "$STAT"
        exit 111 # aborting search in this locus
fi

#coordinates for cutting left flank, SINE and right flank

LFcut=$(( $LF + 1 ))
RFcut=$(( $RF + 1 ))

#gathering statistics. IMPORTANT: esl-alipid cannot process sequences <9bp

#FLANK
seqkit subseq -r 1:$LF -w 0 $bank | esl-alipid - | awk 'NR==2 {print $3,$4,"FLANK"}' > $bank.stats
LFcp=$(awk 'NR==1 {print $2}' $bank.stats) # number of identical nt in left flanks
#FRANK
seqkit subseq -r $RFcut:$len -w 0 $bank | esl-alipid - | awk 'NR==2 {print $3,$4,"FRANK"}' >> $bank.stats
RFcp=$(awk 'NR==2 {print $2}' $bank.stats) # number of identical nt in right flanks
#SINE
seqkit subseq -r $LFcut:$RF -w 0 $bank | esl-alipid - | awk 'NR==2 || NR==3 || NR==4 {print $3,"\n"$4}' >> $bank.stats # percentage identity and number of identities

FL=$(awk -F. 'NR==1 {print $1}' $bank.stats) # identity in flanks
FR=$(awk -F. 'NR==2 {print $1}' $bank.stats) # identity in franks
OneTwo=$(awk -F. 'NR==3 {print $1}' $bank.stats) # identity in SINE between genomes
OneSINE=$(awk -F. 'NR==5 {print $1}' $bank.stats) # identity between first genome locus and SINE consensus
TwoSINE=$(awk -F. 'NR==7 {print $1}' $bank.stats) # identity between second genomelocus and SINE consensus
SL=$(awk -F. 'NR==4 {print $1}' $bank.stats) # number of identical bp between SINEs in two genomes
SO=$(awk -F. 'NR==6 {print $1}' $bank.stats) # number of identical bp between first genome locus and SINE consensus
ST=$(awk -F. 'NR==8 {print $1}' $bank.stats) # number of identical bp between second genome locus and SINE consensus

#check if SINE is + or - or bad, check if flank and/or frank is good

if [[ $LFcp -gt $CP_FLANK_NID && $LF -gt $CP_FLANK_LEN && $FL -gt $CP_FLANK_PID ]] # check if FLANK is good
    then echo left flank good, next check SINE and right flank
	if [[ $OneTwo -gt $CP_SINE_PID && $SL -gt $CP_SINE_NID && $SO -gt $CP_SINE_NID && $ST -gt $CP_SINE_NID ]]
	    then echo SINE in both, next check right flank
		if [[ $RFlength -gt $CP_FLANK_LEN && $FR -gt $CP_FLANK_PID ]]
		    then echo right flank is good, homozygous copy; mv $bank $bank.SINE; status=$(echo "SINE")
		elif masked_rf_ok $RFcut $len $CP_FLANK_PID 0
		    then echo right flank good after masking a large indel; FI=1; mv $bank $bank.SINE; status=$(echo "SINE")
		    else echo right flank bad, next try to fix it?; mv $bank $bank.badRF; status=$(echo "badRF")
		fi
	    else echo SINE can be polymorphic, need to check right flank, preparing it first
		seqkit subseq -r $RFcut:$len -w 0 $bank > $bank.tocut
		cutRF=$(awk '!/^>/ { sub(/-*$/, ""); print}' $bank.tocut | awk 'NR<3 {print length($0)}' | sort -n | head -n 1)   #delete overhang at the tail of minus sequence
		    if [[ $cutRF -le $CP_FLANK_LEN ]]
			then echo Right flank too short; status=shortRF; at_edge && status=$(edge_status); mv $bank $bank.$status; rm $bank.tocut* $bank*.fai $bank.stats
echo $bank $status LF=$LF FL=$FL LFcp=$LFcp RFlength=$RFlength FR=$FR RFcp=$RFcp OneTwo=$OneTwo OneSINE=$OneSINE TwoSINE=$TwoSINE SL=$SL SO=$SO ST=$ST FI=$FI EDGE=$EDGE >> "$STAT"; exit
		    fi
		seqkit subseq -r 1:$cutRF -w 0 $bank.tocut | esl-alipid - | awk 'NR==2 {print $3,$4,"FRANKcut"}' >> $bank.stats; rm $bank.tocut*
		FRcut=$(awk -F. 'NR==9 {print $1}' $bank.stats) # identity in cut franks
		FRcp=$(awk 'NR==9 {print $2}' $bank.stats) # number of identical nt in cut franks
		if [[ $cutRF -gt $CP_FLANK_LEN && $FRcut -gt $CP_CUT_PID && $FRcp -gt $CP_CUT_NID ]] || \
		   { masked_rf_ok $RFcut $((RF+cutRF)) $CP_CUT_PID $CP_CUT_NID && FI=1; }
	    	    then echo right flank is good, heterozygous copy, next find which plus/minus
			trim_right=$((RF+cutRF))
			seqkit subseq -r 1:$trim_right -w 0 $bank > $bank.trim
			if [[ $SO -gt $ST ]]
			    then
				echo plus minus
				mv $bank.trim $bank.PM; rm $bank; status=$(echo "PM")
			    else
				echo minus plus
				mv $bank.trim $bank.MP; rm $bank; status=$(echo "MP")
			fi
	    	    else echo right flank bad, next try to fix it; mv $bank $bank.badRF; status=$(echo "badRF")
		fi
	fi
    else echo left flank bad $LFcp $LF $FL; mv $bank $bank.badLF; status=$(echo "badLF")
fi
if [[ "$status" == bad* ]] && at_edge "$bank.$status"; then new=$(edge_status); mv $bank.$status $bank.$new; status=$new; fi
echo $bank $status LF=$LF FL=$FL LFcp=$LFcp RFlength=$RFlength FR=$FR RFcp=$RFcp OneTwo=$OneTwo OneSINE=$OneSINE TwoSINE=$TwoSINE SL=$SL SO=$SO ST=$ST FI=$FI EDGE=$EDGE >> "$STAT"
rm $bank.seqkit.fai $bank.stats
exit

