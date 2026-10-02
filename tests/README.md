# Simulation tests

    mkdir sim && cd sim
    export PATH="$(cd .. && pwd):$PATH"         # ComPair.sh, rescue_multi.py
    python3 ../tests/simulate_genomes.py          # 4 genomes, truth tables
    for p in aaa-bbb aaa-ccc aaa-ddd bbb-ccc bbb-ddd ccc-ddd; do
        a=${p%-*}; b=${p#*-}; mkdir run_$p; cp $a.bnk $b.bnk $a-SINEX.bed $b-SINEX.bed SINEX.fa run_$p/
        (cd run_$p && bash ../../SINE_orth_loc.bash $a.bnk $b.bnk SINEX.fa > log 2>&1)
        python3 ../tests/score_pairs.py . run_$p $a $b
    done

`score_pairs.py` reports, per pair, rows whose two loci are the same insertion (correct),
different insertions (wrong), a duplicated copy paired with its counterpart (paralog),
and recall overall and for loci with a repetitive left flank. Run with `RESCUE=0` to
compare with the behaviour without the two-flank rescue.

## Genotyping from reads

    python3 ../sine_registry.py build run_aaa-bbb/orth_aaa-bbb.tsv run_aaa-ccc/orth_aaa-ccc.tsv \
        run_bbb-ccc/orth_bbb-ccc.tsv --species aaa,bbb,ccc \
        --copies aaa=aaa-SINEX.bed --copies bbb=bbb-SINEX.bed --copies ccc=ccc-SINEX.bed -o pan3
    samtools faidx aaa.bnk; samtools faidx bbb.bnk; samtools faidx ccc.bnk
    python3 ../sine_genotype.py templates --groups pan3.groups.tsv \
        --genome aaa=aaa.bnk bbb=bbb.bnk ccc=ccc.bnk -o tpl && bwa index tpl.fa
    python3 ../tests/simulate_reads.py --hap aaa.bnk --hap ccc.bnk --cov 20 --type short -o hyb.fq
    bwa mem tpl.fa hyb.fq > hyb.sam
    python3 ../sine_genotype.py call --templates tpl.tsv --sam hyb.sam -o hyb.gt.tsv
    python3 ../tests/score_genotypes.py --min-gq=20 . pan3.groups.tsv hyb.gt.tsv aaa ccc
