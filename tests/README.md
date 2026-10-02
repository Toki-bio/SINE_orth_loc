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
