#!/bin/bash

# Run this script to output files for the convergence study all into one folder.
# MATLAB was used for post-processing all of the field values to compare to the exact solution
# and compute all L1 errors - the file is included in the submission.

make DIM=1 USE_MPI=FALSE -j6


cells=( 64 128 256 512 1024 )
tests=( 1 2 3 4 5 )

for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    ./main1d.gnu.ex inputs amr.n_cell = "$i" "$i" "$i" slopelimiting = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 \
        adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1 timing = 0 conv = 1
    done
done

cells=( 64 128 256 512 1024 )
tests=( 1 3 4 5 )

for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    ./main1d.gnu.ex inputs amr.n_cell = "$i" "$i" "$i" slopelimiting = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 \
        adv.cfl = 0.9 amr.max_level = 1 amr.probin_file = prbn/"$j" timing = 0 conv = 1
    done
done



