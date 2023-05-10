#!/bin/bash

# Run this script to test all of Toro's tests against the MUSCL-Hancock code. Toro test 1 is repeated with a
# domain length of sqrt(2) and discontinuity at sqrt(2)/2 at the end to calculate the exact solution, 
# which can be used later on for 2D code validation.

make DIM=1 USE_MPI=FALSE -j6

cells=( 128 )
tests=( 1 )

for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    ./main1d.gnu.ex inputs amr.n_cell = "$i" "$i" "$i" slopelimiting = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 \
        adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1
    done
done

gnuplot output/scripts/plot1D_MUSCL1NS_test.p

epstopdf output/plot1.eps

rm output/plot1.eps

mv output/plot1.pdf output/pdf/toro/toro1.pdf
