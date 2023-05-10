#!/bin/bash

# Run this script to test all of Toro's tests against the MUSCL-Hancock code. Toro test 1 is repeated with a
# domain length of sqrt(2) and discontinuity at sqrt(2)/2 at the end to calculate the exact solution, 
# which can be used later on for 2D code validation.

rm output/txt/test9/*

make DIM=1 USE_MPI=TRUE -j6

cells=( 128 )
tests=( 9 )

for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 6 ./main1d.gnu.MPI.ex inputs amr.n_cell = "$i" "$i" "$i" euler = 2 viscous = 0 source = 0 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 adv.cfl = 0.9 \
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 30\
        amr.plot_files_output = 1 amr.plot_file = detonation/plt
    done
done

gnuplot output/scripts/plot1D_MUSCL1NS_test.p

epstopdf output/plot1.eps

rm output/plot1.eps

mv output/plot1.pdf output/pdf/toroNS/toro1.pdf

xdg-open output/pdf/toroNS/toro1.pdf
