#!/bin/bash

# Run this script to test all of Toro's tests against the MUSCL-Hancock code. Toro test 1 is repeated with a
# domain length of sqrt(2) and discontinuity at sqrt(2)/2 at the end to calculate the exact solution, 
# which can be used later on for 2D code validation.

make DIM=1 USE_MPI=FALSE -j6

cells=( 32 64 128 )
tests=( 1 2 3 4 5 )

for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    ./main1d.gnu.ex inputs amr.n_cell = "$i" "$i" "$i" slopelimiting = 0\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 \
        adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1
    done
done

gnuplot output/scripts/plot1D_MUSCL1.p
gnuplot output/scripts/plot1D_MUSCL2.p
gnuplot output/scripts/plot1D_MUSCL3.p
gnuplot output/scripts/plot1D_MUSCL4.p
gnuplot output/scripts/plot1D_MUSCL5.p

epstopdf output/plot1.eps
epstopdf output/plot2.eps
epstopdf output/plot3.eps
epstopdf output/plot4.eps
epstopdf output/plot5.eps

rm output/plot1.eps
rm output/plot2.eps
rm output/plot3.eps
rm output/plot4.eps
rm output/plot5.eps

mv output/plot1.pdf output/pdf/toro/toro1.pdf
mv output/plot2.pdf output/pdf/toro/toro2.pdf
mv output/plot3.pdf output/pdf/toro/toro3.pdf
mv output/plot4.pdf output/pdf/toro/toro4.pdf
mv output/plot5.pdf output/pdf/toro/toro5.pdf

./main1d.gnu.ex inputs enIC = 7 amr.n_cell = 32 32 32 \
geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.41421356237 1.41421356237 1.41421356237 \
adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1