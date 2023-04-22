#!/bin/bash

# The 2D implementation is tested using extended 1D tests from Toro. Test 2 can be run in
# Scr_Toro_AMR_Test2.sh

make DIM=2 USE_MPI=FALSE -j6


cells=( 64 )
tests=( 1 3 4 5 )
for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    ./main2d.gnu.ex inputs amr.n_cell = "$i" "$i" "$i" slopelimiting = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 \
        adv.cfl = 0.9 amr.max_level = 1 amr.probin_file = prbn/"$j"
    done
done

gnuplot output/scripts/plot2D_MUSCL1amr.p
gnuplot output/scripts/plot2D_MUSCL3amr.p
gnuplot output/scripts/plot2D_MUSCL4amr.p
gnuplot output/scripts/plot2D_MUSCL5amr.p

epstopdf output/plot1.eps
epstopdf output/plot3.eps
epstopdf output/plot4.eps
epstopdf output/plot5.eps

rm output/plot1.eps
rm output/plot3.eps
rm output/plot4.eps
rm output/plot5.eps

mv output/plot1.pdf output/pdf/toro/toro1amr_slice.pdf
mv output/plot3.pdf output/pdf/toro/toro3amr_slice.pdf
mv output/plot4.pdf output/pdf/toro/toro4amr_slice.pdf
mv output/plot5.pdf output/pdf/toro/toro5amr_slice.pdf