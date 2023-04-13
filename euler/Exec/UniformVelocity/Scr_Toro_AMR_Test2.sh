#!/bin/bash

# Run this script after re-compiling the program with Tagging_ndEuler_test2.f90 instead of 
# Tagging_ndEuler.f90 in 'euler/Source/Src_nd/Make.package' and 'euler/CMakeLists.txt' 

make DIM=1 USE_MPI=FALSE -j6
make DIM=2 USE_MPI=FALSE -j6


cells=( 64 )
tests=( 2 )
for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    ./main1d.gnu.ex inputs amr.n_cell = "$i" "$i" "$i" slopelimiting = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 \
        adv.cfl = 0.9 amr.max_level = 1 amr.probin_file = prbn/"$j"
    done
done

gnuplot output/scripts/plot1D_MUSCL2amr.p
epstopdf output/plot2.eps
rm output/plot2.eps
mv output/plot2.pdf output/pdf/toro/toro2amr.pdf

# Must re-compile with DIM=2 to run the following 2D test

cells=( 64 )
tests=( 2 )
for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    ./main2d.gnu.ex inputs amr.n_cell = "$i" "$i" "$i" slopelimiting = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 \
        adv.cfl = 0.9 amr.max_level = 1 amr.probin_file = prbn/"$j"
    done
done

gnuplot output/scripts/plot1D_MUSCL2amr.p
epstopdf output/plot2.eps
rm output/plot2.eps
mv output/plot2.pdf output/pdf/toro/toro2amr_slice.pdf


# Run test 2 for the convergence study with AMR

cells=( 64 128 256 512 1024 )
tests=( 2 )

for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    ./main1d.gnu.ex inputs amr.n_cell = "$i" "$i" "$i" slopelimiting = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 \
        adv.cfl = 0.9 amr.max_level = 1 amr.probin_file = prbn/"$j" timing = 0 conv = 1
    done
done
