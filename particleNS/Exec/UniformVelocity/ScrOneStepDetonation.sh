#!/bin/bash

make DIM=1 USE_MPI=TRUE -j6

rm -r detonation/*
rm output/txt/test8/*

cells=( 512 )
tests=( 8 )

for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 2 ./main1d.gnu.MPI.ex inputs amr.n_cell = "$i" "$i" "$i" euler = 1 viscous = 0 source = 0 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 adv.cfl = 0.5 \
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 30\
        amr.plot_files_output = 1 amr.plot_file = detonation/plt 
    done
done

gnuplot output/scripts/plot1D_oneStepDet.p
epstopdf output/plot1.eps
rm output/plot1.eps
mv output/plot1.pdf output/pdf/detonation/onestep.pdf
xdg-open output/pdf/detonation/onestep.pdf