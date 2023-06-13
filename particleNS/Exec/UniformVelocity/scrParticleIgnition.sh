#!/bin/bash



rm output/txt/particleIgnition/*data*

make DIM=1 USE_MPI=TRUE -j6

Tp=( 1270 )
# Tp=( 1076 1075 1069 1054 )
tests=( 12 )

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 1 ./main1d.gnu.DEBUG.MPI.ex inputs amr.n_cell = 64 64 64 euler = 0 viscous = 0 source = 0 particle = 1 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 adv.cfl = 0.9 \
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 100\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.00002728 stop_time = 0.04
    done
done

rm -r detonation/*plt*

# python3 output/pyscripts/particleIgnition.py

# gnuplot output/scripts/plot1D_MUSCL1NS_test.p

# epstopdf output/plot1.eps

# rm output/plot1.eps

# mv output/plot1.pdf output/pdf/toroNS/toro1.pdf

# xdg-open output/pdf/toroNS/toro1.pdf


