#!/bin/bash



rm output/txt/flame/isobaric/*iteration*

make DIM=2 USE_MPI=TRUE -j6

Tp=( 300 )
tests=( 14 )

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 1 ./main2d.gnu.DEBUG.MPI.ex inputs amr.n_cell = 256 64 0 euler = 1 viscous = 0 source = 0 particle = 0 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00256 0.00064 0.0 adv.cfl = 0.5 adv.fourier = 0.5\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.1280\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010
    done
done

rm -r detonation/*plt*

# python3 output/pyscripts/particleCombustion.py

# gnuplot output/scripts/plot1D_MUSCL1NS_test.p

# epstopdf output/plot1.eps

# rm output/plot1.eps

# mv output/plot1.pdf output/pdf/toroNS/toro1.pdf

# xdg-open output/pdf/toroNS/toro1.pdf


