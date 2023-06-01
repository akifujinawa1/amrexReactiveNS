#!/bin/bash



rm output/txt/particle/*time*

make DIM=1 USE_MPI=TRUE -j6

cells=( 64 )
tests=( 10 )

for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 1 ./main1d.gnu.DEBUG.MPI.ex inputs amr.n_cell = "$i" "$i" "$i" euler = 0 viscous = 0 source = 0 particle = 1 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 adv.cfl = 0.9 \
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 30\
        amr.max_grid_size = 64 amr.plot_files_output = 1 amr.plot_file = detonation/plt stop_time = 100
    done
done


# python3 output/pyscripts/multiGasDiffusion.py

# gnuplot output/scripts/plot1D_MUSCL1NS_test.p

# epstopdf output/plot1.eps

# rm output/plot1.eps

# mv output/plot1.pdf output/pdf/toroNS/toro1.pdf

# xdg-open output/pdf/toroNS/toro1.pdf


