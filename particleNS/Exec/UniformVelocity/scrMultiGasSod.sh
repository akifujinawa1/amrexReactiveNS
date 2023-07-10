#!/bin/bash



# rm output/txt/multiGasSod/*data*

make DIM=1 USE_MPI=TRUE -j6

cells=( 64 )
tests=( 9 )

for i in "${cells[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 6 ./main1d.gnu.MPI.ex inputs amr.n_cell = "$i" 0 0 euler = 2 viscous = 0 source = 0 particle = 0\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 0.0 0.0 adv.cfl = 0.9 \
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 30\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt conv = 1 stop_time = 0.00078538337
    done
done


# python3 output/pyscripts/multiGasSod.py

# gnuplot output/scripts/plot1D_MUSCL1NS_test.p

# epstopdf output/plot1.eps

# rm output/plot1.eps

# mv output/plot1.pdf output/pdf/toroNS/toro1.pdf

# xdg-open output/pdf/toroNS/toro1.pdf


