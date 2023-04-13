#!/bin/bash

# Run this script after connecting to a 16 core machine via ssh.
# This runs the timing exercise to test the speed-up via use of MPI and AMR.

make DIM=2 USE_MPI=TRUE -j6

cells=( 512 1024 2048 )
ncore=( 1 2 4 8 16 )

echo "Simulation time with different resolutions, with and without MPI recorded here " > output/txt/timing.txt

for i in "${cells[@]}"
do
    echo "$i " >> output/txt/timing.txt
    for j in "${ncore[@]}"
    do
        echo "$j " >> output/txt/timing.txt
	    ./main2d.gnu.MPI.ex inputs amr.n_cell = "$i" "$i" "$i" slopelimiting = 1\
        enIC = 6 geometry.prob_lo = -1.0 -1.0 -1.0 geometry.prob_hi = 1.0 1.0 1.0 \
        adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1 timing = 1
    done
done

echo "Simulation time with the same effective resolution as 2048x2048, with one or two levels of AMR, with and without MPI here " >> output/txt/timing.txt

ncore=( 1 2 4 8 16 )

for j in "${ncore[@]}"
do
    echo "$j " >> output/txt/timing.txt
    echo "512x512x2x2" >> output/txt/timing.txt
    ./main2d.gnu.MPI.ex inputs amr.n_cell = 512 512 512 slopelimiting = 1\
    enIC = 6 geometry.prob_lo = -1.0 -1.0 -1.0 geometry.prob_hi = 1.0 1.0 1.0 \
    adv.cfl = 0.9 amr.max_level = 2 amr.probin_file = prbn/timing timing = 1

    echo "1048x1048x2" >> output/txt/timing.txt
    ./main2d.gnu.MPI.ex inputs amr.n_cell = 1048 1048 1048 slopelimiting = 1\
    enIC = 6 geometry.prob_lo = -1.0 -1.0 -1.0 geometry.prob_hi = 1.0 1.0 1.0 \
    adv.cfl = 0.9 amr.max_level = 1 amr.probin_file = prbn/timing timing = 1
done

