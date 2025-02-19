#!/bin/bash



rm output/txt/1Dflame/isobaric/field/*.*
rm output/txt/1Dflame/isobaric/particle/*.*

make DIM=1 USE_MPI=TRUE -j6

Tp=( 300 )
tests=( 7 )

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 1 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 0 viscous = 1 source = 0 particle = 0 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 0.0 0.0 adv.cfl = 0.8 adv.fourier = 1.0\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.010\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1
    done
done

rm -r detonation/*plt*

# python3 output/pyscripts/flameIsobaric1D.py



