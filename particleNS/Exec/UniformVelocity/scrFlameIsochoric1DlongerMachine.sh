#!/bin/bash


make DIM=1 USE_MPI=TRUE -j6

Tp=( 300 )
tests=( 15 )


rm -r /local/data/public/af793/1DflameConfined/fft512/*

mkdir /local/data/public/af793/1DflameConfined/fft512/field
mkdir /local/data/public/af793/1DflameConfined/fft512/particle

python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.050\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printRate = 1 model = 2
    done
done

sed -i 's/Lx = 0.00512/Lx = 0.00768/g' setupScripts/setup1D_concClosed.py

rm -r /local/data/public/af793/1DflameConfined/fft768/*

mkdir /local/data/public/af793/1DflameConfined/fft768/field
mkdir /local/data/public/af793/1DflameConfined/fft768/particle

python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 24 ./main1d.gnu.MPI.ex inputs amr.n_cell = 768 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00768 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.070\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printRate = 1 model = 2
    done
done


rm -r /local/data/public/af793/1DflameConfined/fft1024/*

mkdir /local/data/public/af793/1DflameConfined/fft1024/field
mkdir /local/data/public/af793/1DflameConfined/fft1024/particle


sed -i 's/Lx = 0.00768/Lx = 0.01024/g' setupScripts/setup1D_concClosed.py
python3 setupScripts/setup1D_concClosed.py


for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 32 ./main1d.gnu.MPI.ex inputs amr.n_cell = 1024 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.01024 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.090\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printRate = 1 model = 2
    done
done



