#!/bin/bash

rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*

make DIM=1 USE_MPI=TRUE -j6

Tp=( 300 )
tests=( 15 )


python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 24 ./main1d.gnu.MPI.ex inputs amr.n_cell = 768 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00768 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.040\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 900
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/domain768/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/domain768/

rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*

sed -i 's/Lx = 0.00768/Lx = 0.01024/g' setupScripts/setup1D_concClosed.py
python3 setupScripts/setup1D_concClosed.py


for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 32 ./main1d.gnu.MPI.ex inputs amr.n_cell = 1024 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.01024 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.040\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 900
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/domain1024/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/domain1024/

rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*




