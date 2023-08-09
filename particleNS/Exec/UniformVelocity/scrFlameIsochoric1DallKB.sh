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
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.060\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printrate = 0 model = 2
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/600kb/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/600kb/

rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*

sed -i 's/conc = 600/conc = 700/g' setupScripts/setup1D_concClosed.py
python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.060\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printrate = 0 model = 2
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/700kb/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/700kb/

rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*

sed -i 's/conc = 700/conc = 800/g' setupScripts/setup1D_concClosed.py
python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.060\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printrate = 0 model = 2
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/800kb/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/800kb/
rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*

sed -i 's/conc = 800/conc = 900/g' setupScripts/setup1D_concClosed.py
python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.060\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printrate = 0 model = 2
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/900kb/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/900kb/
rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*

sed -i 's/conc = 900/conc = 1000/g' setupScripts/setup1D_concClosed.py
python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.060\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printrate = 0 model = 2
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/1000kb/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/1000kb/
rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*

sed -i 's/conc = 1000/conc = 1100/g' setupScripts/setup1D_concClosed.py
python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.060\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printrate = 0 model = 2
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/1100kb/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/1100kb/
rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*

sed -i 's/conc = 1100/conc = 1200/g' setupScripts/setup1D_concClosed.py
python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.060\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printrate = 0 model = 2
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/1200kb/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/1200kb/
rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*

sed -i 's/conc = 1200/conc = 1300/g' setupScripts/setup1D_concClosed.py
python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.060\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printrate = 0 model = 2
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/1300kb/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/1300kb/
rm output/txt/1DflameConfined/isochoric/field/*.*
rm output/txt/1DflameConfined/isochoric/particle/*.*

sed -i 's/conc = 1300/conc = 1400/g' setupScripts/setup1D_concClosed.py
python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.060\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865 printrate = 0 model = 2
    done
done

cp -r output/txt/1DflameConfined/isochoric/field/ output/txt/1DflameConfined/1400kb/
cp -r output/txt/1DflameConfined/isochoric/particle/ output/txt/1DflameConfined/1400kb/

