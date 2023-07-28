#!/bin/bash

rm output/txt/1DflamecloseOpen/isobaric/field/*.*
rm output/txt/1DflamecloseOpen/isobaric/particle/*.*

make DIM=1 USE_MPI=TRUE -j6

Tp=( 300 )
tests=( 16 )

python3 setupScripts/setup1D_concClosed.py

for i in "${Tp[@]}"
do
    for j in "${tests[@]}"
    do
	    mpirun -np 16 ./main1d.gnu.MPI.ex inputs amr.n_cell = 512 0 0 euler = 2 viscous = 1 source = 0 particle = 2 Da = 1\
        enIC = "$j" geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 0.00512 0.0 0.0 adv.cfl = 0.8 adv.fourier = 0.8\
        amr.max_level = 0 amr.ref_ratio = 2 2 4 4 amr.probin_file = prbn/detonation amr.plot_int = 1\
        amr.max_grid_size = 64 amr.plot_files_output = 0 amr.plot_file = detonation/plt stop_time = 0.060\
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865
    done
done

cp -r output/txt/1DflamecloseOpen/isobaric/field/ output/txt/1DflamecloseOpen/600/
cp -r output/txt/1DflamecloseOpen/isobaric/particle/ output/txt/1DflamecloseOpen/600/

rm output/txt/1DflamecloseOpen/isobaric/field/*.*
rm output/txt/1DflamecloseOpen/isobaric/particle/*.*

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
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865
    done
done

cp -r output/txt/1DflamecloseOpen/isobaric/field/ output/txt/1DflamecloseOpen/700/
cp -r output/txt/1DflamecloseOpen/isobaric/particle/ output/txt/1DflamecloseOpen/700/

rm output/txt/1DflamecloseOpen/isobaric/field/*.*
rm output/txt/1DflamecloseOpen/isobaric/particle/*.*

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
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865
    done
done

cp -r output/txt/1DflamecloseOpen/isobaric/field/ output/txt/1DflamecloseOpen/800/
cp -r output/txt/1DflamecloseOpen/isobaric/particle/ output/txt/1DflamecloseOpen/800/
rm output/txt/1DflamecloseOpen/isobaric/field/*.*
rm output/txt/1DflamecloseOpen/isobaric/particle/*.*

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
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865
    done
done

cp -r output/txt/1DflamecloseOpen/isobaric/field/ output/txt/1DflamecloseOpen/900/
cp -r output/txt/1DflamecloseOpen/isobaric/particle/ output/txt/1DflamecloseOpen/900/
rm output/txt/1DflamecloseOpen/isobaric/field/*.*
rm output/txt/1DflamecloseOpen/isobaric/particle/*.*

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
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865
    done
done

cp -r output/txt/1DflamecloseOpen/isobaric/field/ output/txt/1DflamecloseOpen/1000/
cp -r output/txt/1DflamecloseOpen/isobaric/particle/ output/txt/1DflamecloseOpen/1000/
rm output/txt/1DflamecloseOpen/isobaric/field/*.*
rm output/txt/1DflamecloseOpen/isobaric/particle/*.*

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
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865
    done
done

cp -r output/txt/1DflamecloseOpen/isobaric/field/ output/txt/1DflamecloseOpen/1100/
cp -r output/txt/1DflamecloseOpen/isobaric/particle/ output/txt/1DflamecloseOpen/1100/
rm output/txt/1DflamecloseOpen/isobaric/field/*.*
rm output/txt/1DflamecloseOpen/isobaric/particle/*.*

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
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865
    done
done

cp -r output/txt/1DflamecloseOpen/isobaric/field/ output/txt/1DflamecloseOpen/1200/
cp -r output/txt/1DflamecloseOpen/isobaric/particle/ output/txt/1DflamecloseOpen/1200/
rm output/txt/1DflamecloseOpen/isobaric/field/*.*
rm output/txt/1DflamecloseOpen/isobaric/particle/*.*

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
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865
    done
done

cp -r output/txt/1DflamecloseOpen/isobaric/field/ output/txt/1DflamecloseOpen/1300/
cp -r output/txt/1DflamecloseOpen/isobaric/particle/ output/txt/1DflamecloseOpen/1300/
rm output/txt/1DflamecloseOpen/isobaric/field/*.*
rm output/txt/1DflamecloseOpen/isobaric/particle/*.*

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
        TpInitial = "$i" TgInitial = 300 dp0 = 0.000010 Nsub = 1 conc = 865
    done
done

cp -r output/txt/1DflamecloseOpen/isobaric/field/ output/txt/1DflamecloseOpen/1400/
cp -r output/txt/1DflamecloseOpen/isobaric/particle/ output/txt/1DflamecloseOpen/1400/

