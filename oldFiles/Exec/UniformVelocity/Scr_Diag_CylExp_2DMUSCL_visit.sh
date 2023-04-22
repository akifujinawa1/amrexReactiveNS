#!/bin/bash

# This script runs the diagonal discontinuity test (test 7) with three resolutions, and
# creates a VisIt compatible file for each run in output/visit_files (visit_7_resolution)

make DIM=2 USE_MPI=FALSE -j6

./main2d.gnu.ex inputs amr.n_cell = 64 64 64 slopelimiting = 1 \
enIC = 7 geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.41421356237 1.41421356237 1.41421356237 \
adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1 amr.plot_files_output = 1

cd visitEuler
ls -1 plt*/Header | tee euler.visit
cd ..
mkdir output/visit_files/visit_7_64
mv visitEuler output/visit_files/visit_7_64

./main2d.gnu.ex inputs amr.n_cell = 128 128 128 slopelimiting = 1 \
enIC = 7 geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.41421356237 1.41421356237 1.41421356237 \
adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1 amr.plot_files_output = 1

cd visitEuler
ls -1 plt*/Header | tee euler.visit
cd ..
mkdir output/visit_files/visit_7_128
mv visitEuler output/visit_files/visit_7_128

./main2d.gnu.ex inputs amr.n_cell = 256 256 256 slopelimiting = 1 \
enIC = 7 geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.41421356237 1.41421356237 1.41421356237 \
adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1 amr.plot_files_output = 1

cd visitEuler
ls -1 plt*/Header | tee euler.visit
cd ..
mkdir output/visit_files/visit_7_256
mv visitEuler output/visit_files/visit_7_256

# This script runs the cylindrical explosion test (test 6) with three resolutions, and
# creates a VisIt compatible file for each run in output/visit_files (visit_6_resolution)

./main2d.gnu.ex inputs amr.n_cell = 64 64 64 slopelimiting = 1 \
enIC = 6 geometry.prob_lo = -1.0 -1.0 -1.0 geometry.prob_hi = 1.0 1.0 1.0 \
adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1 amr.plot_files_output = 1

cd visitEuler
ls -1 plt*/Header | tee euler.visit
cd ..
mkdir output/visit_files/visit_6_64
mv visitEuler output/visit_files/visit_6_64

./main2d.gnu.ex inputs amr.n_cell = 128 128 128 slopelimiting = 1 \
enIC = 6 geometry.prob_lo = -1.0 -1.0 -1.0 geometry.prob_hi = 1.0 1.0 1.0 \
adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1 amr.plot_files_output = 1

cd visitEuler
ls -1 plt*/Header | tee euler.visit
cd ..
mkdir output/visit_files/visit_6_128
mv visitEuler output/visit_files/visit_6_128

./main2d.gnu.ex inputs amr.n_cell = 256 256 256 slopelimiting = 1 \
enIC = 6 geometry.prob_lo = -1.0 -1.0 -1.0 geometry.prob_hi = 1.0 1.0 1.0 \
adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1 amr.plot_files_output = 1

cd visitEuler
ls -1 plt*/Header | tee euler.visit
cd ..
mkdir output/visit_files/visit_6_256
mv visitEuler output/visit_files/visit_6_256

./main2d.gnu.ex inputs amr.n_cell = 512 512 512 slopelimiting = 1 \
enIC = 6 geometry.prob_lo = -1.0 -1.0 -1.0 geometry.prob_hi = 1.0 1.0 1.0 \
adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1 amr.plot_files_output = 1

cd visitEuler
ls -1 plt*/Header | tee euler.visit
cd ..
mkdir output/visit_files/visit_6_512
mv visitEuler output/visit_files/visit_6_512


# the following script runs the cylindrical explosion test with one AMR level

./main2d.gnu.ex inputs amr.n_cell = 128 128 128 slopelimiting = 1 \
enIC = 6 geometry.prob_lo = -1.0 -1.0 -1.0 geometry.prob_hi = 1.0 1.0 1.0 \
adv.cfl = 0.9 amr.max_level = 1 amr.probin_file = prbn/1 amr.plot_files_output = 1

cd visitEuler
ls -1 plt*/Header | tee euler.visit
cd ..
mkdir output/visit_files/visit_6_128_amr
mv visitEuler output/visit_files/visit_6_128_amr