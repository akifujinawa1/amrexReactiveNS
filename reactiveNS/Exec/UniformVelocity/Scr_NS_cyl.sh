#!/bin/bash


make DIM=2 USE_MPI=FALSE -j6

./main2d.gnu.ex inputs amr.n_cell = 128 128 128 slopelimiting = 1 \
enIC = 6 geometry.prob_lo = -1.0 -1.0 -1.0 geometry.prob_hi = 1.0 1.0 1.0 \
adv.cfl = 0.9 amr.max_level = 0 amr.probin_file = prbn/1 amr.plot_files_output = 1

cd visitEuler
ls -1 plt*/Header | tee euler.visit
cd ..
mkdir output/visit_files/visit_6_64
mv visitEuler output/visit_files/visit_6_64