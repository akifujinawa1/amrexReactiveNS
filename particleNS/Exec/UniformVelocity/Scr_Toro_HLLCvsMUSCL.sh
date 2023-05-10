#!/bin/bash

# Run this to test Toro test 1 with HLLC vs MUSCL-Hancock. To demonstrate
# better spatial resolution via linear boundary reconstruction.

make DIM=1 USE_MPI=FALSE -j6


./main1d.gnu.ex inputs amr.n_cell = 64 64 64 slopelimiting = 0 \
enIC = 1 geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 \
adv.cfl = 0.9 amr.max_level = 0

./main1d.gnu.ex inputs amr.n_cell = 64 64 64 slopelimiting = 1 \
enIC = 1 geometry.prob_lo = 0.0 0.0 0.0 geometry.prob_hi = 1.0 1.0 1.0 \
adv.cfl = 0.9 amr.max_level = 0

gnuplot output/scripts/plot1D_HLLCvsMUSCL.p
epstopdf output/plot1.eps
rm output/plot1.eps
mv output/plot1.pdf output/pdf/toro/HLLCMUSCL.pdf