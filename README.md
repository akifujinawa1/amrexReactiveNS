# MUSCL-Hancock code for the 2-D Euler equations in the AMReX framework
*2023W2*

This code was written for the written assignment for the MPhil in Scientific Computing,
**"Adaptive mesh refinement for the compressible Euler equations using AMReX"**.

---

The code is adapted from the linear advection [AMReX](https://doi.org/10.5281/zenodo.2555438) code of Dr. Stephen Millmore of the Laboratory for Scientific Computing. The code runs using the AMReX source code, which should be placed in the same directory as the 'amrexEuler' folder to run.
The code is structured in two branches:

## Exec
contains:
1. UniformVelocity: A folder where the "make" command is run. The parameters for the make file are set in GNUmakefile, or can be modified in the command line when calling make.

Within UniformVelocity, we have:
- "inputs". This can be freely changed to run different simulations without having to re-compile. Relevant choices are made in each "Scr_..." script that reproduces various plots shown in the report.
- The various scripts that automate the reproduction of plots in the report are included in this folder. Click through them to see which tests will be run, the scripts all start with "Scr_" and the titles may make the contents
self-explanatory.
- The folder "prbn" contains all of the "probin" files used for the AMR tests, containing different tagging criteria, based on the test problem considered.
- "output" contains folders where the plots are stored (pdf), scripts used to run gnuplot (scripts), and text files where the solution output goes (txt).

## Source
contains:
1. "AmrLevelAdv.cpp", where most of the work has been done in transforming the original linear advection code written by Dr. Millmore to a code that can solve the two-dimensional Euler equations. The main sections of the code that have been changed/added have been commented with '-2023W2'.
2. "eulerFunc.cpp", where most functions used to solve the augmented Euler equations with the finite volume solver (the HLLC approximate Riemann solver) are included. This file also contains information on the initial conditions set (via the int enIC), and various functions to convert between primitive and conservative variables. The entirety of this file has been written by 2023W2.
3. "exactFunc.cpp", where functions used to compute the exact Riemann problem solution are included. The entirety of this file has been written by 2023W2.
4. "reconstruct.cpp", where functions used to extend the HLLC approximate Riemann solver to second order via linear boundary reconstruction are included. The entirety of this file has been written by 2023W2.
5. "main.cpp", where the program initializes and finalizes, and various function calls are made in between. The original version of this file has been written by Dr. Millmore, and has been modified/edited to allow for more variables to be parsed via the settings file, and used in "AmrLevelAdv.cpp". The main sections of the code that have been changed/added have been commented with '-2023W2'.
6. The folder "Src_nd" contains files relevant for AMR. The files "Tagging_ndEuler.f90" and "Tagging_ndEuler_test2.f90" have been adapted from the file "Tagging_nd.f90" written by Dr. Millmore, to account for the 4 conserved variables in the two-dimensional compressible Euler equations. Again, modified sections of the code are marked by '-2023W2'.

---

**To run the code,** simply run the "Scr_..." files. You will need to change the access permissions to these files - this can be done by running "chmod 755 Src_... .sh" on the command line. These files include the "make" statement relevant for the tests conducted, and most of the gnuplot scripts and conversion of plots to pdf format. The "epstopdf" tool is used to convert eps files to pdf. This may need to be installed to run those commands. The only plots that require additional post-processing is in "Scr_Toro_conv", where the convergence study is run. This was done on MATLAB - the script has been included in the submission.

The files are as follows:

- Scr_Toro_HLLCvsMUSCL.sh: Runs test 1 from ["Riemann Solvers and Numerical Methods for Fluid Dynamics" by Eleuterio Toro](https://link.springer.com/book/10.1007/b79761) with one grid resolution, using the HLLC method with cell centred values, and the MUSCL-Hancock method, and plots the two solutions against the exact Riemann solution. This is to demonstrate the improved spatial resolution via linear boundary reconstruction.
- Scr_Toro_1DMUSCL.sh: Runs the tests 1-5 from ["Riemann Solvers and Numerical Methods for Fluid Dynamics" by Eleuterio Toro](https://link.springer.com/book/10.1007/b79761) with varying grid resolutions using MUSCL-Hancock, and plots all of the numerical results against the exact Riemann solution. 
- Scr_Toro_1DMUSCL_AMR.sh: Runs the same tests 1-5 with one base grid resolution using MUSCL-Hancock, with one level of AMR, and plots all of the numerical results against the exact Riemann solution. 
- Scr_Toro_2DMUSCL_AMR.sh: Runs the same tests 1-5 with one base grid resolution using MUSCL-Hancock, now in 2D, with one level of AMR, and plots all of the numerical results against the exact Riemann solution (1D slice of 2D domain).
- Scr_Toro_conv.sh: Runs the same tests 1-5 with varying grid resolutions, and outputs all files (numerical solution and exact Riemann solution) to one folder (output/txt/conv). To calculate the L1 error for each case, some post-processing is needed - MATLAB was used by 2023W2, and the relevant script has been included in the submission.
- Scr_Diag_CylExp_2DMUSCL_visit.sh: Runs the diagonal discontinuity problem (adapted from Toro test 1) and the cylindrical explosion test (also from ["Riemann Solvers and Numerical Methods for Fluid Dynamics" by Eleuterio Toro](https://link.springer.com/book/10.1007/b79761)), for varying grid resolutions, with and without AMR. Each time a test is run, the output files are compiled to a .visit file, and are stored in 'output/visit_files'. The user can then use the [VisIt software](https://visit.llnl.gov) to read these files, and either plot the 2D contour plots of the various distributions of conserved variables, or make 1D slices of the domain using the *curve* tool, and create data files that can then be compared to exact Riemann solutions, or reference solutions.
- Scr_Timing.sh: Runs the cylindrical explosion test with a) a range of grid resolutions, with and without MPI, with up to 16 cores, and b) the same effective resolution as 2048x2048, but with one or two levels of refinement, with and without MPI, again with up to 16 cores. If your current machine does not have up to 16 cores, either modify the for loop in this script to suit your machine, or run on another machine via ssh. Simulation times are outputted in a single text file.
- Scr_Toro_AMR_Test2.sh: Runs the AMR tests with Toro test 2. This should be run last, after all the other tests have been completed, as you must re-compile the program using "Tagging_ndEuler_test2.f90" in 'euler/CMakeLists.txt' and 'euler/Source/Src_nd/Make.package'. If you wish to run the other tests with AMR, re-compile with "Tagging_ndEuler.f90" instead. These two files provide either the high/low density criteria or high/low density gradient criteria for tagging.

# amrexReactiveNS
