// The main.cpp file used to run the AMReX simulation. Some parts of the code is added to
// be able to read certain parameters from the inputs file. -2023W2

#include <new>
#include <iostream>
#include <iomanip>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include "eulerFunc.H"
#include "AmrLevelAdv.H"
// #include "particleContainer.H"

using namespace amrex;

amrex::LevelBld* getLevelBld ();

Real stop_time;
int enIC;             // choice of initial condition
int euler;            // =0 if no advection, =1 HLLC, =2 if MUSCL
int viscous;          // =0 if no diffusion, =1 if diffusion via central diff.
int source;           // =0 if no source terms, =1 if subcycled via rk4
int enLimiter;        // choice of limiter
int gCells;           // # of ghost cells
int timing;           // whether to run timing exercise or not
int conv;             // whether to run convergence study or not
int Da;               // Damkohler number, ratio of chemical to flow timescale
int particle;         // =0 if no particles, =1 if enabled
int Nsub;
int printRate;
int model;
double TpInitial;
double TgInitial;
double conc;
// double dp0;           // initial particle size


int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    Real dRunTime1 = amrex::second();

    int  max_step;
    Real strt_time;
    

    {
        ParmParse pp;

        max_step  = -1;
        strt_time =  0.0;
        stop_time = -1.0;

        pp.query("max_step",max_step);
        pp.query("strt_time",strt_time);
        pp.query("enIC",enIC);
        pp.query("euler",euler);
        pp.query("viscous",viscous);
        pp.query("source",source);
        pp.query("enLimiter",enLimiter);
        pp.query("timing",timing);
        pp.query("conv",conv);
        pp.query("Da",Da);
        pp.query("particle",particle);
        // pp.query("dp0",dp0);
        pp.query("TpInitial",TpInitial);
        pp.query("TgInitial",TgInitial);
        pp.query("stop_time",stop_time);
        pp.query("Nsub",Nsub);
        pp.query("conc",conc);
        pp.query("printRate",printRate);
        pp.query("model",model);

        std::cout << "end time is: " << stop_time << std::endl;
    }

    if (euler == 2){
        gCells = 2;
    }
    else {
        gCells = 1;
    }
    
    if (euler < 0 || euler > 3) {
        amrex::Abort("MUST SPECIFY a valid boolean value"); 
    } 
    
    // if (enLimiter < 1 || enLimiter > 4) {
    //     amrex::Abort("MUST SPECIFY a valid limiter value"); 
    // }

    if (strt_time < 0.0) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0) {
	    amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
        Amr amr(getLevelBld());

	amr.init(strt_time,stop_time);



	while ( amr.okToContinue() &&
  	       (amr.levelSteps(0) < max_step || max_step < 0) &&
	       (amr.cumTime() < stop_time || stop_time < 0.0) )

	{
	    //
	    // Do a coarse timestep.  Recursively calls timeStep()
	    //
	    amr.coarseTimeStep(stop_time);
	}

	// Write final checkpoint and plotfile
	if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
	    amr.checkPoint();
	}

	if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
	    amr.writePlotFile();
	}

    }

    Real dRunTime2 = amrex::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;

    std::ofstream timingdata;
    if (timing == 1)
    {
        timingdata.open("output/txt/timing.txt",std::ofstream::app);
        timingdata << dRunTime2 << std::endl;
        timingdata.close();
    }
    

    amrex::Finalize();

    
    return 0;
}
