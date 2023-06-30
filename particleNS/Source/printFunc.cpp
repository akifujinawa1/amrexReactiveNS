#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>

#include <eulerFunc.H>
#include <exactFunc.H>
#include <recon.H>
#include <diffusionFunc.H>
#include <source.H>
// #include <constants.H>
#include <thermoTransport.H>
#include <particleFunc.H>
// #include <particleContainer.H>
#include <math.h>
#include <cmath>
#include <iostream>
#include <algorithm>

// import global variables defined in main.cpp, mostly derived from inputs to the settings file -2023W2
extern   Real      stop_time;
extern   int       enIC;
extern   int       euler;
extern   int       viscous;
extern   int       source;
extern   int       enLimiter;
extern   int       gCells;
extern   int       conv; 
extern   int       Da;
extern   int       particle;
extern   int       NUM_STATE;
extern   int       pfrequency;
extern   int       iter;
extern   int       printlevel;
extern   int       n_cell;
extern   int       counter;

extern   double       Gamma;           // ratio of specific heats -
extern   double       R;               // univeral gas constant   J/K/mol
extern   double       one_atm_Pa;      // one atmosphere          Pa
extern   double       T0;              // initial gas temperature K
extern   double       M_O2;            // O2 molecular weight     kg/mol
extern   double       M_N2;            // N2 molecular weight     kg/mol
extern   double       X_O2;            // mole fraction of O2
extern   double       X_N2;            // mole fraction of N2
extern   double       Y_O2;            // mass fraction of O2
extern   double       Y_N2;            // mass fraction of N2
extern   double       Mavg;            // average molecular weight of the gas mixture
extern   double       TpInitial;


void AmrLevelAdv::writePlotFile()
{
    // AmrLevel::writePlotFile (dir,os,how);
    // const std::string &dir,
    //                             std::ostream &os,
    //                             VisMF::How how
    iter+= 1;

    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();

    //   // If final time, loop through all patches here.
    const Real *dx = geom.CellSize();
    const Real *prob_lo = geom.ProbLo();
    const Real *prob_hi = geom.ProbHi();
    const Real cur_time = state[Phi_Type].curTime();

    int NxCell = int(prob_hi[0]/dx[0]);
    int NyCell = int(prob_hi[1]/dx[1]);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if ((enIC == 14)&&(cur_time > 1e-5*counter)) { // 1D isobaric flame case
        counter += 1;
        const MultiFab &S_plot = get_new_data(Phi_Type);

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        Vector<double> Tg_x(NxCell);
        Vector<double> YO2_x(NxCell);

        std::cout << "Printing data on AMR level: " << level << " with CPU rank: " << rank << std::endl;

        const Real dX = dx[0];
        const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);

        const Real probLoX = prob_lo[0];
        const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0);

        std::ofstream approx;
        std::string iteration = std::to_string(iter);
        int microTime         = std::floor(cur_time * 1.0e5);
        std::string curtime   = std::to_string(microTime);

        approx.open("output/txt/1Dflame/isobaric/field/" + curtime + "0.txt", std::ofstream::app);

        for (MFIter mfi(S_plot); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            const Dim3 lo = lbound(bx);
            const Dim3 hi = ubound(bx);
            const auto &arr = S_plot.array(mfi);
            
            int j = 0;
            int k = 0;

            for (int i = lo.x; i <= hi.x; i++){
                const Real x = probLoX + (double(i) + 0.5) * dX;
            
                double rho = arr(i, j, k, gasVar::rho);
                double xmomentum = arr(i, j, k, gasVar::rhou);
                double ymomentum = arr(i, j, k, gasVar::rhov);
                double ener = arr(i, j, k, gasVar::E);
                double rhoYO2 = arr(i, j, k, gasVar::rhoYO2);
                double rhoYN2 = arr(i, j, k, gasVar::rhoYN2);

                double vx = xmomentum / rho;
                double vy = ymomentum / rho;
                double YO2 = rhoYO2 / rho;
                double YN2 = rhoYN2 / rho;
                double Tgas = Tg(rho, vx, vy, YO2, YN2, ener);
                double p = pressure(rho, YO2, YN2, Tgas);
                double eps = specIntEner(rho, vx, vy, ener, p);
                double a = soundSpeed(p, rho, Tgas, YO2, YN2);
                double gammaval = a * a * rho / p;

                // std::cout << "rho rhou rhov e o2c n2c: " << rho << " " << xmomentum << " " << ymomentum << " " \
                // << ener << " " << rhoYO2 << " " << rhoYN2 << std::endl;

                // std::cout << "Tg, p: " << Tgas << " " << p << std::endl;

                AllPrint(approx) << x << " " << Tgas << " " << YO2 << " " << p << " " << vx << " " << rho << std::endl;
            }
        }
        approx.close();
    }
    if ((enIC == 14)&&(cur_time > 5e-6*counter)) { // 1D isobaric flame case
        const MultiFab &S_plot = get_new_data(Phi_Type);
        const int lev = 0;

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        Real x, y=0, up, vp, wp, mFe, mFeO, mFe3O4, Hp, Tp, LFe, LFeO, LFe3O4;
        int  phaseFe, phaseFeO, phaseFe3O4;

        // for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        for (MFIter mfi(S_plot); mfi.isValid(); ++mfi)
        {

            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
            auto& particles = particle_tile.GetArrayOfStructs();
            const int np = particles.numParticles();
            // std::cout << "number of particles in grid in getParticleInfo, within advance: " << np << std::endl;

            for(int pindex = 0; pindex < np; ++pindex) {
                ParticleType& p = particles[pindex];
                const IntVect& iv = this->Index(p, lev);

                std::ofstream parOut;
                std::string id        = std::to_string(p.idata(IntData::pIter));
                std::string cpu        = std::to_string(p.cpu());

                parOut.open("output/txt/1Dflame/isobaric/particle/" + id + "-" + cpu + ".txt", std::ofstream::app);

                x  = p.pos(0);   //given in meters, lo.x is in cell count, so scale by dx to get meters
                if (spacedim == 2){
                    y  = p.pos(1);
                }
                up     = p.rdata(RealData::u);
                vp     = p.rdata(RealData::v);
                wp     = p.rdata(RealData::w);
                mFe    = p.rdata(RealData::mFe);
                mFeO   = p.rdata(RealData::mFeO);
                mFe3O4 = p.rdata(RealData::mFe3O4);
                Hp     = p.rdata(RealData::Hp);
                LFe    = p.rdata(RealData::LFe);
                LFeO   = p.rdata(RealData::LFeO);
                LFe3O4 = p.rdata(RealData::LFe3O4);
                phaseFe    = p.idata(IntData::Fe);
                phaseFeO   = p.idata(IntData::FeO);
                phaseFe3O4 = p.idata(IntData::Fe3O4);
                int regime     = p.idata(IntData::regime);
                int particleID = p.idata(IntData::pIter);
                int regime = p.idata(IntData::regime);
                Tp = Tparticle(mFe,mFeO,mFe3O4,Hp,phaseFe,phaseFeO,phaseFe3O4,LFe,LFeO,LFe3O4);
                
                AllPrint(parOut) << cur_time << " " << x << " " << mFe << " " << mFeO << " " << mFe3O4 << " " << Tp << " " << regime << std::endl;
                
                parOut.close();

            }
        }
    }


    if ((enIC == 12)){ // if testing for particle ignition
        Vector<double> pPosTp(3);
        Vector<double> pReal(RealData::ncomps);
        Vector<int> pInt(IntData::ncomps);

        pPosTp = getParticleInfo(pReal,pInt);

        std::string Tp0 = std::to_string((int)TpInitial);
        
        double x = pPosTp[0];
        double y = pPosTp[1];
        double Tp = pPosTp[2];
        double mFe = pReal[RealData::mFe];
        double mFeO = pReal[RealData::mFeO];
        double mFe3O4 = pReal[RealData::mFe3O4];


        
        std::ofstream approx;
        approx.open("output/txt/particleIgnition/data"+Tp0+".txt", std::ofstream::app);

        AllPrint(approx) << cur_time << " " << Tp << std::endl;
    }
    if ((enIC == 13)){ // if testing for particle combustion
        Vector<double> pPosTp(3);
        Vector<double> pReal(RealData::ncomps);
        Vector<int> pInt(IntData::ncomps);

        pPosTp = getParticleInfo(pReal,pInt);

        std::string Tp0 = std::to_string((int)TpInitial);
        
        double x = pPosTp[0];
        double y = pPosTp[1];
        double Tp = pPosTp[2];
        double mFe = pReal[RealData::mFe];
        double mFeO = pReal[RealData::mFeO];
        double mFe3O4 = pReal[RealData::mFe3O4];


        
        std::ofstream approx;
        approx.open("output/txt/particleCombustion/data"+Tp0+".txt", std::ofstream::app);

        AllPrint(approx) << cur_time << " " << Tp << std::endl;
    }

    if ((enIC == 10)||(enIC == 11)){ // if testing for particle drag or heat
        Vector<double> pPosTp(3);
        Vector<double> pReal(RealData::ncomps);
        Vector<int> pInt(IntData::ncomps);

        pPosTp = getParticleInfo(pReal,pInt);
        
        double x = pPosTp[0];
        double y = pPosTp[1];
        double Tp = pPosTp[2];

        double u = pReal[0];
        
        std::ofstream approx;
        approx.open("output/txt/particleDrag/data.txt", std::ofstream::app);

        AllPrint(approx) << cur_time << " " << x << " " << u << " " << Tp << " " << pReal[RealData::w] << std::endl;
    }



    if (((enIC == 8) || (enIC == 9)) && (cur_time == stop_time))
    {

        const MultiFab &S_plot = get_new_data(Phi_Type);

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        std::cout << "Printing data on AMR level: " << level << " with CPU rank: " << rank << std::endl;

        const Real dX = dx[0];
        const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);

        const Real probLoX = prob_lo[0];
        const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0);

        std::string method, test;

        std::string resolution = std::to_string(n_cell);
        std::string dimension = std::to_string(amrex::SpaceDim);
        std::string iteration = std::to_string(iter);
        std::string stoptime = std::to_string((int)stop_time);

        if (enIC == 9)
        {
            test = "multiGasSod";
        }
        if (enIC == 8)
        {
            test = "multiGasDiffusion";
        }

        std::ofstream approx;

        if (conv == 1)
        {
            approx.open("output/txt/" + test + "/time" + iteration + resolution + ".txt", std::ofstream::app);
        }
        if (enIC == 8)
        {
            approx.open("output/txt/" + test + "/time" + iteration + stoptime + ".txt", std::ofstream::app);
        }
        else
        {
            approx.open("output/txt/" + test + "/time" + iteration + ".txt", std::ofstream::app);
        }

        for (MFIter mfi(S_plot); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            const Dim3 lo = lbound(bx);
            const Dim3 hi = ubound(bx);
            const auto &arr = S_plot.array(mfi);

            for (int k = lo.z; k <= hi.z; k++)
            {
                for (int j = lo.y; j <= hi.y; j++)
                {
                    const Real y = probLoY + (double(j) + 0.5) * dY;
                    for (int i = lo.x; i <= hi.x; i++)
                    {
                        const Real x = probLoX + (double(i) + 0.5) * dX;

                        double rho = arr(i, j, k, 0);
                        double xmomentum = arr(i, j, k, 1);
                        double ymomentum = arr(i, j, k, 2);
                        double ener = arr(i, j, k, 3);
                        double rhoYO2 = arr(i, j, k, 4);
                        double rhoYN2 = arr(i, j, k, 5);

                        double vx = xmomentum / rho;
                        double vy = ymomentum / rho;
                        double YO2 = rhoYO2 / rho;
                        double YN2 = rhoYN2 / rho;
                        double Tgas = Tg(rho, vx, vy, YO2, YN2, ener);
                        double p = pressure(rho, YO2, YN2, Tgas);
                        double eps = specIntEner(rho, vx, vy, ener, p);
                        double a = soundSpeed(p, rho, Tgas, YO2, YN2);
                        double gammaval = a * a * rho / p;

                        AllPrint(approx) << x << " " << rho << " " << vx / sqrt(one_atm_Pa) << " " << p / one_atm_Pa << " "
                                         << eps / one_atm_Pa << " " << YO2 << " " << YN2 << " " << Tgas << " " << gammaval << std::endl;
                    }
                }
            }
        }
        if ((level == 0) && (rank == 0))
        {
            // iter += 1;
            std::cout << "wrote output data to text file" << std::endl;
        }

        approx.close();
    }
}