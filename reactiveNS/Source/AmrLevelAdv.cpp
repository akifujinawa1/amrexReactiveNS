
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

#include <eulerFunc.H>
#include <exactFunc.H>
#include <recon.H>
#include <diffusionFunc.H>
#include <source.H>
#include <constants.H>
#include <math.h>
#include <iostream>
#include <algorithm> 

// using namespace amrex;

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
extern   double       Gamma;
extern   double       R;
extern   double       one_atm_Pa;
extern   double       M;
extern   double       T0;

// define the remaining global variables here. NUM_GROW should be defined based on the value of slope limiting.
int      AmrLevelAdv::verbose         = 0;
Real     AmrLevelAdv::cfl             = 0.9; // Default value - can be overwritten in settings file
int      AmrLevelAdv::do_reflux       = 1;  
int      AmrLevelAdv::NUM_STATE       = 5;  // set this to 5 for one-step reactive NS. -2023W2
int      AmrLevelAdv::NUM_GROW        = 2;  // number of ghost cells, set gCells from main.cpp file. -2023W2
int      n_cell;
int      max_level;
const int      spacedim               = amrex::SpaceDim;
int      NUM_STATE                    = AmrLevelAdv::NUM_STATE;
int      pfrequency                   = 20;
int      iter                         = 0;
int      printlevel                   = 0;

// A few vector quantities are defined here to be used in the slope reconstruction / flux calculations
Vector<double> qL(NUM_STATE);
Vector<double> qR(NUM_STATE);
Vector<double> qLlo(NUM_STATE);
Vector<double> qRlo(NUM_STATE);
Vector<double> qLhi(NUM_STATE);
Vector<double> qRhi(NUM_STATE);
Vector<double> fluxvals(NUM_STATE);
Vector <double> printCount(pfrequency);
// Vector <Vector<double> > u0; 
// Vector <double> slopeCells(3);
// Vector <double> boundLslice(NUM_STATE);
// Vector <double> boundRslice(NUM_STATE);
// Vector <double> boundLsliceOld(NUM_STATE);
// Vector <double> boundRsliceOld(NUM_STATE);
Vector <double> viscSlice(NUM_STATE);

//
//Default constructor.  Builds invalid object.
//
AmrLevelAdv::AmrLevelAdv ()
{
  // Flux registers store fluxes at patch boundaries to ensure fluxes are conservative between AMR levels
  flux_reg = 0;
}

//
//The basic constructor.
//
AmrLevelAdv::AmrLevelAdv (Amr&  papa,
     	                    int   lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real  time)
  :
  AmrLevel(papa,lev,level_geom,bl,dm,time) 
{
  // Flux registers are only required if AMR is actually used, and if flux fix up is being done (recommended)
  flux_reg = 0;
  if (level > 0 && do_reflux)
  {
    flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
  }
}

//
//The destructor.
//
AmrLevelAdv::~AmrLevelAdv () 
{
    delete flux_reg;
}

//
//Restart from a checkpoint file.
//
// AMReX can save simultion state such
// that if the code crashes, it can be restarted, with different
// settings files parameters if necessary (e.g. to output about the
// point of the crash).
//
void
AmrLevelAdv::restart (Amr&          papa,
	                    std::istream& is,
                      bool          bReadSpecial)
{
  AmrLevel::restart(papa,is,bReadSpecial);
  
  BL_ASSERT(flux_reg == 0);
  if (level > 0 && do_reflux)
  {
    flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
  }
  
}

//
// Write a checkpoint file - format is handled automatically by AMReX
void 
AmrLevelAdv::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old) 
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
}

//
//Write a plotfile to specified directory - format is handled automatically by AMReX.
//

void
AmrLevelAdv::writePlotFile (const std::string& dir,
	 	                        std::ostream&      os,
                            VisMF::How         how)
{
  // AmrLevel::writePlotFile (dir,os,how);
  
  int finest_level = parent->finestLevel();
  
  if (do_reflux && level < finest_level)
    reflux();
  
  if (level < finest_level)
    avgDown();
  
//   // If final time, loop through all patches here.
  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();
  const Real* prob_hi = geom.ProbHi();
  const Real cur_time = state[Phi_Type].curTime();
  const MultiFab& S_plot = get_new_data(Phi_Type);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if ((level == 0)&&(rank == 0)){
    iter += 1;
  }

  std::cout << "Printing data on AMR level: " << level << std::endl;
  
  const Real dX = dx[0];
  const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);

  const Real probLoX = prob_lo[0];
  const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0);


  std::string method;

  std::string resolution = std::to_string(n_cell);
  std::string dimension  = std::to_string(amrex::SpaceDim);
  std::string test       = std::to_string(enIC);
  std::string iteration  = std::to_string(iter);

  std::ofstream approx;

  approx.open("output/txt/test8/time"+iteration+".txt",std::ofstream::app);

  for (MFIter mfi(S_plot); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);
    const auto& arr = S_plot.array(mfi);
    
    for(int k = lo.z; k <= hi.z; k++)
    {
      for(int j = lo.y; j <= hi.y; j++)
      {
        const Real y     = probLoY + (double(j)+0.5) * dY;
        for(int i = lo.x; i <= hi.x; i++)
        {
          const Real x     = probLoX + (double(i)+0.5) * dX;

          double rho       = arr(i,j,k,0);
          double xmomentum = arr(i,j,k,1);
          double ymomentum = arr(i,j,k,2);
          double ener      = arr(i,j,k,3);
          double rholambda = arr(i,j,k,4);

          double vx  = xmomentum/rho;
          double vy  = ymomentum/rho;
          double lambda = rholambda/rho;
          double p   = pressure(rho,vx,vy,ener);  //print in bar
          double T = p/(rho*R/M);
          double eps = specIntEner(rho,vx,vy,ener);

          approx << x << " " << rho << " " << vx << " " << p/101325 << " " << eps << " " << lambda << " " << T << std::endl;

        }
      }
    }
  }
  std::cout << "wrote output data to text file" << std::endl;
  approx.close();

}

//
//Define data descriptors.
//
// This is how the variables in a simulation are defined.  In the case
// of the advection equation, a single variable, phi, is defined.
//
void
AmrLevelAdv::variableSetUp ()
{
  BL_ASSERT(desc_lst.size() == 0);

  // A function which contains all processing of the settings file,
  // setting up initial data, choice of numerical methods and
  // boundary conditions
  read_params();
  
  const int storedGhostZones = 0;
    
  // Setting up a container for a variable, or vector of variables:
  // Phi_Type: Enumerator for this variable type
  // IndexType::TheCellType(): AMReX can support cell-centred and vertex-centred variables (cell centred here)
  // StateDescriptor::Point: Data can be a point in time, or an interval over time (point here)
  // storedGhostZones: Ghost zones can be stored (e.g. for output).  Generally set to zero.
  // NUM_STATE: Number of variables in the variable vector (1 in the case of advection equation,3 for 1-D euler)
  // cell_cons_interp: Controls interpolation between levels - cons_interp is good for finite volume
  desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
			                   StateDescriptor::Point,storedGhostZones,NUM_STATE,
			                   &cell_cons_interp);

  //Set up boundary conditions, all boundaries can be set
  //independently, including for individual variables, but lo (left) and hi (right) are useful ways to
  //store them, for consistent access notation for the boundary
  //locations
  int lo_bc[amrex::SpaceDim];
  int hi_bc[amrex::SpaceDim];
  // AMReX has pre-set BCs, including periodic (int_dir) and transmissive (foextrap)
  for (int i = 0; i < amrex::SpaceDim; ++i) {
    lo_bc[i] = hi_bc[i] = BCType::foextrap;   // we want to use transmissive BCs -2023W2
  }

  // Object for storing all the boundary conditions
  BCRec bc(lo_bc, hi_bc);

  // Set up variable-specific information; needs to be done for each variable in NUM_STATE
  // Phi_Type: Enumerator for the variable type being set
  // 0: Position of the variable in the variable vector.  Single variable for advection.
  // phi: Name of the variable - appears in output to identify what is being plotted
  // bc: Boundary condition object for this variable (defined above)
  // BndryFunc: Function for setting boundary conditions.  For basic BCs, AMReX can handle these automatically

  // We set up variable-specific information for 4 variables, density, x-momentum, y-momentum, and energy. This
  // will be the default - when running 1-D test cases, simply set all y-velocities to zero and do not solve for 
  // the y-flux. -2023W2

  desc_lst.setComponent(Phi_Type, 0, "rho", bc, 
			                  StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 1, "xmomentum", bc, 
			                  StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 2, "ymomentum", bc, 
			                  StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 3, "energy", bc, 
			                  StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, 4, "rhoY", bc, 
			                  StateDescriptor::BndryFunc(nullfill));
}

//
//Cleanup data descriptors at end of run.
//
void
AmrLevelAdv::variableCleanUp () 
{
    desc_lst.clear();
}

//
//Initialize grid data at problem start-up.
//
void
AmrLevelAdv::initData ()
{
  //
  // Loop over grids, call FORTRAN function to init with data.
  //
  const Real* dx  = geom.CellSize();
  // Position of the bottom left corner of the domain
  const Real* prob_lo = geom.ProbLo();
  const Real* prob_hi = geom.ProbHi();
  // Create a multifab which can store the initial data
  MultiFab& S_new = get_new_data(Phi_Type);
  Real cur_time   = state[Phi_Type].curTime();

  // amrex::Print works like std::cout, but in parallel only prints from the root processor
  if (verbose) {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

  // Set Gamma to required value based on initial condition here


  // Slightly messy way to ensure uninitialised data is not used.
  // AMReX has an XDim3 object, but a function needs to be written to
  // convert Real* to XDim3
  const Real dX = dx[0];
  const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);

  const Real probLoX = prob_lo[0];
  const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0);
  
  // Initialize vector to store left and right state for Riemann problem
  Vector<double> RPLeftRight(14);
  RPLeftRight = setIC(amrex::SpaceDim);
  double xDisc = RPLeftRight[10];

  // Loop over all the patches at this level
  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto& arr = S_new.array(mfi);

    for(int k = lo.z; k <= hi.z; k++)
    {
      for(int j = lo.y; j <= hi.y; j++)
      {
        const Real y = probLoY + (double(j)+0.5) * dY;
        for(int i = lo.x; i <= hi.x; i++)
        {
          const Real x = probLoX + (double(i)+0.5) * dX;
          // here we can assign initial conditions to the variables we created. 
          
          // for 1-D tests:
          if(amrex::SpaceDim == 1)
          {
            // Implement new initial conditions here, take IC from FV_func.cpp,
            // which are converted to conserved variables, apply to arr(i,j,k,NUM_STATES) -2023W2
            if (enIC == 8){
              double xDomain = prob_hi-prob_lo;
              double Tstar   = 1500.0;
              double L       = (3.0/4.0)*xDomain;
              double T       = std::max(T0,Tstar - (Tstar - T0)*(x/L));
              
              // density
              arr(i,j,k,0) = one_atm_Pa/((R/M)*T);

              // x-momentum
              arr(i,j,k,1) = 0.0;

              // y-momentum
              arr(i,j,k,2) = 0.0;

              // energy
              arr(i,j,k,3) = energy(arr(i,j,k,0),0.0,0.0,one_atm_Pa);

              // density*unburned fuel mass fraction
              arr(i,j,k,4) = arr(i,j,k,0)*1.0;

              // std::cout << "x: " << x << "Density: " << arr(i,j,k,0) << " Temperature: " << T << " Energy: " << arr(i,j,k,3) << std::endl; 
            }
            else {
              if (x<xDisc){
                arr(i,j,k,0) = RPLeftRight[0];
                arr(i,j,k,1) = RPLeftRight[1];
                arr(i,j,k,2) = RPLeftRight[2];
                arr(i,j,k,3) = RPLeftRight[3];
                arr(i,j,k,4) = RPLeftRight[4];
              }
              else {
                arr(i,j,k,0) = RPLeftRight[5];
                arr(i,j,k,1) = RPLeftRight[6];
                arr(i,j,k,2) = RPLeftRight[7];
                arr(i,j,k,3) = RPLeftRight[8];
                arr(i,j,k,4) = RPLeftRight[9];
              }
            }
          }
          else // for 2-D tests:
          {
            if (enIC<6){  //Toro's tests in 2-D
              if (x<xDisc){
                arr(i,j,k,0) = RPLeftRight[0];
                arr(i,j,k,1) = RPLeftRight[1];
                arr(i,j,k,2) = RPLeftRight[2];
                arr(i,j,k,3) = RPLeftRight[3];
                arr(i,j,k,4) = RPLeftRight[4];
              }
              else {
                arr(i,j,k,0) = RPLeftRight[5];
                arr(i,j,k,1) = RPLeftRight[6];
                arr(i,j,k,2) = RPLeftRight[7];
                arr(i,j,k,3) = RPLeftRight[8];
                arr(i,j,k,4) = RPLeftRight[9];
              }
            }
            else if (enIC == 6){  //cylindrical explosion
              if (x*x+y*y<=xDisc*xDisc){
                arr(i,j,k,0) = RPLeftRight[0];
                arr(i,j,k,1) = RPLeftRight[1];
                arr(i,j,k,2) = RPLeftRight[2];
                arr(i,j,k,3) = RPLeftRight[3];
                arr(i,j,k,4) = RPLeftRight[4];
              }
              else {
                arr(i,j,k,0) = RPLeftRight[5];
                arr(i,j,k,1) = RPLeftRight[6];
                arr(i,j,k,2) = RPLeftRight[7];
                arr(i,j,k,3) = RPLeftRight[8];
                arr(i,j,k,4) = RPLeftRight[9];
              }
              //std::cout << "x,y=" << x << "," << y << ", rho=" << arr(i,j,k,0) << std::endl;
            }
            else {  //case where discontinuity is not aligned with the grid
              if ((y-x)>0){
                arr(i,j,k,0) = RPLeftRight[0];
                arr(i,j,k,1) = RPLeftRight[1];
                arr(i,j,k,2) = RPLeftRight[2];
                arr(i,j,k,3) = RPLeftRight[3];
                arr(i,j,k,4) = RPLeftRight[4];
              }
              else {
                arr(i,j,k,0) = RPLeftRight[5];
                arr(i,j,k,1) = RPLeftRight[6];
                arr(i,j,k,2) = RPLeftRight[7];
                arr(i,j,k,3) = RPLeftRight[8];
                arr(i,j,k,4) = RPLeftRight[9];
              }
            }
            
          }
        }
      }
    }
  }

  // A few vector quantities are resized here, only once in the code
  

  if (verbose) {
    amrex::Print() << "Done initializing the level " << level 
		   << " data " << std::endl;
  }  
  
}

//
//Initialize data on this level from another AmrLevelAdv (during regrid).
// These are standard AMReX commands which are unlikely to need altering
//
void
AmrLevelAdv::init (AmrLevel &old)
{
  
  AmrLevelAdv* oldlev = (AmrLevelAdv*) &old;
  //
  // Create new grid data by fillpatching from old.
  //
  Real dt_new    = parent->dtLevel(level);
  Real cur_time  = oldlev->state[Phi_Type].curTime();
  Real prev_time = oldlev->state[Phi_Type].prevTime();
  Real dt_old    = cur_time - prev_time;
  setTimeLevel(cur_time,dt_old,dt_new);
  
  MultiFab& S_new = get_new_data(Phi_Type);

  const int zeroGhosts = 0;
  // FillPatch takes the data from the first argument (which contains
  // all patches at a refinement level) and fills (copies) the
  // appropriate data onto the patch specified by the second argument:
  // old: Source data
  // S_new: destination data
  // zeroGhosts: If this is non-zero, ghost zones could be filled too - not needed for init routines
  // cur_time: AMReX can attempt interpolation if a different time is specified - not recommended for advection eq.
  // Phi_Type: Specify the type of data being set
  // 0: This is the first data index that is to be copied
  // NUM_STATE: This is the number of states to be copied
  FillPatch(old, S_new, zeroGhosts, cur_time, Phi_Type, 0, NUM_STATE);

  // Note: In this example above, the all states in Phi_Type (which is
  // only 1 to start with) are being copied.  However, the FillPatch
  // command could be used to create a velocity vector from a
  // primitive variable vector.  In this case, the `0' argument is
  // replaced with the position of the first velocity component in the
  // primitive variable vector, and the NUM_STATE arguement with the
  // dimensionality - this argument is the number of variables that
  // are being filled/copied, and NOT the position of the final
  // component in e.g. the primitive variable vector.
}

//
// Initialize data on this level after regridding if old level did not previously exist
// These are standard AMReX commands which are unlikely to need altering
//
void
AmrLevelAdv::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);

    // See first init function for documentation
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Advance grids at this level in time.
//  This function is the one that actually calls the flux functions.
//
Real
AmrLevelAdv::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{

  MultiFab& S_mm = get_new_data(Phi_Type);

  // Note that some useful commands exist - the maximum and minumum
  // values on the current level can be computed directly - here the
  // max and min of variable 0 are being calculated, and output.
  Real maxval = S_mm.max(0);
  Real minval = S_mm.min(0);
  amrex::Print() << "phi max = " << maxval << ", min = " << minval  << std::endl;

  // This ensures that all data computed last time step is moved from
  // `new' data to `old data' - this should not need changing If more
  // than one type of data were declared in variableSetUp(), then the
  // loop ensures that all of it is updated appropriately
  for (int k = 0; k < NUM_STATE_TYPE; k++) {
    state[k].allocOldData();
    state[k].swapTimeLevels(dt);
  }

  // S_new is the MultiFab that will be operated upon to update the data
  MultiFab& S_new = get_new_data(Phi_Type);

  const Real prev_time = state[Phi_Type].prevTime();
  const Real cur_time = state[Phi_Type].curTime();
  const Real ctr_time = 0.5*(prev_time + cur_time);

  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();

  //
  // Get pointers to Flux registers, or set pointer to zero if not there.
  //
  FluxRegister *fine    = 0;
  FluxRegister *current = 0;
    
  int finest_level = parent->finestLevel();

  // If we are not on the finest level, fluxes may need correcting
  // from those from finer levels.  To start this process, we set the
  // flux register values to zero
  if (do_reflux && level < finest_level) {
    fine = &getFluxReg(level+1);
    fine->setVal(0.0);
  }

  // If we are not on the coarsest level, the fluxes are going to be
  // used to correct those on coarser levels.  We get the appropriate
  // flux level to include our fluxes within
  if (do_reflux && level > 0)
  {
    current = &getFluxReg(level);
  }

  // Set up a dimensional multifab that will contain the fluxes
  // MultiFab fluxes[amrex::SpaceDim];
  Array<MultiFab, AMREX_SPACEDIM> fluxes;

  // Define the appropriate size for the flux MultiFab.
  // Fluxes are defined at cell faces - this is taken care of by the
  // surroundingNodes(j) command, ensuring the size of the flux
  // storage is increased by 1 cell in the direction of the flux.
  // This is only needed if refluxing is happening, otherwise fluxes
  // don't need to be stored, just used
  for (int j = 0; j < amrex::SpaceDim; j++)
  {
    BoxArray ba = S_new.boxArray();
    ba.surroundingNodes(j);
    fluxes[j].define(ba, dmap, NUM_STATE, 0);
  }

  // Advection velocity - AMReX allows the defintion of a vector
  // object (similar functionality to C++ std::array<N>, since its size must
  // be known, but was implemented before array was added to C++)

  // The vel vector code has been commented out as it is no longer necessary - we 
  // calculate the conservative fluxes using the HLLC method for the Euler eqns. -2023W2
  //const Vector<Real> vel{1.0,0.0,0.0};

  // State with ghost cells - this is used to compute fluxes and perform the update.
  MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
  // See init function for details about the FillPatch function
  FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  

  // We need to implement transmissive boundary conditions, this is defined as 
  // 'foextrap' in AMReX (strange name...). AMReX_BCUtil.H contains header files required
  // for setting boundary conditions (setBC, setHi, setLo, etc.). Following the convention 
  // in AMReX_BCRec.cpp, we define a BCRec vector of size NUM_STATE for Euler Eqn. compatibility.
  // Note that the line of code with geometry.is_periodic has been commented out in the inputs file. -2023W2

  Vector<BCRec> BCVec(NUM_STATE);

  // We now loop through these variables, as well as the total number of spatial dimensions
  // to update the lo and hi (each side of the boundary) boundary conditions
  // The setBC function contains the following lines of code:
  // bcr[dc].setLo(dir, ( bxlo[dir]<=dlo[dir]
  //                      ? bc_dom[sc].lo(dir) : BCType::int_dir ));
  // bcr[dc].setHi(dir, ( bxhi[dir]>=dhi[dir]
  //                      ? bc_dom[sc].hi(dir) : BCType::int_dir ));
  // This is now adapted to our code using BCType::foextrap -2023W2

  // if (enIC < 8) // For the Euler equation tests, use transmissive BCs everywhere
  // {
    // for (int nVar = 0; nVar < NUM_STATE; nVar++){
    //   for (int nDim = 0; nDim < amrex::SpaceDim; nDim++){
    //     BCVec[nVar].setLo(nDim,BCType::foextrap);
    //     BCVec[nVar].setHi(nDim,BCType::foextrap);
    //   }
    // }
  // }
  // else // For the shock-induced igntion problem, use the reflective BC at the left wall
  // {
    for (int nVar = 0; nVar < NUM_STATE; nVar++){
      for (int nDim = 0; nDim < amrex::SpaceDim; nDim++){
        if ((nVar == 1)){
          BCVec[nVar].setLo(nDim,BCType::reflect_odd);
          BCVec[nVar].setHi(nDim,BCType::foextrap);
        }
        else{
          BCVec[nVar].setLo(nDim,BCType::reflect_even);
          BCVec[nVar].setHi(nDim,BCType::foextrap);
        }
      }
    }
  // }

  // Fill periodic boundaries where they exist.  More accurately, the
  // FillBoundary call will fill overlapping boundaries (with periodic
  // domains effectively being overlapping).  It also takes care of
  // AMR patch and CPU boundaries.
  Sborder.FillBoundary(geom.periodicity());

  // We need to fill the domain boundary using FillDomainBoundary
  // in order to fill cell-centered data outside the physical domain, following
  // instructions on AMReX_BCUtil.H -2023W2
  FillDomainBoundary(Sborder,geom,BCVec);

  // The following section of the code calculates the conservative fluxes, and updates the
  // new cell-centred values. 
  
  // This loops through all spatial dimensions. We use this variable 'd' to track 
  // which direction we are computing fluxes along. -2023W2

  double dX = dx[0], dY;
  if (amrex::SpaceDim == 1){
    dY = dx[0];
  }


  // We update the cell-centred values in a dimensionally-split manner, using an inviscid-viscous-source splitting 

  // ____ EULER ____ //

  for (int d = 0; d < amrex::SpaceDim ; d++)   
  {
    updateEuler(Sborder, fluxes, qL, qR, fluxvals, d, dt, dX, dY, euler);

    Sborder.FillBoundary(geom.periodicity());
    FillDomainBoundary(Sborder,geom,BCVec);   // use this to fill cell-centred data outside domain (AMReX_BCUtil.H) -2023W2
    
    // Flux scaling for reflux. By the size of the boundary through which the flux passes, e.g. the x-flux needs scaling by the dy, dz and dt
    if(do_reflux)
    {
      Real scaleFactor = dt;
      for(int scaledir = 0; scaledir < amrex::SpaceDim; ++scaledir)
      {
        // Fluxes don't need scaling by dx[d]
        if(scaledir == d)
        {
          continue;
        }
        scaleFactor *= dx[scaledir];
      }
      fluxes[d].mult(scaleFactor, 0, NUM_STATE);
    }
  } //this closes the d=0 to d=spacedim loop for the EULER update


  // ____ VISCOUS ____ //

  for (int d = 0; d < amrex::SpaceDim ; d++)   
  {
    updateViscous(Sborder, fluxes, qL, qR, qLlo, qRlo, qLhi, qRhi, viscSlice, \
                  d, dt, dX, dY, amrex::SpaceDim, viscous);

    Sborder.FillBoundary(geom.periodicity());
    FillDomainBoundary(Sborder,geom,BCVec);   // use this to fill cell-centred data outside domain (AMReX_BCUtil.H) -2023W2
    
    // Flux scaling for reflux. By the size of the boundary through which the flux passes, e.g. the x-flux needs scaling by the dy, dz and dt
    if(do_reflux)
    {
      Real scaleFactor = dt;
      for(int scaledir = 0; scaledir < amrex::SpaceDim; ++scaledir)
      {
        // Fluxes don't need scaling by dx[d]
        if(scaledir == d)
        {
          continue;
        }
        scaleFactor *= dx[scaledir];
      }
      fluxes[d].mult(scaleFactor, 0, NUM_STATE);
    }
  } //this closes the d=0 to d=spacedim loop for the EULER update

  // ____ SOURCE ____ //

  if ((source == 1)||(source == 2))
  {
    double dtSource = dt/Da;
    for (int i = 0; i < Da; i++)
    {
      updateSource(Sborder, qL, fluxvals, dtSource, source);
    }
  }

  // The updated data is now copied to the S_new multifab.  This means
  // it is now accessible through the get_new_data command, and AMReX
  // can automatically interpolate or extrapolate between layers etc.
  // S_new: Destination
  // Sborder: Source
  // Third entry: Starting variable in the source array to be copied (the zeroth variable in this case)
  // Fourth entry: Starting variable in the destination array to receive the copy (again zeroth here)
  // NUM_STATE: Total number of variables being copied
  // Sixth entry: Number of ghost cells to be included in the copy (zero in this case, since only real
  //              data is needed for S_new)
  MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, 0);

  // Refluxing at patch boundaries.  Amrex automatically does this
  // where needed, but you need to state a few things to make sure it
  // happens correctly:
  // FineAdd: If we are not on the coarsest level, the fluxes at this level will form part of the correction
  //          to a coarse level
  // CrseInit:  If we are not the finest level, the fluxes at patch boundaries need correcting.  Since we
  //            know that the coarse level happens first, we initialise the boundary fluxes through this
  //            function, and subsequently FineAdd will modify things ready for the correction
  // Both functions have the same arguments:
  // First: Name of the flux MultiFab (this is done dimension-by-dimension
  // Second: Direction, to ensure the correct vertices are being corrected
  // Third: Source component - the first entry of the flux MultiFab that is to be copied (it is possible that
  //        some variables will not need refluxing, or will be computed elsewhere (not in this example though)
  // Fourth: Destinatinon component - the first entry of the flux register that this call to FineAdd sends to
  // Fifth: NUM_STATE - number of states being added to the flux register
  // Sixth: Multiplier - in general, the least accurate (coarsest) flux is subtracted (-1) and the most
  //        accurate (finest) flux is added (+1)
  if (do_reflux) {
    if (current) {
      for (int i = 0; i < amrex::SpaceDim ; i++)
      current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
    }
    if (fine) {
      for (int i = 0; i < amrex::SpaceDim ; i++)
      fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
    }
  }
  return dt;
}

//
//Estimate time step.
// This function is called by all of the other time step functions in AMReX, and is the only one that should
// need modifying
//
Real
AmrLevelAdv::estTimeStep (Real)
{
  // This is just a dummy value to start with 
  Real dt_est  = 1.0e+20;

  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();
  const Real cur_time = state[Phi_Type].curTime();
  const MultiFab& S_new = get_new_data(Phi_Type);
  double sL=0.0,sR=0.0,sStar=0.0,sMax = 0.0,sDiff = 0.0,sMaxDiff=0.0;
  Vector<double> qL(NUM_STATE);
  Vector<double> qR(NUM_STATE);

  // This should not really be hard coded

  // We must loop through all cells in the domain to determine the maximum stable time-step
  // permitted. Loop through cells in the same way we did for fluxes and conservative update. -2023W2


    // Loop over all the patches at this level
    // Note we switch the input to mfi from Sborder to S_new -2023W2
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
      const auto& arr = S_new.array(mfi);
      int dir = 0;
      for(int k = lo.z; k <= hi.z; k++)
      {
        for(int j = lo.y; j <= hi.y; j++)
        {
          for(int i = lo.x; i <= hi.x-1; i++)
          {
            // calculate velocity in this cell
            for(int h = 0; h < NUM_STATE; h++){
              qL[h] = arr(i,j,k,h);
              qR[h] = arr(i+1,j,k,h);
            }
            wavespeedEstimate(qL,qR,sL,sR,sStar,dir);
            sDiff = diffusiveSpeed(qL,qR);
            sMaxDiff = std::max(sMaxDiff,sDiff);
            sMax = std::max(sMax,std::max(std::abs(sL),std::abs(sR)));

            // std::cout << "Wavespeed: " << sMax << ", Diffusive speed: " << sMaxDiff << std::endl;
                   
            i++; // we skip a cell each iteration in the for loop as we use two cells to compute the max wavespeed -2023W2
          }
        }
      }
    }

  //const Real velMag = sqrt(2.);
  for(unsigned int d = 0; d < amrex::SpaceDim; ++d)
  {
    dt_est = std::min(dt_est, cfl*dx[d]/sMax);
    if (viscous == 1){
      dt_est = std::min(dt_est, dx[d]*dx[d]/sMaxDiff);
    }
    
  }
  
  // Ensure that we really do have the minimum across all processors
  ParallelDescriptor::ReduceRealMin(dt_est);
  // dt_est *= cfl;

    
  if (verbose) {
    amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level 
		   << ":  dt_est = " << dt_est << std::endl;
  }

  return dt_est;
}

//
//Compute initial time step.
//
Real
AmrLevelAdv::initialTimeStep ()
{
  return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
AmrLevelAdv::computeInitialDt (int                    finest_level,
	  	                         int                    sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                   stop_time)
{
  //
  // Grids have been constructed, compute dt for all levels.
  //
  // AMReX's AMR Level mode assumes that the time step only needs
  // calculating on the coarsest level - all subsequent time steps are
  // reduced by the refinement factor
  if (level > 0)
    return;

  // Initial guess
  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor   *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_level[i]);
  }

  //
  // Limit dt's by the value of stop_time.
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[Phi_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

//
//Compute new `dt'.
//
void
AmrLevelAdv::computeNewDt (int                   finest_level,
		                       int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
  //
  // We are at the end of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  //
  if (level > 0)
    return;

  // Although we only compute the time step on the finest level, we
  // need to take information from all levels into account.  The
  // sharpest features may be smeared out on coarse levels, so not
  // using finer levels could cause instability
  for (int i = 0; i <= finest_level; i++)
  {
    AmrLevelAdv& adv_level = getLevel(i);
    dt_min[i] = adv_level.estTimeStep(dt_level[i]);
  }

  // A couple of things are implemented to ensure that time step's
  // don't suddenly grow by a lot, as this could lead to errors - for
  // sensible mesh refinement choices, these shouldn't really change
  // anything
  if (post_regrid_flag == 1) 
  {
    //
    // Limit dt's by pre-regrid dt
    //
    for (int i = 0; i <= finest_level; i++)
    {
      dt_min[i] = std::min(dt_min[i],dt_level[i]);
    }
  }
  else 
  {
    //
    // Limit dt's by change_max * old dt
    //
    static Real change_max = 1.1;
    for (int i = 0; i <= finest_level; i++)
    {
      dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
    }
  }
    
  //
  // Find the minimum over all levels
  //
  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_min[i]);
  }

  //
  // Limit dt's by the value of stop_time.
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[Phi_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }
  
  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

//
//Do work after timestep().
// If something has to wait until all processors have done their advance function, the post_timestep function
// is the place to put it.  Refluxing and averaging down are two standard examples for AMR
//
void
AmrLevelAdv::post_timestep (int iteration)
{
  //
  // Integration cycle on fine level grids is complete
  // do post_timestep stuff here.
  //
  // std::cout << "level is " << level << std::endl;
  // To validate the test cases against exact Riemann solutions, we wish to output only the final time-step information
  // in 1-D. Initialize the step-size, lower bound of the domain, current time of simulation, and a new set of patch data,
  // S_plot. We follow the same convention as was used in the time-step and flux calculation loops. -2023W2

  
  int finest_level = parent->finestLevel();
  
  if (do_reflux && level < finest_level)
    reflux();
  
  if (level < finest_level)
    avgDown();
  
//   // If final time, loop through all patches here.
  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();
  const Real* prob_hi = geom.ProbHi();
  const Real cur_time = state[Phi_Type].curTime();
  const MultiFab& S_plot = get_new_data(Phi_Type);

// if (enIC < 8){

//   if (cur_time == stop_time){

//     std::string method;

//     if (euler==1){
//       method = "HLLC";
//     }
//     else if (euler==2) {
//       method = "MUSCL";
//     }

//     std::cout << enIC << std::endl;

//     std::string resolution = std::to_string(n_cell);
//     std::string dimension  = std::to_string(amrex::SpaceDim);
//     std::string test       = std::to_string(enIC);
//     std::string AMR;

//     std::ofstream approx;

//     if (conv == 0) // if we are not testing for convergence, i.e., just outputting spatial distributions to plot
//     {
//       if (max_level == 0)
//       {
//         AMR = "noAMR";
//         approx.open("output/txt/test"+test+"/field"+dimension+method+resolution+AMR+".txt",std::ofstream::out | std::ofstream::trunc);
//       }
//       else
//       {
//         AMR = "AMR";
//         if (level == 1)
//         {
//           approx.open("output/txt/test"+test+"/field"+dimension+method+resolution+AMR+".txt",std::ofstream::out | std::ofstream::trunc);
//         }
//         else
//         {
//           approx.open("output/txt/test"+test+"/field"+dimension+method+resolution+AMR+".txt",std::ofstream::app);
//         }
//       }
//     }
//     else // if we are testing for convergence, i.e., running for various resolutions and calculating error via post processing
//     {
//       if (max_level == 0)
//       {
//         AMR = "noAMR";
//         approx.open("output/txt/conv/field"+dimension+method+resolution+AMR+test+".txt",std::ofstream::out | std::ofstream::trunc);
//       }
//       else
//       {
//         AMR = "AMR";
//         if (level == 1)
//         {
//           approx.open("output/txt/conv/field"+dimension+method+resolution+AMR+test+".txt",std::ofstream::out | std::ofstream::trunc);
//         }
//         else
//         {
//           approx.open("output/txt/conv/field"+dimension+method+resolution+AMR+test+".txt",std::ofstream::app);
//         }
//       }
//     }
    
    
//     const Real dX = dx[0];
//     const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);

//     const Real probLoX = prob_lo[0];
//     const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0);


//     // After AMReX computations are completed, solve for the exact solution using the exact Riemann solver 
//     // contained within the exactFunc.cpp file, where headers are located in exactFunc.H
//     // The exact solver is only run if the test is 1D.

//     if (amrex::SpaceDim == 1)
//     {
//       std::ofstream exact;
//       if (conv == 0)
//       {
//         exact.open("output/txt/test"+test+"/field"+dimension+method+"exact.txt",std::ofstream::out | std::ofstream::trunc);
//       }
//       else 
//       {
//         exact.open("output/txt/conv/field"+dimension+method+test+"exact.txt",std::ofstream::out | std::ofstream::trunc);
//       }
      
//       Vector <Vector<double> > arrExact; 
//       int n_exact = pow(2,12);
      
//       double dXexact = (prob_hi[0]-prob_lo[0])/n_exact;
//       arrExact.resize(n_exact, Vector<double> (NUM_STATE));
//       updateExact(n_exact, dXexact, amrex::SpaceDim, exact, arrExact);
//       exact.close();
//     }
    

//     int uppery;
//     int lowery;

//     std::cout << "start data collection for output here" << std::endl;

//     for (MFIter mfi(S_plot); mfi.isValid(); ++mfi)
//     {
//       Box bx = mfi.tilebox();
//       const Dim3 lo = lbound(bx);
//       const Dim3 hi = ubound(bx);
//       const auto& arr = S_plot.array(mfi);

//       if (amrex::SpaceDim == 1)
//       {
//         lowery = lo.y;
//         uppery = hi.y;
//       }
//       else
//       {
//         if (enIC == 6)  // if cyl expl
//         {
//           lowery = (n_cell/2)-1;
//           uppery = lowery;
//         }
//         else
//         {
//           if (lo.y > 4) // only print if in patches along bottom boundary
//           {
//             continue;
//           }
//           lowery = lo.y;
//           uppery = lo.y; // take 1-D slice at lower boundary
//         }
//       }
      
//       for(int k = lo.z; k <= hi.z; k++)
//       { 
//         for(int j = lowery; j <= uppery; j++)
//           {
//             for(int i = lo.x; i <= hi.x; i++)
//             {
//               const Real x     = probLoX + (double(i)+0.5) * dX;
//               double rho       = arr(i,j,k,0);
//               double xmomentum = arr(i,j,k,1);
//               double ymomentum = arr(i,j,k,2);
//               double energy    = arr(i,j,k,3);
//               double rhoY      = arr(i,j,k,4);

//               double vx  = xmomentum/rho;
//               double vy  = ymomentum/rho;
//               double p   = pressure(rho,vx,vy,energy);
//               double eps = specIntEner(rho,vx,vy,energy);
//               double Y   = rhoY/rho;
//               double T   = p/((R/M)*rho);

//               approx << x << " " << rho << " " << vx << " " << p << " " << eps << " " << Y << " " << T << std::endl;

//             }
//           }
//       }
//     }
 
//     std::cout << "wrote output data to text file" << std::endl;
//     approx.close();

//   }
// }

// // else {
// //     double val = cur_time/stop_time;

// //     if (val > (1.0/pfrequency)*(iter-1)){

// //       if (level == 1){
// //         printIter += 1;
// //       }
      
// //       std::cout << "Printing data on AMR level: " << level << std::endl;
      
// //       const Real dX = dx[0];
// //       const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);

// //       const Real probLoX = prob_lo[0];
// //       const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0);


// //       std::string method;

// //       if (euler==0){
// //         method = "HLLC";
// //       }
// //       else {
// //         method = "MUSCL";
// //       }

// //       std::string resolution = std::to_string(n_cell);
// //       std::string dimension  = std::to_string(amrex::SpaceDim);
// //       std::string test       = std::to_string(enIC);
// //       std::string iteration  = std::to_string(iter);

// //       std::ofstream approx;

// //       // if (level == 0){
// //       //   approx.open("output/txt/test8/time"+iteration+".txt",std::ofstream::out | std::ofstream::trunc);
// //       // }
// //       // else{
// //         approx.open("output/txt/test8/time"+iteration+".txt",std::ofstream::app);
// //       // }

// //       // approx.open("output/txt/test8/time"+iteration+".txt",std::ofstream::app);
// //       //std::ofstream::out | std::ofstream::trunc

// //       for (MFIter mfi(S_plot); mfi.isValid(); ++mfi)
// //       {
// //         Box bx = mfi.tilebox();
// //         const Dim3 lo = lbound(bx);
// //         const Dim3 hi = ubound(bx);
// //         const auto& arr = S_plot.array(mfi);
        
// //         for(int k = lo.z; k <= hi.z; k++)
// //         {
// //           for(int j = lo.y; j <= hi.y; j++)
// //           {
// //             const Real y     = probLoY + (double(j)+0.5) * dY;
// //             for(int i = lo.x; i <= hi.x; i++)
// //             {
// //               const Real x     = probLoX + (double(i)+0.5) * dX;

// //               double rho       = arr(i,j,k,0);
// //               double xmomentum = arr(i,j,k,1);
// //               double ymomentum = arr(i,j,k,2);
// //               double ener      = arr(i,j,k,3);
// //               double rholambda = arr(i,j,k,4);

// //               double vx  = xmomentum/rho;
// //               double vy  = ymomentum/rho;
// //               double lambda = rholambda/rho;
// //               double p   = pressure(rho,vx,vy,ener);  //print in bar
// //               double T = p/(rho*R/M);
// //               double eps = specIntEner(rho,vx,vy,ener);

// //               approx << x << " " << rho << " " << vx << " " << p/101325 << " " << eps << " " << lambda << " " << T << std::endl;

// //             }
// //           }
// //         }
// //       }
// //       std::cout << "wrote output data to text file" << std::endl;
// //       approx.close();

// //       if (level == 0){
// //         iter+=1;
// //       }
      

// //     }
// //   } 

  
}

//
//Do work after regrid().
// Nothing normally needs doing here, but if something was calculated on a per-patch basis, new patches might
// this to be calcuated immediately
//
void
AmrLevelAdv::post_regrid (int lbase, int new_finest)
{

}

//
//Do work after a restart().
// Similar to post_regrid, nothing normally needs doing here
//
void
AmrLevelAdv::post_restart() 
{

}

//
//Do work after init().
// Once new patches have been initialised, work may need to be done to ensure consistency, for example,
// averaging down - though for linear interpolation, this probably won't change anything
//
void
AmrLevelAdv::post_init (Real stop_time)
{
  if (level > 0)
    return;
  //
  // Average data down from finer levels
  // so that conserved data is consistent between levels.
  //
  int finest_level = parent->finestLevel();
  for (int k = finest_level-1; k>= 0; k--)
    getLevel(k).avgDown();
}

//
//Error estimation for regridding.
//  Determine which parts of the domain need refinement
//
void
AmrLevelAdv::errorEst (TagBoxArray& tags,
	               int          clearval,
                       int          tagval,
                       Real         time,
                       int          n_error_buf,
                       int          ngrow)
{
  const Real* dx        = geom.CellSize();
  const Real* prob_lo   = geom.ProbLo();

  MultiFab& S_new = get_new_data(Phi_Type);

  Vector<int> itags;
	
  for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
  {
    const Box&  tilebx  = mfi.tilebox();

    // An AMReX construction, effectively a boolean array which is true in positions that are valid for refinement
    TagBox&     tagfab  = tags[mfi];

    // Traditionally, a lot of the array-based operations in AMReX happened in Fortran.  The standard template
    // for these is short and easy to read, flagging on values or gradients (first order calculation)
    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
    // So we are going to get a temporary integer array.
    tagfab.get_itags(itags, tilebx);
	    
    // data pointer and index space
    int*        tptr    = itags.dataPtr();
    const int*  tlo     = tilebx.loVect();
    const int*  thi     = tilebx.hiVect();

    // Various macros exist to convert the C++ data structures to Fortran
    state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
		BL_TO_FORTRAN_3D(S_new[mfi]),
		&tagval, &clearval, 
		AMREX_ARLIM_3D(tilebx.loVect()), AMREX_ARLIM_3D(tilebx.hiVect()), 
		AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &level);
    //
    // Now update the tags in the TagBox.
    //
    tagfab.tags_and_untags(itags, tilebx);
  }
}

//
// This function reads the settings file
//
void
AmrLevelAdv::read_params ()
{
  // Make sure that this is only done once
  static bool done = false;

  if (done) return;

  done = true;

  // A ParmParse object allows settings, with the correct prefix, to be read in from the settings file
  // The prefix can help identify what a settings parameter is used for
  // AMReX has some default ParmParse names, amr and geometry are two commonly needed ones
  ParmParse pp("adv");   

  // ParmParse has two options; query and get.  Query will only alter
  // a parameter if it can be found (if these aren't in the settings
  // file, then the values at the top of this file will be used).  Get
  // will throw an error if the parameter is not found in the settings
  // file.
  pp.query("v",verbose);
  pp.query("cfl",cfl);
  pp.query("do_reflux",do_reflux);
  // pp.query("slopelimiting",slopelimiting);

  
  // Vector variables can be read in; these require e.g.\ pp.queryarr
  // and pp.getarr, so that the ParmParse object knows to look for
  // more than one variable

  // Geometries can be Cartesian, cylindrical or spherical - some
  // functions (e.g. divergence in linear solvers) are coded with this
  // geometric dependency
  Geometry const* gg = AMReX::top()->getDefaultGeometry();

  // This tutorial code only supports Cartesian coordinates.
  if (! gg->IsCartesian()) {
    amrex::Abort("Please set geom.coord_sys = 0");
  }

  // This tutorial code only supports periodic boundaries.
  // The periodicity is read from the settings file in AMReX source code, but can be accessed here

  // The following lines are commented out since we are not using periodic boundary conditions -2023W2

  // if (! gg->isAllPeriodic()) {
  //   amrex::Abort("Please set geometry.is_periodic = 1 1 1");
  // }

  //
  // read tagging parameters from probin file
  //
  // Tradtionally, the inputs file with ParmParse functionality is handled by C++.  However, a Fortran settings
  // file, by default named probin, can also supply variables.  Mostly used for mesh refinement (tagging) critera
  std::string probin_file("probin");

  ParmParse ppa("amr");
  ppa.query("probin_file",probin_file);
  ppa.query("n_cell",n_cell);
  ppa.query("max_level",max_level);


  int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  // use a fortran routine to
  // read in tagging parameters from probin file
  get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);

}

//
// AMReX has an inbuilt reflux command, but we still have the freedom
// to decide what goes into it (for example, which variables are
// actually refluxed).  This also gives a little flexibility as to
// where flux registers are stored.  In this example, they are stored
// on levels [1,fine] but not level 0.  
//
void
AmrLevelAdv::reflux ()
{
  BL_ASSERT(level<parent->finestLevel());

  const Real strt = amrex::second();

  // Call the reflux command with the appropriate data.  Because there
  // are no flux registers on the coarse level, they start from the
  // first level.  But the coarse level to the (n-1)^th are the ones
  // that need refluxing, hence the `level+1'.  
  getFluxReg(level+1).Reflux(get_new_data(Phi_Type),1.0,0,0,NUM_STATE,geom);
    
  if (verbose)
  {
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    Real      end    = amrex::second() - strt;
    
    ParallelDescriptor::ReduceRealMax(end,IOProc);
    
    amrex::Print() << "AmrLevelAdv::reflux() at level " << level 
		   << " : time = " << end << std::endl;
  }
}

//
// Generic function for averaging down - in this case it just makes sure it doesn't happen on the finest level
//
void
AmrLevelAdv::avgDown ()
{
  if (level == parent->finestLevel())
  {
    return;
  }
  // Can select which variables averaging down will happen on - only one to choose from in this case!
  avgDown(Phi_Type);
}

//
// Setting up the call to the AMReX-implemented average down function
//
void
AmrLevelAdv::avgDown (int state_indx)
{
  // For safety, again make sure this only happens if a finer level exists
  if (level == parent->finestLevel()) return;

  // You can access data at other refinement levels, use this to
  // specify your current data, and the finer data that is to be
  // averaged down
  AmrLevelAdv& fine_lev = getLevel(level+1);
  MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
  MultiFab&  S_crse   = get_new_data(state_indx);

  // Call the AMReX average down function:
  // S_fine: Multifab with the fine data to be averaged down
  // S_crse: Multifab with the coarse data to receive the fine data where necessary
  // fine_lev.geom:  Geometric information (cell size etc.) for the fine level
  // geom: Geometric information for the coarse level (i.e. this level)
  // 0: First variable to be averaged (as not all variables need averaging down
  // S_fine.nComp(): Number of variables to average - this can be computed automatically from a multifab
  // refRatio: The refinement ratio between this level and the finer level
  amrex::average_down(S_fine,S_crse,
		      fine_lev.geom,geom,
		      0,S_fine.nComp(),parent->refRatio(level));
}