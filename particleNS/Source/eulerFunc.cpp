#include "eulerFunc.H"
#include "recon.H"
#include "thermoTransport.H"

#include "AMReX_Vector.H"
#include "AMReX_Array.H"
#include "AMReX_REAL.H"
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>


#include <math.h>
#include <iostream>
#include <algorithm> 

extern int enIC;
extern double Gamma;
extern double M_O2;
extern double M_N2;
extern double Y_O2;
extern double Y_N2;
extern double R;
extern double one_atm_Pa;      // one atmosphere in Pa
extern double Mavg;            // average molecular weight of the gas mixture
extern double TgInitial;
extern double TpInitial;

extern int NUM_STATE;
extern int gCells;
extern const int spacedim;

using namespace amrex;

void updateEuler(MultiFab& Sborder, Array<MultiFab, SpaceDim>& fluxes, Vector<double> &qL, Vector<double> &qR, \
                 Vector<double> &fluxvals, const int &d, const double &dt, \
                 const double& dx, const double& dy, const int &euler){
    switch (euler)
    {
        case 0: 
        {
            break;
        }
        case 1: // Use HLLC
        {
            const int iOffset = ( d == 0 ? 1 : 0);
            const int jOffset = ( d == 1 ? 1 : 0);
            const int kOffset = ( d == 2 ? 1 : 0);

            for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Dim3 lo = lbound(bx);
                const Dim3 hi = ubound(bx);

                // Indexable arrays for the data, and the directional flux
                // Based on the vertex-centred definition of the flux array, the
                // data array runs from e.g. [0,N] and the flux array from [0,N+1]
                const auto& arr = Sborder.array(mfi);
                const auto& fluxArr = fluxes[d].array(mfi);
                const auto& fluxArrY = fluxes[d].array(mfi);

                // If slopelimiting has been enabled in the settings file, we run through the following set of loops
                // to compute boundary-reconstructed values -2023W2
                if (d == 0) // x-direction flux calculation using HLLC
                {
                    for(int k = lo.z; k <= hi.z+kOffset; k++)
                    {
                        for(int j = lo.y; j <= hi.y+jOffset; j++)
                        {
                            for(int i = lo.x; i <= hi.x+iOffset; i++)
                            {
                                // std::cout << "\n" << "x location: " << i/128.0 << std::endl;
                                // Conservative flux using HLLC scheme -2023W2
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    qL[h] = arr(i-iOffset,j,k,h);
                                    qR[h] = arr(i,j,k,h);
                                }
                                // Call HLLCflux using the left and right states, qL and qR, and the flux direction
                                // tracked by the for loop variable 'd' to compute the correct fluxes. d=0 corresponds 
                                // to the x-direction, and d=1 corresponds to the y-direction. -2023W2
                                fluxvals = HLLCflux(qL,qR,d);                        
                                
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    fluxArr(i,j,k,h) = fluxvals[h];
                                }
                                // std::cout << "mass flux = " << fluxvals[0] << std::endl;
                                // std::cout << "momentum flux = " << fluxvals[1] << std::endl;
                                // std::cout << "energy flux = " << fluxvals[3] << std::endl;
                                // std::cout << "O2 flux = " << fluxvals[4] << std::endl;
                                // std::cout << "N2 flux = " << fluxvals[5] << std::endl;
                            }
                        }
                    }  
                }
                else // y-direction flux calculation using HLLC
                {
                    for(int k = lo.z; k <= hi.z+kOffset; k++)
                    {
                        for(int i = lo.x; i <= hi.x+iOffset; i++)
                        {
                            for(int j = lo.y; j <= hi.y+jOffset; j++)
                            {
                                // Conservative flux using HLLC scheme -2023W2
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    qL[h] = arr(i,j-jOffset,k,h);
                                    qR[h] = arr(i,j,k,h);
                                }
                                // Call HLLCflux using the left and right states, qL and qR, and the flux direction
                                // tracked by the for loop variable 'd' to compute the correct fluxes. d=0 corresponds 
                                // to the x-direction, and d=1 corresponds to the y-direction. -2023W2
                                fluxvals = HLLCflux(qL,qR,d);
                                
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    fluxArrY(i,j,k,h) = fluxvals[h];
                                }
                                
                            }
                        }
                    }  
                } // close y-direction flux calculation

                for(int k = lo.z; k <= hi.z; k++)
                {
                    for(int j = lo.y; j <= hi.y; j++)
                    {
                        for(int i = lo.x; i <= hi.x; i++)
                        {
                            // Conservative update formula
                            if (d == 0){ // x-direction update
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    arr(i,j,k,h) = arr(i,j,k,h) - (dt / dx) * (fluxArr(i+iOffset, j, k,h) - fluxArr(i,j,k,h));
                                }
                                // std::cout << "density = " << arr(i,j,k,0) << std::endl;
                                // std::cout << "momentum = " << arr(i,j,k,1) << std::endl;
                                // std::cout << "energy = " << arr(i,j,k,3) << std::endl;
                                // std::cout << "O2 conc = " << arr(i,j,k,4) << std::endl;
                                // std::cout << "N2 conc = " << arr(i,j,k,5) << std::endl;
                            }
                            else { // y=direction update
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    arr(i,j,k,h) = arr(i,j,k,h) - (dt / dy) * (fluxArrY(i, j+jOffset, k,h) - fluxArrY(i,j,k,h));
                                }
                            }
                            
                        }
                    }
                } // close loop to update cell-centred values
            }
            break;
        } // close case 1 HLLC
        case 2: // Use MUSCL
        {
            const int iOffset = ( d == 0 ? 1 : 0);
            const int jOffset = ( d == 1 ? 1 : 0);
            const int kOffset = ( d == 2 ? 1 : 0);

            Vector <Vector<double> > u0; 
            u0.resize(3, Vector<double> (NUM_STATE));
            Vector <double> slopeCells(3);
            Vector <double> boundLslice(NUM_STATE);
            Vector <double> boundRslice(NUM_STATE);
            Vector <double> boundLsliceOld(NUM_STATE);
            Vector <double> boundRsliceOld(NUM_STATE);

            double rho,T1,T2,T3;

            std::cout << "in muscl" << std::endl;
            
            for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
            {
                std::cout << "in mfi loop" << std::endl;
                const Box& bx = mfi.tilebox();
                const Dim3 lo = lbound(bx);
                const Dim3 hi = ubound(bx);

                // Indexable arrays for the data, and the directional flux
                // Based on the vertex-centred definition of the flux array, the
                // data array runs from e.g. [0,N] and the flux array from [0,N+1]
                const auto& arr = Sborder.array(mfi);
                const auto& fluxArr = fluxes[d].array(mfi);
                const auto& fluxArrY = fluxes[d].array(mfi);

                if (d==0) // reconstruction and flux calculation procedure in the x-direction
                {
                    for(int k = lo.z; k <= hi.z; k++)
                    {
                        for(int j = lo.y; j <= hi.y; j++)
                        {
                            for(int i = lo.x; i <= hi.x+gCells; i++)
                            {
                            // std::cout << "y index is j = " << j << std::endl;
                            // Take values of energy in i, i+1, i+2 to compute the slope ratio for cell i+1.
                            // Use this slope ratio to reconstruct the left and right boundary
                            // states for cell i+1. Store this information in boundL and boundR
                            // starting with entry i+1. -2023W2

                                // For real gas mixture, apply slope limiting on temperature instead of energy,
                                // since depending on the gas-phase enthalpy, energy will give a negative value,
                                // which may give non-physical humps when slope limiting. -2023W2

                                rho = arr(i-gCells,j,k,0);
                                T1  = Tg(rho,arr(i-gCells,j,k,1)/rho,arr(i-gCells,j,k,2)/rho,\
                                         arr(i-gCells,j,k,4)/rho,arr(i-gCells,j,k,5)/rho,arr(i-gCells,j,k,3));
                                rho = arr(i-gCells+1,j,k,0);
                                T2  = Tg(rho,arr(i-gCells+1,j,k,1)/rho,arr(i-gCells+1,j,k,2)/rho,\
                                         arr(i-gCells+1,j,k,4)/rho,arr(i-gCells+1,j,k,5)/rho,arr(i-gCells+1,j,k,3));
                                rho = arr(i-gCells+2,j,k,0);
                                T3  = Tg(rho,arr(i-gCells+2,j,k,1)/rho,arr(i-gCells+2,j,k,2)/rho,\
                                         arr(i-gCells+2,j,k,4)/rho,arr(i-gCells+2,j,k,5)/rho,arr(i-gCells+2,j,k,3));

                                slopeCells = {T1,T2,T3};

                                for(int h = 0; h < NUM_STATE; h++) // fill matrix u0 with values from arr to calculate slope ratio. -2023W2
                                {
                                    u0[0][h] = arr(i-gCells,j,k,h);
                                    u0[1][h] = arr(i-gCells+1,j,k,h);
                                    u0[2][h] = arr(i-gCells+2,j,k,h);
                                }

                                double w = 0;
                                reconstruct(boundLslice,boundRslice,slopeCells,u0,w,dx,dt,d);

                                // if i>0, i.e., the first value has been computed already, we can proceed with the
                                // flux calculation within the same for loops

                                if ((i>lo.x) && (i<=hi.x+gCells))
                                {
                                    // Conservative flux using HLLC scheme -2023W2
                                    for(int h = 0; h < NUM_STATE; h++)
                                    {
                                        qL[h] = boundRsliceOld[h];
                                        qR[h] = boundLslice[h];
                                    }
                                    // Call HLLCflux using the left and right states, qL and qR, and the flux direction
                                    // tracked by the for loop variable 'd' to compute the correct fluxes. d=0 corresponds 
                                    // to the x-direction, and d=1 corresponds to the y-direction.  -2023W2
                                    fluxvals = HLLCflux(qL,qR,d);
                                    
                                    for(int h = 0; h < NUM_STATE; h++)
                                    {
                                        fluxArr(i-1,j,k,h) = fluxvals[h];
                                        // std::cout << "flux val " << fluxvals[h] << std::endl;
                                    }
                                }

                                // The boundary values of the previous cell are stored to be used in the flux calculations.
                                // The boundary values of the entire domain are not kept while sweeping the domain as we are
                                // only looking to calculate fluxes, which we do store. -2023W2
                                boundLsliceOld = boundLslice;
                                boundRsliceOld = boundRslice;
                            }
                        }
                    }
                } // close d=0 (x-direction flux calculations)

                else // reconstruction and flux calculation procedure in the y-direction
                    // Note that the order of variables to loop through are changed - we start at a given x-coordinate, and
                    // sweep in the y-direction (hence the j loop contained within the i loop), in order to calculate fluxes
                    // in the y-direction. -2023W2
                {
                    for(int k = lo.z; k <= hi.z; k++)
                    {
                        for(int i = lo.x; i <= hi.x; i++)
                        {
                            for(int j = lo.y; j <= hi.y+gCells; j++)
                            {
                                // std::cout << "cell number " << i << " array value " << arr(i-gCells,j,k,0) << std::endl;
                                // Take values of energy in i, i+1, i+2 to compute the slope ratio for cell i+1.
                                // Use this slope ratio to reconstruct the left and right boundary
                                // states for cell i+1. Store this information in boundL and boundR
                                // starting with entry i+1

                                slopeCells = {arr(i,j-gCells,k,3),\
                                              arr(i,j-gCells+1,k,3),\
                                              arr(i,j-gCells+2,k,3)};

                                // std::cout << "Y entry i,j " << i << "," << j << " " << slopeCells[0] << slopeCells[1] << slopeCells[2] << std::endl;

                                for(int h = 0; h < NUM_STATE; h++) // fill matrix u0 with values from arr to calculate slope ratio
                                {
                                    u0[0][h] = arr(i,j-gCells,k,h);
                                    u0[1][h] = arr(i,j-gCells+1,k,h);
                                    u0[2][h] = arr(i,j-gCells+2,k,h);
                                }

                                double w = 0;
                                reconstruct(boundLslice,boundRslice,slopeCells,u0,w,dy,dt,d);


                                // if j>0, i.e., the first value has been computed already, we can proceed with the
                                // flux calculation within the same for loops

                                if ((j>lo.y) && (j<=hi.y+gCells))
                                {
                                    // Conservative flux using HLLC scheme -2023W2
                                    for(int h = 0; h < NUM_STATE; h++)
                                    {
                                        qL[h] = boundRsliceOld[h];
                                        qR[h] = boundLslice[h];
                                    }
                                    // Call HLLCflux using the left and right states, qL and qR, and the flux direction
                                    // tracked by the for loop variable 'd' to compute the correct fluxes. d=0 corresponds 
                                    // to the x-direction, and d=1 corresponds to the y-direction.  -2023W2
                                    fluxvals = HLLCflux(qL,qR,d);
                                    
                                    for(int h = 0; h < NUM_STATE; h++)
                                    {
                                        fluxArrY(i,j-1,k,h) = fluxvals[h];
                                    }
                                }
                                boundLsliceOld = boundLslice;
                                boundRsliceOld = boundRslice;
                            }
                        }
                    }
                } // close y-direction flux calculations

                for(int k = lo.z; k <= hi.z; k++)
                {
                    for(int j = lo.y; j <= hi.y; j++)
                    {
                        for(int i = lo.x; i <= hi.x; i++)
                        {
                            // Conservative update formula
                            if (d == 0){ // x-direction update
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    arr(i,j,k,h) = arr(i,j,k,h) - (dt / dx) * (fluxArr(i+iOffset, j, k,h) - fluxArr(i,j,k,h));
                                    // std::cout << "arr val " << arr(i,j,k,h) << std::endl;
                                }
                            }
                            else { // y=direction update
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    arr(i,j,k,h) = arr(i,j,k,h) - (dt / dy) * (fluxArrY(i, j+jOffset, k,h) - fluxArrY(i,j,k,h));
                                }
                            }
                            
                        }
                    }
                } // close loop to update cell-centred values
            } // close loop through patches
            break;
        } // close MUSCL
        case 3: // Use MUSCL for MPI
        {
            const int iOffset = ( d == 0 ? 1 : 0);
            const int jOffset = ( d == 1 ? 1 : 0);
            const int kOffset = ( d == 2 ? 1 : 0);

            Vector <Vector<double> > u0; 
            u0.resize(3, Vector<double> (NUM_STATE));
            Vector <double> slopeCells(3);
            Vector <double> boundLslice(NUM_STATE);
            Vector <double> boundRslice(NUM_STATE);
            
            for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Dim3 lo = lbound(bx);
                const Dim3 hi = ubound(bx);

                // Indexable arrays for the data, and the directional flux
                // Based on the vertex-centred definition of the flux array, the
                // data array runs from e.g. [0,N] and the flux array from [0,N+1]
                const auto& arr      = Sborder.array(mfi);
                const auto& fluxArr  = fluxes[d].array(mfi);
                const auto& fluxArrY = fluxes[d].array(mfi);
                const auto& boundL   = fluxes[d].array(mfi);
                const auto& boundR   = fluxes[d].array(mfi);

                if (d==0) // reconstruction and flux calculation procedure in the x-direction
                {
                    for(int k = lo.z; k <= hi.z; k++)
                    {
                        for(int j = lo.y; j <= hi.y; j++)
                        {
                            for(int i = lo.x; i <= hi.x+gCells; i++)
                            {
                            // std::cout << "y index is j = " << j << std::endl;
                            // Take values of energy in i, i+1, i+2 to compute the slope ratio for cell i+1.
                            // Use this slope ratio to reconstruct the left and right boundary
                            // states for cell i+1. Store this information in boundL and boundR
                            // starting with entry i+1. -2023W2

                                slopeCells = {arr(i-gCells,j,k,3),arr(i-gCells+1,j,k,3),arr(i-gCells+2,j,k,3)};

                                for(int h = 0; h < NUM_STATE; h++) // fill matrix u0 with values from arr to calculate slope ratio. -2023W2
                                {
                                    u0[0][h] = arr(i-gCells,j,k,h);
                                    u0[1][h] = arr(i-gCells+1,j,k,h);
                                    u0[2][h] = arr(i-gCells+2,j,k,h);
                                }

                                double w = 0;
                                reconstruct(boundLslice,boundRslice,slopeCells,u0,w,dx,dt,d);

                                for(int h = 0; h < NUM_STATE; h++) //
                                {
                                    boundL(i-1,j,k,h) = boundLslice[h];
                                    boundR(i-1,j,k,h) = boundRslice[h];
                                }                                
                            }
                        }
                    }
                    for(int k = lo.z; k <= hi.z; k++)
                    {
                        for(int j = lo.y; j <= hi.y; j++)
                        {
                            for(int i = lo.x; i <= hi.x+iOffset; i++)
                            {
                                for(int h = 0; h < NUM_STATE; h++) // fill matrix u0 with values from arr to calculate slope ratio. -2023W2
                                {
                                    qL[h] = boundR(i-1,j,k,h);
                                    qR[h] = boundL(i,  j,k,h);
                                }                
                                fluxvals = HLLCflux(qL,qR,d);    
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    fluxArr(i,j,k,h) = fluxvals[h];
                                }            
                            }
                        }
                    }
                } // close d=0 (x-direction flux calculations)

                else // reconstruction and flux calculation procedure in the y-direction
                    // Note that the order of variables to loop through are changed - we start at a given x-coordinate, and
                    // sweep in the y-direction (hence the j loop contained within the i loop), in order to calculate fluxes
                    // in the y-direction. -2023W2
                {
                    for(int k = lo.z; k <= hi.z; k++)
                    {
                        for(int i = lo.x; i <= hi.x; i++)
                        {
                            for(int j = lo.y; j <= hi.y+gCells; j++)
                            {
                                // std::cout << "cell number " << i << " array value " << arr(i-gCells,j,k,0) << std::endl;
                                // Take values of energy in i, i+1, i+2 to compute the slope ratio for cell i+1.
                                // Use this slope ratio to reconstruct the left and right boundary
                                // states for cell i+1. Store this information in boundL and boundR
                                // starting with entry i+1

                                slopeCells = {arr(i,j-gCells,k,3),arr(i,j-gCells+1,k,3),arr(i,j-gCells+2,k,3)};
                                // std::cout << "Y entry i,j " << i << "," << j << " " << slopeCells[0] << slopeCells[1] << slopeCells[2] << std::endl;

                                for(int h = 0; h < NUM_STATE; h++) // fill matrix u0 with values from arr to calculate slope ratio
                                {
                                    u0[0][h] = arr(i,j-gCells,k,h);
                                    u0[1][h] = arr(i,j-gCells+1,k,h);
                                    u0[2][h] = arr(i,j-gCells+2,k,h);
                                }

                                double w = 0;
                                reconstruct(boundLslice,boundRslice,slopeCells,u0,w,dy,dt,d);

                                for(int h = 0; h < NUM_STATE; h++) //
                                {
                                    boundL(i-1,j,k,h) = boundLslice[h];
                                    boundR(i-1,j,k,h) = boundRslice[h];
                                }  


                                // if j>0, i.e., the first value has been computed already, we can proceed with the
                                // flux calculation within the same for loops

                                if ((j>lo.y) && (j<=hi.y+gCells))
                                {
                                    // Conservative flux using HLLC scheme -2023W2
                                    for(int h = 0; h < NUM_STATE; h++)
                                    {
                                        qL[h] = boundRslice[h];
                                        qR[h] = boundLslice[h];
                                    }
                                    // Call HLLCflux using the left and right states, qL and qR, and the flux direction
                                    // tracked by the for loop variable 'd' to compute the correct fluxes. d=0 corresponds 
                                    // to the x-direction, and d=1 corresponds to the y-direction.  -2023W2
                                    fluxvals = HLLCflux(qL,qR,d);
                                    
                                    for(int h = 0; h < NUM_STATE; h++)
                                    {
                                        fluxArrY(i,j-1,k,h) = fluxvals[h];
                                    }
                                }
                                // boundLsliceOld = boundLslice;
                                // boundRsliceOld = boundRslice;
                            }
                        }
                    }
                } // close y-direction flux calculations

                for(int k = lo.z; k <= hi.z; k++)
                {
                    for(int j = lo.y; j <= hi.y; j++)
                    {
                        for(int i = lo.x; i <= hi.x; i++)
                        {
                            // Conservative update formula
                            if (d == 0){ // x-direction update
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    arr(i,j,k,h) = arr(i,j,k,h) - (dt / dx) * (fluxArr(i+iOffset, j, k,h) - fluxArr(i,j,k,h));
                                    // std::cout << "arr val " << arr(i,j,k,h) << std::endl;
                                }
                            }
                            else { // y=direction update
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    arr(i,j,k,h) = arr(i,j,k,h) - (dt / dy) * (fluxArrY(i, j+jOffset, k,h) - fluxArrY(i,j,k,h));
                                }
                            }
                            
                        }
                    }
                } // close loop to update cell-centred values
            } // close loop through patches
        } // close MUSCL adapted for MPI
        break;
    } // close euler switch case    

}

Vector<double> setIC(const int dim, const double& probLoX, const double& probHiX, const double& probLoY, const double& probHiY){

    double rhoL, vxL, vyL, pL, YO2L, YN2L, rhoR, vxR, vyR, pR, YO2R, YN2R, TgL, TgR, EL, ER;
    double x0,xEnd,xDisc,tEnd;
    Vector<double> RPLeftRight(16);

    if (dim==1) {
        vyL = 0;
        vyR = 0;
        x0 = probLoX;
        xEnd = probHiX;
        switch (enIC) {
            case 1: //toro test 1
            {
                xDisc = 0.5; tEnd=0.25;
                rhoL = 1.0;   vxL = 0.0;  pL = 1.0;
                rhoR = 0.125; vxR = 0.0;  pR = 0.1;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 2: //toro test 2
            {
                xDisc = 0.5; tEnd=0.25;
                rhoL = 1.0; vxL = -2.0; pL = 0.4;      
                rhoR = 1.0; vxR = 2.0;  pR = 0.4;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 3: //toro test 3
            {
                xDisc = 0.5; tEnd=0.012;
                rhoL = 1.0; vxL = 0.0; pL = 1000.0;    
                rhoR = 1.0; vxR = 0.0; pR = 0.01;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 4: //toro test 4
            {
                xDisc = 0.5; tEnd=0.035;
                rhoL = 1.0; vxL = 0.0; pL = 0.01;      
                rhoR = 1.0; vxR = 0.0; pR = 100.0;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 5: //toro test 5
            {
                xDisc = 0.5; tEnd=0.035;
                rhoL = 5.99924; vxL = 19.5975;  pL = 460.894;      
                rhoR = 5.99242; vxR = -6.19633; pR = 46.0950;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 6: //cyl exp
            {
                xDisc = 0.4; tEnd=0.25;
                rhoL = 1; vxL = 0;  pL = 1;      
                rhoR = 0.125; vxR = 0; pR = 0.1;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 7: //offset for 2-D, revert to toro1
            {
                xDisc = 0.4; tEnd=0.25;
                rhoL = 1; vxL = 0;  pL = 1;      
                rhoR = 0.125; vxR = 0; pR = 0.1;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 8: // Diffusion validation test
            {
                xDisc = xEnd*0.5; tEnd=0.25;
                rhoL = 8.7345e-4; vxL = 0.0; pL = 1.0; 
                rhoR = 8.7345e-4; vxR = 0.0; pR = 1.0;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 9: // Sods shock tube problem with O2-N2 diatomic gas mixture and detailed thermodynamics
            {
                xDisc = xEnd*0.5; tEnd=0.25;
                rhoL = 1.0;   vxL = 0.0;  pL = 1.0*one_atm_Pa; //*one_atm_Pa;
                rhoR = 0.125; vxR = 0.0;  pR = 0.1*one_atm_Pa; //*one_atm_Pa;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                
                // std::cout << "initial condition temperature test" << std::endl;
                TgL = pL*Mavg/(R*rhoL);
                TgR = pR*Mavg/(R*rhoR);; //Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 10: // Particle drag test
            {
                xDisc = 0.5; tEnd=0.25;
                vxL = 5.0;  pL = 1.0*one_atm_Pa; 
                vxR = 5.0;  pR = 1.0*one_atm_Pa; 
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                TgL = TgInitial;
                TgR = TgInitial;
                rhoL = pL*Mavg/(R*TgL);
                rhoR = pR*Mavg/(R*TgR);
                break;
            }
            case 11: // Particle heat test
            {
                xDisc = 0.5; tEnd=0.25;
                vxL = 5.0;  pL = 1.0*one_atm_Pa; 
                vxR = 5.0;  pR = 1.0*one_atm_Pa; 
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                TgL = TgInitial;
                TgR = TgInitial;
                rhoL = pL*Mavg/(R*TgL);
                rhoR = pR*Mavg/(R*TgR);
                break;
            }
            case 12: // Particle ignition test
            {
                xDisc = 0.5; tEnd=0.25;
                vxL = 0.0;  pL = 1.0*one_atm_Pa; 
                vxR = 0.0;  pR = 1.0*one_atm_Pa; 
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                TgL = TgInitial;
                TgR = TgInitial;
                rhoL = pL*Mavg/(R*TgL);
                rhoR = pR*Mavg/(R*TgR);
                break;
            }
            case 13: // Particle combustion test
            {
                xDisc = 0.5; tEnd=0.25;
                vxL = 0.0;  pL = 1.0*one_atm_Pa; 
                vxR = 0.0;  pR = 1.0*one_atm_Pa; 
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                TgL = TgInitial;
                TgR = TgInitial;
                rhoL = pL*Mavg/(R*TgL);
                rhoR = pR*Mavg/(R*TgR);
                break;
            }
            case 14: // Particle combustion test
            {
                xDisc = xEnd*0.05; tEnd=0.25;
                vxL = 0.0;  pL = 1.0*one_atm_Pa; 
                vxR = 0.0;  pR = 1.0*one_atm_Pa; 
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                TgL = 1270;
                TgR = 300;
                rhoL = pL*Mavg/(R*TgL);
                rhoR = pR*Mavg/(R*TgR);
                break;
            }
        }
    } //closes dimension=0 if statement
    else {
        vyL = 0;
        vyR = 0;
        switch (enIC) {
            case 1: //toro test 1
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.25;
                rhoL = 1.0;   vxL = 0.0;  pL = 1.0;
                rhoR = 0.125; vxR = 0.0;  pR = 0.1;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 2: //toro test 2
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.25;
                rhoL = 1.0; vxL = -2.0; pL = 0.4;      
                rhoR = 1.0; vxR = 2.0;  pR = 0.4;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 3: //toro test 3
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.012;
                rhoL = 1.0; vxL = 0.0; pL = 1000;    
                rhoR = 1.0; vxR = 0.0; pR = 0.01;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 4: //toro test 4
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.035;
                rhoL = 1.0; vxL = 0.0; pL = 0.01;      
                rhoR = 1.0; vxR = 0.0; pR = 100.0;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 5: //toro test 5
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.035;
                rhoL = 5.99924; vxL = 19.5975;  pL = 460.894;      
                rhoR = 5.99242; vxR = -6.19633; pR = 46.0950;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 6: //cyl exp
            {
                x0  = -1.0; xEnd = 1.0; xDisc = 0.4; tEnd=0.25;
                rhoL = 1; vxL = 0;  pL = 1;      
                rhoR = 0.125; vxR = 0; pR = 0.1;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 7: //offset for 2-D, revert to toro1
            {
                x0  = 0; xEnd = sqrt(2.0); xDisc = sqrt(2.0)/2; tEnd=0.25;
                rhoL = 1; vxL = 0;  pL = 1;      
                rhoR = 0.125; vxR = 0; pR = 0.1;
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                EL = pL/(Gamma-1);
                ER = pR/(Gamma-1);
                TgL = Tg(rhoL,0.0,0.0,Y_O2,Y_N2,EL);
                TgR = Tg(rhoR,0.0,0.0,Y_O2,Y_N2,ER);
                break;
            }
            case 14:
            {
                x0  = 0; xEnd = 0.00256; xDisc =  0.00064; tEnd=1;
                vxL = 0.0;  pL = 1.0*one_atm_Pa; 
                vxR = 0.0;  pR = 1.0*one_atm_Pa; 
                YO2L = Y_O2; YN2L = Y_N2; YO2R = Y_O2; YN2R = Y_N2;
                TgL = 600;
                TgR = 300;
                rhoL = pL*Mavg/(R*TgL);
                rhoR = pR*Mavg/(R*TgR);
                break;
            }
        }
    }

    // double TgL = (pL*Mavg)/(rhoL*R);
    // double TgR = (pR*Mavg)/(rhoR*R);
    
    RPLeftRight[0] = rhoL;
    RPLeftRight[1] = rhoL*vxL;
    RPLeftRight[2] = rhoL*vyL;
    RPLeftRight[3] = energy(rhoL,vxL,vyL,TgL,YO2L,YN2L,pL);
    RPLeftRight[4] = rhoL*YO2L;
    RPLeftRight[5] = rhoL*YN2L;
    RPLeftRight[6] = rhoR;
    RPLeftRight[7] = rhoR*vxR;
    RPLeftRight[8] = rhoR*vyR;
    RPLeftRight[9] = energy(rhoR,vxR,vyR,TgR,YO2R,YN2R,pR);
    RPLeftRight[10] = rhoR*YO2R;
    RPLeftRight[11] = rhoR*YN2R;
    RPLeftRight[12] = xDisc;
    RPLeftRight[13] = x0;
    RPLeftRight[14] = xEnd;
    RPLeftRight[15] = tEnd;

    return RPLeftRight;    
}

void getStopTime(int enIC, Real& stop_time){
    switch (enIC) {
        case 1: //toro test 1
        {
            stop_time=0.25;
            break;
        }
        case 2: //toro test 2
        {
            stop_time=0.25;
            break;
        }
        case 3: //toro test 3
        {
            stop_time=0.012;
            break;
        }
        case 4: //toro test 4
        {
            stop_time=0.035;
            break;
        }
        case 5: //toro test 5
        {
            stop_time=0.035;
            break;
        }
        case 6: //cyl exp
        {
            stop_time=0.25;
            break;
        }
        case 7: //offset for 2-D, revert to toro1
        {
            stop_time=0.25;
            break;
        }
        case 8: //offset for 2-D, revert to toro1
        {
            stop_time=100;
            break;
        }
        case 9: //offset for 2-D, revert to toro1
        {
            stop_time=0.25/sqrt(one_atm_Pa);
            break;
        }
        case 10:
        {
            stop_time=0.1;
            break;
        }
        case 11:
        {
            stop_time=0.1;
            break;
        }
        case 12:
        {
            stop_time=0.04;
            break;
        }
        case 13:
        {
            stop_time=0.04;
            break;
        }
        case 14:
        {
            stop_time=0.04;
            break;
        }
    }
}


double soundSpeed(const double& p, const double& rho, const double& Tg,\
                  const double& YO2, const double& YN2){
    double cpval,gammaval,a,Mavgval;
    cpval = cpMix(cpO2(Tg),cpN2(Tg),YO2,YN2);
    Mavgval = getMavg(YO2,YN2);
    gammaval = cpval/(cpval-R/Mavgval);
    a = sqrt(fabs(gammaval*p/rho));

    // std::cout << "Tg is: " << Tg << ", gamma is: " << gammaval << std::endl;
    return a;
}

void wavespeedEstimate(const Vector<double>& qL, const Vector<double>& qR,\
                       double& sL, double& sR, double& sStar, const int dir){

    // This function computes wavespeed estimates for the HLLC solver, retrieved via pass by reference variables sL, sR, and sStar.

    // std::cout << "in wavespeedEstimate calculator" << std::endl;

    double pL,pR,rhoL,rhoR,vxL,vxR,vyL,vyR,vL,vR,aL,aR,TL,TR,YO2R,YN2R,YO2L,YN2L,enerL,enerR;
    double rhobar, abar, pPVRS, pStar, qqL, qqR;    
    Vector<double> primL(NUM_STATE);
    Vector<double> primR(NUM_STATE);

    primL = getPrim(qL);
    primR = getPrim(qR);

    rhoL = qL[0];
    rhoR = qR[0];
    vxL  = primL[1];
    vxR  = primR[1];
    vyL  = primL[2];
    vyR  = primR[2];
    pL   = primL[3];
    pR   = primR[3];
    enerL= qL[3];
    enerR= qR[3];
    YO2L = primL[4];
    YO2R = primR[4];
    YN2L = primL[5];
    YN2R = primR[5];

    if (dir==0){
        vL = vxL; vR = vxR;
    }
    else {
        vL = vyL; vR = vyR;
    }

    TL = Tg(rhoL,vxL,vyL,YO2L,YN2L,enerL);
    TR = Tg(rhoR,vxR,vyR,YO2R,YN2R,enerR);

    aL = soundSpeed(pL,rhoL,TL,YO2L,YN2L);
    aR = soundSpeed(pR,rhoR,TR,YO2R,YN2R);

    // // Here, we use a pressure-based wave speed estimate proposed by Toro et al. (1994, Shock Waves)
    // // The star-region pressure estimate follows the form of the PVRS approximate Riemann solver by Toro (1991, Proc. Roy. Soc. London)

    // rhobar = 0.5*(rhoL+rhoR);
    // abar   = 0.5*(aL+aR);
    // pPVRS  = 0.5*(pL+pR)-0.5*(vR-vL)*rhobar*abar;  
    // pStar  = std::max(0.0,pPVRS);

    // if (pStar <= pL){
    //     qqL = 1;
    // }
    // else {
    //     qqL = sqrt(std::abs(1+((Gamma+1)/(2*Gamma))*((pStar/pL) - 1)));
    // }
    // if (pStar <= pR){
    //     qqR = 1;
    // }
    // else {
    //     qqR = sqrt(std::abs(1+((Gamma+1)/(2*Gamma))*((pStar/pR) - 1)));
    // }

    // sL = vL - aL*qqL;
    // sR = vR + aR*qqR;
    // sStar = (pR-pL+rhoL*vL*(sL-vL)-rhoR*vR*(sR-vR))/(rhoL*(sL-vL)-rhoR*(sR-vR));

    // Here, we use the wavespeed estimate proposed by Davis in "Simplified second-order godunov-type methods,
    // SIAM Journal on Scientific and Statistical Computing  (1988)

    double splus = std::max(std::max(fabs(vL-aL),fabs(vR-aR)),std::max(fabs(vL+aL),fabs(vR+aR)));
    sL    = -splus;
    sR    = splus;
    sStar = (pR-pL+rhoL*vL*(sL-vL)-rhoR*vR*(sR-vR))/(rhoL*(sL-vL)-rhoR*(sR-vR));

    if (fabs(sStar) < 1e-12){
        sStar = 0.0;
    }

    // std::cout << "sL, sR, sStar calculated:" << std::endl;
    // std::cout << sL << " " << sR << " " << sStar << std::endl;
}

Vector<double> HLLCflux(const Vector<double>& qL, const Vector<double>& qR, const int& dir){

    // This function calculates the HLLC flux based on left and right states, qL and qR, respectively

    // std::cout << "in HLLCflux calculator" << std::endl;
    
    double sL,sR,sStar;
    Vector<double> fL(NUM_STATE);
    Vector<double> fR(NUM_STATE);
    Vector<double> fLstar(NUM_STATE);
    Vector<double> fRstar(NUM_STATE);
    Vector<double> fluxvals(NUM_STATE);

    wavespeedEstimate(qL,qR,sL,sR,sStar,dir);

    fL = getEulerFlux(qL,dir);
    fR = getEulerFlux(qR,dir);
    fLstar = HLLCstarFlux(fL,qL,sL,sStar,dir);
    fRstar = HLLCstarFlux(fR,qR,sR,sStar,dir);

    if(sL >= 0)
    {	
        fluxvals = fL;
        // std::cout << "energy flux: " << fluxvals[3] << std::endl;
        // std::cout << "left state" << std::endl;
    }

    else if((sL < 0)&&(sStar >= 0))
    {
        fluxvals = fLstar;
        // std::cout << "energy flux: " << fluxvals[3] << std::endl;
        // std::cout << "left star state" << std::endl;
    }

    else if((sStar < 0)&&(sR >= 0))
    {
        fluxvals = fRstar;
        // std::cout << "energy flux: " << fluxvals[3] << std::endl;
        // std::cout << "right star state" << std::endl;
    }
    
    else //(sR < 0)
    {	
        fluxvals = fR;
        // std::cout << "energy flux: " << fluxvals[3] << std::endl;
        // std::cout << "right state" << std::endl;
    }

    // #ifdef DEBUG
        // std::cout << "fluxes calculated" << std::endl;
    // #endif

    return fluxvals;

}

Vector<double> HLLCstarFlux(const Vector<double>& f, const Vector<double>& q,\
                                 const double& s, const double& sStar, const int& dir){

    // Calculates HLLC star-state flux for HLLC approximate Riemann solver.

	double rho,vx,vy,p,ener,YO2,YN2,temp;
    Vector<double> starFlux(NUM_STATE);
	Vector<double> prim(NUM_STATE);
	
	prim = getPrim(q);
	
	rho  = prim[0];
	vx   = prim[1];
    vy   = prim[2];
	p    = prim[3];
    YO2  = prim[4];
    YN2  = prim[5];
    ener = q[3];
	
    if (dir == 0){  // x-direction flux
        temp = rho*(s-vx)/(s-sStar);
        starFlux[0] = f[0]+s*(temp-q[0]);
        starFlux[1] = f[1]+s*(sStar*temp-q[1]);
        starFlux[2] = f[2]+s*(vy*temp-q[2]);
        starFlux[3] = f[3]+s*(temp*((ener/rho)+(sStar-vx)*(sStar+p/(rho*(s-vx))))-q[3]);
        starFlux[4] = f[4]+s*(temp*YO2-q[4]);
        starFlux[5] = f[5]+s*(temp*YN2-q[5]);
    }
    else {  // y-direction flux
        temp = rho*(s-vy)/(s-sStar);
        starFlux[0] = f[0]+s*(temp-q[0]);
        starFlux[1] = f[1]+s*(vx*temp-q[1]);
        starFlux[2] = f[2]+s*(sStar*temp-q[2]);
        starFlux[3] = f[3]+s*(temp*((ener/rho)+(sStar-vy)*(sStar+p/(rho*(s-vy))))-q[3]);
        starFlux[4] = f[4]+s*(temp*YO2-q[4]);
        starFlux[5] = f[5]+s*(temp*YN2-q[5]);
    }
	
	
	return starFlux;
	
}

Vector<double> HLLflux(const Vector<double>& qL, const Vector<double>& qR, const int& dir){

    // This function calculates the HLLC flux based on left and right states, qL and qR, respectively
    // Although this function is not used, it was initially implemented to test the functionality of
    // the approximate Riemann solvers.

    double sL,sR,sStar;
    Vector<double> fL(NUM_STATE);
    Vector<double> fR(NUM_STATE);
    Vector<double> fstar(NUM_STATE);
    Vector<double> fluxvals(NUM_STATE);

    wavespeedEstimate(qL,qR,sL,sR,sStar,dir);

    fL = getEulerFlux(qL,dir);
    fR = getEulerFlux(qR,dir);
    fstar = HLLstarFlux(fL,fR,qL,qR,sL,sR,dir);
    

    if(sL > 0)
    {	
        fluxvals = fL;
        // std::cout << "Eflux: " << fluxvals[3] << std::endl;
    }

    else if((sL <= 0)&&(sR > 0))
    {
        fluxvals = fstar;
        // std::cout << "Eflux: " << fluxvals[3] << std::endl;
    }
    
    else //(sR < 0)
    {	
        fluxvals = fR;
        // std::cout << "Eflux: " << fluxvals[3] << std::endl;
    }

    return fluxvals;

}

Vector<double> HLLstarFlux(const Vector<double>& fL, const Vector<double>& fR, const Vector<double>& qL,\
                           const Vector<double>& qR, const double& sL, const double& sR, const int& dir){

    // Calculates star-state flux for HLL approximate Riemann solver

    Vector<double> starFlux(NUM_STATE);
	
    if (dir == 0){  // x-direction flux
        for (int h=0; h<NUM_STATE; h++)
        {
            starFlux[h] = (sR*fL[h]-sL*fR[h]+sL*sR*(qR[h]-qL[h]))/(sR-sL);
        }
    }
    else {  // y-direction flux
        for (int h=0; h<NUM_STATE; h++)
        {
            starFlux[h] = (sR*fL[h]-sL*fR[h]+sL*sR*(qR[h]-qL[h]))/(sR-sL);
        }
    }
	
	return starFlux;
	
}

Vector<double> getEulerFlux(const Vector<double>& q, const int& dir){

    // Takes in conservative variable vector and directional input, returns conservative flux vector
    // based on conservation form of Euler equations.
    
    Vector<double> EulerFlux(NUM_STATE);
    Vector<double> prim(NUM_STATE);
    prim = getPrim(q);

    double rho  = prim[0];
    double vx   = prim[1];
    double vy   = prim[2];
    double p    = prim[3];
    double YO2  = prim[4];
    double YN2  = prim[5];
    
    double ener = q[3];

    
    if (dir == 0){
        EulerFlux[0] = rho*vx;
        EulerFlux[1] = rho*vx*vx+p;
        EulerFlux[2] = rho*vx*vy;
        EulerFlux[3] = (ener+p)*vx;
        EulerFlux[4] = rho*YO2*vx;
        EulerFlux[5] = rho*YN2*vx;
    }
    else {
        EulerFlux[0] = rho*vy;
        EulerFlux[1] = rho*vx*vy;
        EulerFlux[2] = rho*vy*vy+p;
        EulerFlux[3] = (ener+p)*vy;
        EulerFlux[4] = rho*YO2*vy;
        EulerFlux[5] = rho*YN2*vy;
    }
    
    return EulerFlux;

}

Vector<double> getPrim(const Vector<double>& q){

    // Convert conserved -> primitive variables
    
    Vector<double> prim(NUM_STATE);

    double rho   = q[0];
    double rhovx = q[1];
    double rhovy = q[2];
    double ener  = q[3];
    double rhoYO2 = q[4];
    double rhoYN2 = q[5];

    double vx = rhovx/rho;
    double vy = rhovy/rho;
    double YO2 = rhoYO2/rho;
    double YN2 = rhoYN2/rho;

    double Tgas = Tg(rho,vx,vy,YO2,YN2,ener);
    double p = pressure(rho,YO2,YN2,Tgas);

    // std::cout << "Tgas and pressure are: " << Tgas << " " << p << std::endl;

    prim[0] = rho;
    prim[1] = vx;
    prim[2] = vy;
    prim[3] = p;
    prim[4] = YO2;
    prim[5] = YN2;

    return prim;
}

Vector<double> getCons(const Vector<double>& prim){

    // Convert primitive -> conserved variables

    Vector<double> q(NUM_STATE);

    double rho = prim[0];
    double vx  = prim[1];
    double vy  = prim[2];
    double p   = prim[3];
    double YO2 = prim[4];
    double YN2 = prim[5];

    double Mavg = getMavg(YO2,YN2);
    double Tg   = p*Mavg/(rho*R);

    q[0] = rho;
    q[1] = rho*vx;
    q[2] = rho*vy;
    q[3] = energy(rho,vx,vy,Tg,YO2,YN2,p);
    q[4] = rho*YO2;
    q[5] = rho*YN2;

    return q;
}

// Below are functions that convert between various primitive and conserved variables,
// primarily used when outputting values to text files for plotting.

double energy(const double& rho, const double& u, const double& v, const double& Tg, \
              const double& YO2, const double& YN2, const double& p){
    double vT  = sqrt(u*u + v*v);
    double ene = Hgas(rho,YO2,YN2,Tg) + 0.5*rho*vT*vT - p;    // total energy = internal energy + kinetic energy
    return ene;
}

// double energy(const double& rho, const double& u, const double& v, const double& Tg, \
//               const double& YO2, const double& YN2, const double& p){
//     double vT  = sqrt(u*u + v*v);
//     double ene = Hgas(rho,YO2,YN2,Tg) + 0.5*rho*vT*vT;    // total energy = internal energy + kinetic energy
//     return ene;
// }

double pressure(const double& rho, const double& YO2, const double& YN2, const double& Tg){
    double Mavg = getMavg(YO2,YN2);
    double pres = rho*(R/Mavg)*Tg;
    return pres;
}

// double pressure(const double& rho, const double& u, const double& v, const double& Tg, \
//                 const double& YO2, const double& YN2, const double& ener){
//     double vT  = sqrt(u*u + v*v);
//     double pres = Hgas(rho,YO2,YN2,Tg) + 0.5*rho*vT*vT - ener;
//     return pres;
// }

double specIntEner(const double& rho, const double& vx, const double& vy, const double& energy, const double& p){
    double v       = pow(vx*vx + vy*vy,0.5);
    double specInt = (energy-0.5*rho*v*v)/(rho);
    return specInt;
}