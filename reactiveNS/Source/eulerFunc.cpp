#include "eulerFunc.H"
#include "recon.H"
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
        {}
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
                                    fluxvals = HLLflux(qL,qR,d);
                                    
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
        } // close MUSCL
    } // close euler switch case    

}

Vector<double> setIC(const int dim) {

    double rhoL, vxL, vyL, pL, YL, rhoR, vxR, vyR, pR, YR;
    double x0,xEnd,xDisc,tEnd;
    Vector<double> RPLeftRight(14);

    if (dim==1) {
        vyL = 0;
        vyR = 0;
        switch (enIC) {
            case 1: //toro test 1
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.25;
                rhoL = 1.0;   vxL = 0.0;  pL = 1.0;
                rhoR = 0.125; vxR = 0.0;  pR = 0.1;
                break;
            }
            case 2: //toro test 2
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.25;
                rhoL = 1.0; vxL = -2.0; pL = 0.4;      
                rhoR = 1.0; vxR = 2.0;  pR = 0.4;
                break;
            }
            case 3: //toro test 3
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.012;
                rhoL = 1.0; vxL = 0.0; pL = 1000.0;    
                rhoR = 1.0; vxR = 0.0; pR = 0.01;
                break;
            }
            case 4: //toro test 4
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.035;
                rhoL = 1.0; vxL = 0.0; pL = 0.01;      
                rhoR = 1.0; vxR = 0.0; pR = 100.0;
                break;
            }
            case 5: //toro test 5
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.035;
                rhoL = 5.99924; vxL = 19.5975;  pL = 460.894;      
                rhoR = 5.99242; vxR = -6.19633; pR = 46.0950;
                break;
            }
            case 6: //cyl exp
            {
                x0  = 0; xEnd = 1; xDisc = 0.4; tEnd=0.25;
                rhoL = 1; vxL = 0;  pL = 1;      
                rhoR = 0.125; vxR = 0; pR = 0.1;
                break;
            }
            case 7: //offset for 2-D, revert to toro1
            {
                x0  = 0; xEnd = 1; xDisc = 0.4; tEnd=0.25;
                rhoL = 1; vxL = 0;  pL = 1;      
                rhoR = 0.125; vxR = 0; pR = 0.1;
                break;
            }
            case 8: //one step ns Detonation
            {
                x0  = 0; xEnd = 4.0; xDisc = 0.5; tEnd=0.25;
                rhoL = 8.7345e-4; vxL = 0.0; pL = 1.0; YL = 1;
                rhoR = 8.7345e-4; vxR = 0.0; pR = 1.0; YR = 1;
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
                break;
            }
            case 2: //toro test 2
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.25;
                rhoL = 1.0; vxL = -2.0; pL = 0.4;      
                rhoR = 1.0; vxR = 2.0;  pR = 0.4;
                break;
            }
            case 3: //toro test 3
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.012;
                rhoL = 1.0; vxL = 0.0; pL = 1000;    
                rhoR = 1.0; vxR = 0.0; pR = 0.01;
                break;
            }
            case 4: //toro test 4
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.035;
                rhoL = 1.0; vxL = 0.0; pL = 0.01;      
                rhoR = 1.0; vxR = 0.0; pR = 100.0;
                break;
            }
            case 5: //toro test 5
            {
                x0  = 0; xEnd = 1; xDisc = 0.5; tEnd=0.035;
                rhoL = 5.99924; vxL = 19.5975;  pL = 460.894;      
                rhoR = 5.99242; vxR = -6.19633; pR = 46.0950;
                break;
            }
            case 6: //cyl exp
            {
                x0  = -1.0; xEnd = 1.0; xDisc = 0.4; tEnd=0.25;
                rhoL = 1; vxL = 0;  pL = 1;      
                rhoR = 0.125; vxR = 0; pR = 0.1;
                break;
            }
            case 7: //offset for 2-D, revert to toro1
            {
                x0  = 0; xEnd = sqrt(2.0); xDisc = sqrt(2.0)/2; tEnd=0.25;
                rhoL = 1; vxL = 0;  pL = 1;      
                rhoR = 0.125; vxR = 0; pR = 0.1;
                break;
            }
        }
    }

    YL = 1; YR = 1; 

    RPLeftRight[0] = rhoL;
    RPLeftRight[1] = rhoL*vxL;
    RPLeftRight[2] = rhoL*vyL;
    RPLeftRight[3] = energy(rhoL,vxL,vyL,pL);
    RPLeftRight[4] = rhoL*YL;
    RPLeftRight[5] = rhoR;
    RPLeftRight[6] = rhoR*vxR;
    RPLeftRight[7] = rhoR*vyR;
    RPLeftRight[8] = energy(rhoR,vxR,vyR,pR);
    RPLeftRight[9] = rhoR*YR;
    RPLeftRight[10] = xDisc;
    RPLeftRight[11] = x0;
    RPLeftRight[12] = xEnd;
    RPLeftRight[13] = tEnd;

    return RPLeftRight;    
}

void getStopTime(int enIC, Real& stop_time){
    switch (enIC) {
        case 1: //toro test 1
        {
            stop_time=0.05;
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
            stop_time=6e-4;
            break;
        }
    }
}


double soundSpeed(const double& p, const double& rho){
    double a = sqrt(fabs(Gamma*p/rho));
    return a;
}

void wavespeedEstimate(const Vector<double>& qL, const Vector<double>& qR,\
                       double& sL, double& sR, double& sStar, const int dir){

    // This function computes wavespeed estimates for the HLLC solver, retrieved via pass by reference variables sL, sR, and sStar.

    double pL,pR,rhoL,rhoR,vxL,vxR,vyL,vyR,vL,vR,aL,aR;
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


    if (dir==0){
        vL = vxL; vR = vxR;
    }
    else {
        vL = vyL; vR = vyR;
    }
    
    

    aL = soundSpeed(pL,rhoL);
    aR = soundSpeed(pR,rhoR);

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


}

Vector<double> HLLCflux(const Vector<double>& qL, const Vector<double>& qR, const int& dir){

    // This function calculates the HLLC flux based on left and right states, qL and qR, respectively
    
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
        //std::cout << "left state" << std::endl;
    }

    else if((sL < 0)&&(sStar >= 0))
    {
        fluxvals = fLstar;
        //std::cout << "left star state" << std::endl;
    }

    else if((sStar < 0)&&(sR >= 0))
    {
        fluxvals = fRstar;
        //std::cout << "right star state" << std::endl;
    }
    
    else //(sR < 0)
    {	
        fluxvals = fR;
        //std::cout << "right state" << std::endl;
    }

    #ifdef DEBUG
        std::cout << "fluxes calculated" << std::endl;
    #endif

    return fluxvals;

}

Vector<double> HLLCstarFlux(const Vector<double>& f, const Vector<double>& q,\
                                 const double& s, const double& sStar, const int& dir){

    // Calculates HLLC star-state flux for HLLC approximate Riemann solver.

	double rho,vx,vy,p,ener,Y,temp;
    Vector<double> starFlux(NUM_STATE);
	Vector<double> prim(NUM_STATE);
	
	prim = getPrim(q);
	
	rho  = prim[0];
	vx   = prim[1];
    vy   = prim[2];
	p    = prim[3];
    Y    = prim[4];
    ener = q[3];
	
    if (dir == 0){  // x-direction flux
        temp = rho*(s-vx)/(s-sStar);
        starFlux[0] = f[0]+s*(temp-q[0]);
        starFlux[1] = f[1]+s*(sStar*temp-q[1]);
        starFlux[2] = f[2]+s*(vy*temp-q[2]);
        starFlux[3] = f[3]+s*(temp*((ener/rho)+(sStar-vx)*(sStar+p/(rho*(s-vx))))-q[3]);
        starFlux[4] = f[4]+s*(temp*Y-q[4]);
    }
    else {  // y-direction flux
        temp = rho*(s-vy)/(s-sStar);
        starFlux[0] = f[0]+s*(temp-q[0]);
        starFlux[1] = f[1]+s*(vx*temp-q[1]);
        starFlux[2] = f[2]+s*(sStar*temp-q[2]);
        starFlux[3] = f[3]+s*(temp*((ener/rho)+(sStar-vy)*(sStar+p/(rho*(s-vy))))-q[3]);
        starFlux[4] = f[4]+s*(temp*Y-q[4]);
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
    }

    else if((sL <= 0)&&(sR > 0))
    {
        fluxvals = fstar;
    }
    
    else //(sR < 0)
    {	
        fluxvals = fR;
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
    double Y    = prim[4];
    double ener = q[3];

    
    if (dir == 0){
        EulerFlux[0] = rho*vx;
        EulerFlux[1] = rho*vx*vx+p;
        EulerFlux[2] = rho*vx*vy;
        EulerFlux[3] = (ener+p)*vx;
        EulerFlux[4] = rho*Y*vx;
    }
    else {
        EulerFlux[0] = rho*vy;
        EulerFlux[1] = rho*vx*vy;
        EulerFlux[2] = rho*vy*vy+p;
        EulerFlux[3] = (ener+p)*vy;
        EulerFlux[4] = rho*Y*vy;
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
    double rhoY  = q[4];

    prim[0] = rho;
    prim[1] = rhovx/rho;
    prim[2] = rhovy/rho;
    prim[3] = pressure(rho,rhovx/rho,rhovy/rho,ener);
    prim[4] = rhoY/rho;

    return prim;

}

Vector<double> getCons(const Vector<double>& prim){

    // Convert primitive -> conserved variables

    
    Vector<double> q(NUM_STATE);

    double rho = prim[0];
    double vx  = prim[1];
    double vy  = prim[2];
    double p   = prim[3];
    double Y   = prim[4];

    q[0] = rho;
    q[1] = rho*vx;
    q[2] = rho*vy;
    q[3] = energy(rho,vx,vy,p);
    q[4] = rho*Y;

    return q;

}

// Below are functions that convert between various primitive and conserved variables,
// primarily used when outputting values to text files for plotting.

double energy(double rho, double vx, double vy, double p){
    double v   = pow(vx*vx + vy*vy,0.5);
    double ene = p/(Gamma-1) + 0.5*rho*v*v;
    return ene;
}

double pressure(double rho, double vx, double vy, double energy){
    double v    = pow(vx*vx + vy*vy,0.5);
    double pres = (energy-0.5*rho*v*v)*(Gamma-1);
    return pres;
}

double specIntEner(double rho, double vx, double vy, double energy){
    double v       = pow(vx*vx + vy*vy,0.5);
    double specInt = (energy-0.5*rho*v*v)/(rho);
    return specInt;
}