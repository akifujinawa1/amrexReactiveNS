#include "eulerFunc.H"
#include "diffusionFunc.H"
#include "thermoTransport.H"
#include "AmrLevelAdv.H"
#include "AMReX_Vector.H"
#include "AMReX_Array.H"
#include "AMReX_REAL.H"
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>

#include <cmath>
#include <iostream>
#include <algorithm> 

extern int enIC;
extern double Gamma;
extern int NUM_STATE;
extern const int spacedim;
extern double R;
extern double M;
extern double Pr;
extern double one_atm_Pa;
extern double Dpre;

using namespace amrex;

// Use this file to write functions required for diffusive flux calculations

void updateViscous(MultiFab& Sborder, Array<MultiFab, SpaceDim>& fluxes, Vector<double> &qL, Vector<double> &qR, \
                   Vector<double> &qLlo, Vector<double> &qRlo, Vector<double> &qLhi, Vector<double> &qRhi, \
                   Vector<double> &viscSlice, const int &d, const double &dt, \
                   const double& dx, const double& dy, const int& SpaceDim, const int &viscous){
    switch (viscous)
    {
        case 0: 
        {
            break;
        }
        case 1: // Use central
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
                const auto& fluxArrVisc = fluxes[d].array(mfi);
                const auto& fluxArrViscY = fluxes[d].array(mfi);

                // If slopelimiting has been enabled in the settings file, we run through the following set of loops
                // to compute boundary-reconstructed values -2023W2
                if (d == 0)
                {
                    for(int k = lo.z; k <= hi.z; k++)
                    {
                        for(int j = lo.y; j <= hi.y; j++)
                        {
                            for(int i = lo.x; i <= hi.x+iOffset; i++) // at each i we calculate the left interface flux 
                            {           
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    qL[h]   = arr(i-1,j,k,h);   // cell to the left of cell i
                                    qR[h]   = arr(i  ,j,k,h);   // cell i 
                                    if (amrex::SpaceDim == 2)
                                    {
                                        qLlo[h] = arr(i-1,j-1,k,h);   // cell to the left lower diagonal of cell i
                                        qRlo[h] = arr(i  ,j-1,k,h);   // cell to the right lower diagonal of cell i 
                                        qLhi[h] = arr(i-1,j+1,k,h);   // cell to the left higher diagonal of cell i
                                        qRhi[h] = arr(i  ,j+1,k,h);   // cell to the right higher diagonal of cell i 
                                    }
                                }
                                
                                if (amrex::SpaceDim == 1){ // x-direction viscous flux calculation in 1-D
                                    getViscFlux1D(viscSlice,qL,qR,dx);
                                    for(int h = 0; h < NUM_STATE; h++)
                                    {
                                        fluxArrVisc(i,j,k,h) = viscSlice[h]; // this is the viscous flux function for the left interface of cell i
                                    }
                                }
                                else { // x-direction viscous flux calculation if domain is 2-D
                                    getViscFlux2D(viscSlice,qL,qR,qLlo,qRlo,qLhi,qRhi,d,dx,dy);
                                    for(int h = 0; h < NUM_STATE; h++) {
                                    fluxArrVisc(i,j,k,h) = viscSlice[h];
                                    }
                                }
                            }
                        }
                    }
                }
                else // if solving for y-direction viscous fluxes (only applicable to 2-D)
                { 
                    for(int k = lo.z; k <= hi.z; k++)
                    {
                        for(int i = lo.x; i <= hi.x; i++)
                        {
                            for(int j = lo.y; j <= hi.y+jOffset; j++)
                            {           
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    qL[h]   = arr(i,  j-1,k,h);   // cell below cell i,j
                                    qR[h]   = arr(i,  j  ,k,h);   // cell i,j 
                                    qLlo[h] = arr(i-1,j-1,k,h);   // cell to the left lower diagonal of cell i
                                    qRlo[h] = arr(i-1,j  ,k,h);   // cell to the right lower diagonal of cell i 
                                    qLhi[h] = arr(i+1,j-1,k,h);   // cell to the left higher diagonal of cell i
                                    qRhi[h] = arr(i+1,j  ,k,h);   // cell to the right higher diagonal of cell i 
                                }
                                getViscFlux2D(viscSlice,qL,qR,qLlo,qRlo,qLhi,qRhi,d,dx,dy);
                                for(int h = 0; h < NUM_STATE; h++) {
                                    fluxArrViscY(i,j,k,h) = viscSlice[h];
                                }
                            }
                        }
                    }
                }
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
                                    arr(i,j,k,h) = arr(i,j,k,h) + (dt / dx) * (fluxArrVisc(i,j,k,h) - fluxArrVisc(i+iOffset,j,k,h));
                                    if (arr(i,j,k,h) != arr(i,j,k,h)){
                                        std::cout << "Nan found in diffusion calculation, variable h: " << h << std::endl;
                                        Abort("nan found in diffusion calculation");
                                    }
                                }
                                if (i == 7){
                                    std::cout << "Cell 7, New rho rhou rhov e o2 n2\n" << arr(i,j,k,0) << " " << arr(i,j,k,1) << \
                                    " " << arr(i,j,k,2) << " " << arr(i,j,k,3) << " " << arr(i,j,k,4) << " " << arr(i,j,k,5) << std::endl;
                                }
                                if (i == 8){
                                    std::cout << "Cell 8, New rho rhou rhov e o2 n2\n" << arr(i,j,k,0) << " " << arr(i,j,k,1) << \
                                    " " << arr(i,j,k,2) << " " << arr(i,j,k,3) << " " << arr(i,j,k,4) << " " << arr(i,j,k,5) << std::endl;
                                }
                                if (i == 9){
                                    std::cout << "Cell 9, New rho rhou rhov e o2 n2\n" << arr(i,j,k,0) << " " << arr(i,j,k,1) << \
                                    " " << arr(i,j,k,2) << " " << arr(i,j,k,3) << " " << arr(i,j,k,4) << " " << arr(i,j,k,5) << std::endl;
                                }
                            }
                            else { // y=direction update
                                for(int h = 0; h < NUM_STATE; h++)
                                {
                                    arr(i,j,k,h) = arr(i,j,k,h) - (dt / dy) * (fluxArrViscY(i, j+jOffset, k,h) - fluxArrViscY(i,j,k,h));
                                }
                            }
                           
                        }
                    }
                }
            } // close loop for patches
        } // close case 1 central
        break;
    } // close viscous switch case    

}

void getViscFlux1D(Vector<double>& viscSlice, const Vector<double>& qL,\
                   const Vector<double>& qR, const double& dx){

    double rhoL, uL, vL, enerL, epsL, TL, YO2L, YN2L, pL, rhoR, uR, vR, enerR, epsR, TR, YO2R, YN2R, pR;
    double dudx, dTdx, dYO2dx, dYN2dx, uAvg, TAvg, pAvg, rhoAvg, YO2Avg, YN2Avg;
    double mu_avg, k_avg;

    Vector<double> mixDiffCoeffs(gases::ncomps);
    Vector<double> primL(NUM_STATE);
    Vector<double> primR(NUM_STATE);

    primL   = getPrim(qL);
    primR   = getPrim(qR);
    
    rhoL  = primL[0];
    uL    = primL[1];
    vL    = primL[2];
    enerL = qL[3];
    pL    = primL[3];
    YO2L  = primL[4];
    YN2L  = primL[5];
    TL    = Tg(rhoL,uL,vL,YO2L,YN2L,enerL);
    
    rhoR  = primR[0];
    uR    = primR[1];
    vR    = primR[2];
    enerR = qR[3];
    pR    = primR[3];
    YO2R  = primR[4];
    YN2R  = primR[5];
    TR    = Tg(rhoR,uR,vR,YO2R,YN2R,enerR);

    dudx = (uL-uR)/dx;
    dTdx = (TL-TR)/dx;
    dYO2dx = (YO2L-YO2R)/dx;
    dYN2dx = (YN2L-YN2R)/dx;
    uAvg = (uR+uL)/2.0;
    TAvg = (TR+TL)/2.0;
    pAvg = (pR+pL)/2.0;
    rhoAvg = (rhoR+rhoL)/2.0;
    YO2Avg = (YO2R+YO2L)/2.0;
    YN2Avg = (YN2R+YN2L)/2.0;
    mixDiffCoeffs = getMixDiffCoeffs(TAvg,pAvg,YO2Avg,YN2Avg,0.0,0.0);
    mu_avg = muMix_O2N2(muO2(TAvg),muN2(TAvg),YO2Avg,YN2Avg); // consider only the O2-N2 gas mixture
    k_avg  = kMix_O2N2(kO2(TAvg),kN2(TAvg),YO2Avg,YN2Avg);

    // std::cout << "rho, T, p, YO2, YN2: " << rhoAvg << " " << TAvg << " " << pAvg << " " << YO2Avg << " " << YN2Avg << std::endl; 
    // std::cout << "mu, k, DO2, DN2: " << mu_avg << " " << k_avg << " " << mixDiffCoeffs[gases::O2] << " " << mixDiffCoeffs[gases::N2] << std::endl;

    viscSlice[0] = 0;
    viscSlice[1] = (4.0/3.0)*mu_avg*dudx;
    viscSlice[2] = 0;
    viscSlice[3] = uAvg*viscSlice[1] + k_avg*dTdx;
    viscSlice[4] = rhoAvg*mixDiffCoeffs[gases::O2]*dYO2dx;
    viscSlice[5] = rhoAvg*mixDiffCoeffs[gases::N2]*dYN2dx;

}

void getViscFlux2D(Vector<double>& viscSlice, const Vector<double>& qL, const Vector<double>& qR, \
                   const Vector<double>& qLlo, const Vector<double>& qRlo, const Vector<double>& qLhi, \
                   const Vector<double>& qRhi, const int& d, const double& dx, const double& dy){

    double cV = 3.5*287;
    double rhoL, uL, vL, enerL, epsL, TL, pL, YO2L, YN2L, rhoR, uR, vR, enerR, epsR, TR, pR, YO2R, YN2R;
    double rhoLlo, uLlo, vLlo, enerLlo, epsLlo, TLlo, pLlo, YO2Llo, YN2Llo;
    double rhoRlo, uRlo, vRlo, enerRlo, epsRlo, TRlo, pRlo, YO2Rlo, YN2Rlo;
    double rhoLhi, uLhi, vLhi, enerLhi, epsLhi, TLhi, pLhi, YO2Lhi, YN2Lhi;
    double rhoRhi, uRhi, vRhi, enerRhi, epsRhi, TRhi, pRhi, YO2Rhi, YN2Rhi;
    double dudx, dudy, dvdx, dvdy, dTdx, dTdy, dYO2dx, dYN2dx, dYO2dy, dYN2dy, uAvg, vAvg, TAvg, rhoAvg, pAvg, YO2Avg, YN2Avg;
    // double mu = 3e-5;
    // double k  = 1.4;
    double mu_avg, k_avg;
    Vector<double> mixDiffCoeffs(gases::ncomps);

    Vector<double> primL(NUM_STATE);
    Vector<double> primR(NUM_STATE);
    Vector<double> primLlo(NUM_STATE);
    Vector<double> primRlo(NUM_STATE);
    Vector<double> primLhi(NUM_STATE);
    Vector<double> primRhi(NUM_STATE);

    primL   = getPrim(qL);
    primR   = getPrim(qR);
    primLlo = getPrim(qLlo);
    primRlo = getPrim(qRlo);
    primLhi = getPrim(qLhi);
    primRhi = getPrim(qRhi);
    
    rhoL  = primL[0];
    uL    = primL[1];
    vL    = primL[2];
    enerL = qL[3];
    pL    = primL[3];
    YO2L  = primL[4];
    YN2L  = primL[5];
    TL    = Tg(rhoL,uL,vL,YO2L,YN2L,enerL);
    
    rhoR  = primR[0];
    uR    = primR[1];
    vR    = primR[2];
    enerR = qR[3];
    pR    = primR[3];
    YO2R  = primR[4];
    YN2R  = primR[5];
    TR    = Tg(rhoR,uR,vR,YO2R,YN2R,enerR);

    rhoLlo  = primLlo[0];
    uLlo    = primLlo[1];
    vLlo    = primLlo[2];
    enerLlo = qLlo[3];
    pLlo    = primLlo[3];
    YO2Llo  = primLlo[4];
    YN2Llo  = primLlo[5];
    TLlo    = Tg(rhoLlo,uLlo,vLlo,YO2Llo,YN2Llo,enerLlo);
    
    rhoRlo  = primRlo[0];
    uRlo    = primRlo[1];
    vRlo    = primRlo[2];
    enerRlo = qRlo[3];
    pRlo    = primRlo[3];
    YO2Rlo  = primRlo[4];
    YN2Rlo  = primRlo[5];
    TRlo    = Tg(rhoRlo,uRlo,vRlo,YO2Rlo,YN2Rlo,enerRlo);

    rhoLhi  = primLhi[0];
    uLhi    = primLhi[1];
    vLhi    = primLhi[2];
    enerLhi = qLhi[3];
    pLhi    = primLhi[3];
    YO2Lhi  = primLhi[4];
    YN2Lhi  = primLhi[5];
    TLhi    = Tg(rhoLhi,uLhi,vLhi,YO2Lhi,YN2Lhi,enerLhi);
    
    rhoRhi  = primRhi[0];
    uRhi    = primRhi[1];
    vRhi    = primRhi[2];
    enerRhi = qRhi[3];
    pRhi    = primRhi[3];
    YO2Rhi  = primRhi[4];
    YN2Rhi  = primRhi[5];
    TRhi    = Tg(rhoRhi,uRhi,vRhi,YO2Rhi,YN2Rhi,enerRhi);

    if (d == 0) { // finding viscous fluxes for x-direction update

        dudx = (uR-uL)/dx;
        dvdx = (uR-uL)/dx;
        dTdx = (TR-TL)/dx;
        dYO2dx = (YO2R-YO2L)/dx;
        dYN2dx = (YN2R-YN2L)/dx;
        dudy = 0.5*((uRhi-uRlo)/(2*dy)+(uLhi-uLlo)/(2.0*dy));
        dvdy = 0.5*((vRhi-vRlo)/(2*dy)+(vLhi-vLlo)/(2.0*dy));
        uAvg = (uR+uL)/2.0;
        vAvg = (uR+uL)/2.0;
        TAvg = (TR+TL)/2.0;
        pAvg = (pR+pL)/2.0;
        rhoAvg = (rhoR+rhoL)/2.0;
        YO2Avg = (YO2R+YO2L)/2.0;
        YN2Avg = (YN2R+YN2L)/2.0;
        mixDiffCoeffs = getMixDiffCoeffs(TAvg,pAvg,YO2Avg,YN2Avg,0.0,0.0);
        mu_avg = muMix_O2N2(muO2(TAvg),muN2(TAvg),YO2Avg,YN2Avg);
        k_avg  = kMix_O2N2(kO2(TAvg),kN2(TAvg),YO2Avg,YN2Avg);

        std::cout << "p: " << pAvg << ", T: " << TAvg << ", D: " << mixDiffCoeffs[gases::O2] << ", mu: " << mu_avg << ", k: " << k_avg << std::endl;

        viscSlice[0] = 0;
        viscSlice[1] = 2.0*mu_avg*dudx - (2.0/3.0)*mu_avg*(dudx+dvdy);
        viscSlice[2] = mu_avg*(dvdx+dudy);
        viscSlice[3] = uAvg*viscSlice[1] + vAvg*viscSlice[2] + k_avg*dTdx;
        viscSlice[4] = rhoAvg*mixDiffCoeffs[gases::O2]*dYO2dx;
        viscSlice[5] = rhoAvg*mixDiffCoeffs[gases::N2]*dYN2dx;

    }

    else { // finding viscous fluxes for y-direction update

        dudy = (uR-uL)/dy;
        dvdy = (vR-vL)/dy;
        dTdy = (TR-TL)/dy;
        dYO2dy = (YO2R-YO2L)/dy;
        dYN2dy = (YN2R-YN2L)/dy;
        dudx = 0.5*((uRhi-uRlo)/(2.0*dx)+(uLhi-uLlo)/(2.0*dx));
        dvdx = 0.5*((vRhi-vRlo)/(2.0*dx)+(vLhi-vLlo)/(2.0*dx));
        uAvg = (uR+uL)/2.0;
        vAvg = (uR+uL)/2.0;
        TAvg = (TR+TL)/2.0;
        pAvg = (pR+pL)/2.0;
        rhoAvg = (rhoR+rhoL)/2.0;
        YO2Avg = (YO2R+YO2L)/2.0;
        YN2Avg = (YN2R+YN2L)/2.0;
        mixDiffCoeffs = getMixDiffCoeffs(TAvg,pAvg,YO2Avg,YN2Avg,0.0,0.0);
        mu_avg = muMix_O2N2(muO2(TAvg),muN2(TAvg),YO2Avg,YN2Avg);
        k_avg  = kMix_O2N2(kO2(TAvg),kN2(TAvg),YO2Avg,YN2Avg);

        std::cout << "p: " << pAvg << ", T: " << TAvg << ", D: " << mixDiffCoeffs[gases::O2] << ", mu: " << mu_avg << ", k: " << k_avg << std::endl;

        viscSlice[0] = 0.0;
        viscSlice[1] = mu_avg*(dvdx+dudy);
        viscSlice[2] = 2.0*mu_avg*dvdy - (2.0/3.0)*mu_avg*(dudx+dvdy);
        viscSlice[3] = uAvg*viscSlice[1] + vAvg*viscSlice[2] + k_avg*dTdy;
        viscSlice[4] = rhoAvg*mixDiffCoeffs[gases::O2]*dYO2dy;
        viscSlice[5] = rhoAvg*mixDiffCoeffs[gases::N2]*dYN2dy;
    }
    
}


double diffusiveSpeed(const Vector<double>& qL, const Vector<double>& qR){

    double rhoL, uL, vL, enerL, epsL, TL, YO2L, YN2L, pL, rhoR, uR, vR, enerR, epsR, TR, YO2R, YN2R, pR;
    double uAvg, TAvg, rhoAvg, pAvg, YO2Avg, YN2Avg;
    double mu_avg,k_avg,cp_avg;
    double maxLeft, maxRight;

    Vector<double> mixDiffCoeffs(gases::ncomps);
    Vector<double> primL(NUM_STATE);
    Vector<double> primR(NUM_STATE);

    primL   = getPrim(qL);
    primR   = getPrim(qR);
    
    rhoL  = primL[0];
    uL    = primL[1];
    vL    = primL[2];
    enerL = qL[3];
    pL    = primL[3];
    YO2L  = primL[4];
    YN2L  = primL[5];
    TL    = Tg(rhoL,uL,vL,YO2L,YN2L,enerL);

    // std::cout << "temp: " << TL << std::endl;
    
    rhoR  = primR[0];
    uR    = primR[1];
    vR    = primR[2];
    enerR = qR[3];
    pR    = primR[3];
    YO2R  = primR[4];
    YN2R  = primR[5];
    TR    = Tg(rhoR,uR,vR,YO2R,YN2R,enerR);

    uAvg = (uR+uL)/2.0;
    TAvg = (TR+TL)/2.0;
    pAvg = (pR+pL)/2.0;
    rhoAvg = (rhoR+rhoL)/2.0;
    YO2Avg = (YO2R+YO2L)/2.0;
    YN2Avg = (YN2R+YN2L)/2.0;
    mixDiffCoeffs = getMixDiffCoeffs(TAvg,pAvg,YO2Avg,YN2Avg,0.0,0.0);
    double maxD   = std::max(mixDiffCoeffs[gases::O2],mixDiffCoeffs[gases::N2]);
    mu_avg = muMix_O2N2(muO2(TAvg),muN2(TAvg),YO2Avg,YN2Avg); // consider only the O2-N2 gas mixture
    k_avg  = kMix_O2N2(kO2(TAvg),kN2(TAvg),YO2Avg,YN2Avg);
    cp_avg = cpMix(cpO2(TAvg),cpN2(TAvg),YO2Avg,YN2Avg);
    double heat = k_avg/(rhoAvg*cp_avg);
    
    double maxSpeed = 2*std::max(std::max(mu_avg/rhoAvg,k_avg/(rhoAvg*cp_avg)),maxD);

    // if (maxSpeed == mu_avg/rhoAvg){
    //     std::cout << "momentum limiting" << std::endl;
    // }
    // else if (maxSpeed == k_avg/(rhoAvg*cp_avg)){
    //     std::cout << "heat limiting" << std::endl;
    // }
    // else {
    //     std::cout << "species limiting" << std::endl;
    // }

    if (enIC == 8){
        return 2.0*2.0e-5;
    }

    // std::cout << "mu/rho: " << mu_avg/rhoAvg << ", heat: " << heat << ", D: " << maxD << ", max speed" << maxSpeed << std::endl;

    return maxSpeed;

}

