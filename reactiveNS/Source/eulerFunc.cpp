#include "eulerFunc.H"
#include "AMReX_Vector.H"
#include "AMReX_REAL.H"

#include <math.h>
#include <iostream>
#include <algorithm> 

extern int enIC;
extern double Gamma;
extern int NUM_STATE;

amrex::Vector<double> setIC(const int dim) {

    double rhoL, vxL, vyL, pL, rhoR, vxR, vyR, pR;
    double x0,xEnd,xDisc,tEnd;
    amrex::Vector<double> RPLeftRight(12);

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

    RPLeftRight[0] = rhoL;
    RPLeftRight[1] = rhoL*vxL;
    RPLeftRight[2] = rhoL*vyL;
    RPLeftRight[3] = energy(rhoL,vxL,vyL,pL);
    RPLeftRight[4] = rhoR;
    RPLeftRight[5] = rhoR*vxR;
    RPLeftRight[6] = rhoR*vyR;
    RPLeftRight[7] = energy(rhoR,vxR,vyR,pR);
    RPLeftRight[8] = xDisc;
    RPLeftRight[9] = x0;
    RPLeftRight[10] = xEnd;
    RPLeftRight[11] = tEnd;

    return RPLeftRight;    
}

void getStopTime(int enIC, amrex::Real& stop_time){
    switch (enIC) {
        case 1: //toro test 1
        {
            stop_time=0.45;
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
    }
}


double soundSpeed(const double& p, const double& rho){
    double a = sqrt(fabs(Gamma*p/rho));
    return a;
}

void wavespeedEstimate(const amrex::Vector<double>& qL, const amrex::Vector<double>& qR,\
                       double& sL, double& sR, double& sStar, const int dir){

    // This function computes wavespeed estimates for the HLLC solver, retrieved via pass by reference variables sL, sR, and sStar.

    double pL,pR,rhoL,rhoR,vxL,vxR,vyL,vyR,vL,vR,aL,aR;
    double rhobar, abar, pPVRS, pStar, qqL, qqR;    
    amrex::Vector<double> primL(NUM_STATE);
    amrex::Vector<double> primR(NUM_STATE);

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

    // Here, we use a pressure-based wave speed estimate proposed by Toro et al. (1994, Shock Waves)
    // The star-region pressure estimate follows the form of the PVRS approximate Riemann solver by Toro (1991, Proc. Roy. Soc. London)

    rhobar = 0.5*(rhoL+rhoR);
    abar   = 0.5*(aL+aR);
    pPVRS  = 0.5*(pL+pR)-0.5*(vR-vL)*rhobar*abar;  
    pStar  = std::max(0.0,pPVRS);

    if (pStar <= pL){
        qqL = 1;
    }
    else {
        qqL = sqrt(std::abs(1+((Gamma+1)/(2*Gamma))*((pStar/pL) - 1)));
    }
    if (pStar <= pR){
        qqR = 1;
    }
    else {
        qqR = sqrt(std::abs(1+((Gamma+1)/(2*Gamma))*((pStar/pR) - 1)));
    }

    sL = vL - aL*qqL;
    sR = vR + aR*qqR;
    sStar = (pR-pL+rhoL*vL*(sL-vL)-rhoR*vR*(sR-vR))/(rhoL*(sL-vL)-rhoR*(sR-vR));


}

amrex::Vector<double> HLLCflux(const amrex::Vector<double>& qL, const amrex::Vector<double>& qR, const int& dir){

    // This function calculates the HLLC flux based on left and right states, qL and qR, respectively
    
    double sL,sR,sStar;
    amrex::Vector<double> fL(NUM_STATE);
    amrex::Vector<double> fR(NUM_STATE);
    amrex::Vector<double> fLstar(NUM_STATE);
    amrex::Vector<double> fRstar(NUM_STATE);
    amrex::Vector<double> fluxvals(NUM_STATE);

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

amrex::Vector<double> HLLCstarFlux(const amrex::Vector<double>& f, const amrex::Vector<double>& q,\
                                 const double& s, const double& sStar, const int& dir){

    // Calculates HLLC star-state flux for HLLC approximate Riemann solver.

	double rho,vx,vy,p,ener,temp;
    amrex::Vector<double> starFlux(NUM_STATE);
	amrex::Vector<double> prim(NUM_STATE);
	
	prim = getPrim(q);
	
	rho  = prim[0];
	vx   = prim[1];
    vy   = prim[2];
	p    = prim[3];
    ener = q[3];
	
    if (dir == 0){  // x-direction flux
        temp = rho*(s-vx)/(s-sStar);
        starFlux[0] = f[0]+s*(temp-q[0]);
        starFlux[1] = f[1]+s*(sStar*temp-q[1]);
        starFlux[2] = f[2]+s*(vy*temp-q[2]);
        starFlux[3] = f[3]+s*(temp*((ener/rho)+(sStar-vx)*(sStar+p/(rho*(s-vx))))-q[3]);
    }
    else {  // y-direction flux
        temp = rho*(s-vy)/(s-sStar);
        starFlux[0] = f[0]+s*(temp-q[0]);
        starFlux[1] = f[1]+s*(vx*temp-q[1]);
        starFlux[2] = f[2]+s*(sStar*temp-q[2]);
        starFlux[3] = f[3]+s*(temp*((ener/rho)+(sStar-vy)*(sStar+p/(rho*(s-vy))))-q[3]);
    }
	
	
	
	return starFlux;
	
}

amrex::Vector<double> getEulerFlux(const amrex::Vector<double>& q, const int& dir){

    // Takes in conservative variable vector and directional input, returns conservative flux vector
    // based on conservation form of Euler equations.
    
    amrex::Vector<double> EulerFlux(NUM_STATE);
    amrex::Vector<double> prim(NUM_STATE);
    prim = getPrim(q);

    double rho  = prim[0];
    double vx   = prim[1];
    double vy   = prim[2];
    double p    = prim[3];
    double ener = q[3];

    
    if (dir == 0){
        EulerFlux[0] = rho*vx;
        EulerFlux[1] = rho*vx*vx+p;
        EulerFlux[2] = rho*vx*vy;
        EulerFlux[3] = (ener+p)*vx;
    }
    else {
        EulerFlux[0] = rho*vy;
        EulerFlux[1] = rho*vx*vy;
        EulerFlux[2] = rho*vy*vy+p;
        EulerFlux[3] = (ener+p)*vy;
    }
    
    return EulerFlux;

}

amrex::Vector<double> getPrim(const amrex::Vector<double>& q){

    // Convert conserved -> primitive variables
    
    amrex::Vector<double> prim(NUM_STATE);

    double rho   = q[0];
    double rhovx = q[1];
    double rhovy = q[2];
    double ener  = q[3];

    prim[0] = rho;
    prim[1] = rhovx/rho;
    prim[2] = rhovy/rho;
    prim[3] = pressure(rho,rhovx/rho,rhovy/rho,ener);

    return prim;

}

amrex::Vector<double> getCons(const amrex::Vector<double>& prim){

    // Convert primitive -> conserved variables
    
    amrex::Vector<double> q(NUM_STATE);

    double rho = prim[0];
    double vx  = prim[1];
    double vy  = prim[2];
    double p   = prim[3];

    q[0] = rho;
    q[1] = rho*vx;
    q[2] = rho*vy;
    q[3] = energy(rho,vx,vy,p);

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