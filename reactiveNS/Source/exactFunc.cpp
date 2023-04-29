// This file contains various functions that are used to compute the exact
// solution to Riemann problems for the one-dimensional Euler equations.


// include header files here

#include "eulerFunc.H"
#include "exactFunc.H"
#include "AMReX_Vector.H"

// include dependencies here

#include <math.h>
#include <iostream>
#include <algorithm>
#include <stdlib.h>

extern int enIC;
extern double Gamma;
extern int NUM_STATE;


void updateExact(const int nSteps, const double dx, const int dim, std::ofstream &exact, Vector <Vector<double> >& arrExact)
{

    // This function is called to compute the exact solution based on the exact Riemann solver. Uses
    // the exact Riemann solver algorithm outlined in:
    // E. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction, 2009

    int i, k;
    double sL, sR, p0, p1, pL, pR, vL, vR, rhoL, rhoR, TOL, vStar, mfluxL, mfluxR, AL, BL, AR, BR;
    double rhoLstar, aL, aLstar, sHL, sTL, rhoRstar, aR, aRstar, sHR, sTR, rhoLfan, rhoRfan;
    double lambdaLstar, lambdaRstar, lambdaLfan, lambdaRfan, lambdaL, lambdaR;
    double x0,xEnd,tEnd,xJump;
    int iter;

    amrex::Vector<double> RPLeftRight(12);
    amrex::Vector<double> uExact(NUM_STATE);
    amrex::Vector<double> qL(NUM_STATE);
    amrex::Vector<double> qR(NUM_STATE);
    amrex::Vector<double> primL(NUM_STATE);
    amrex::Vector<double> primR(NUM_STATE);
    amrex::Vector<double> primBound(NUM_STATE);

    // Retrive the Riemann problem based on what initial condition has been selected in the inputs file.
    RPLeftRight = setIC(dim);
    
    qL[0] = RPLeftRight[0];
    qL[1] = RPLeftRight[1];
    qL[2] = RPLeftRight[2];
    qL[3] = RPLeftRight[3];
    qL[4] = RPLeftRight[4];
    qR[0] = RPLeftRight[5];
    qR[1] = RPLeftRight[6];
    qR[2] = RPLeftRight[7];
    qR[3] = RPLeftRight[8];
    qR[4] = RPLeftRight[9];
    xJump = RPLeftRight[10];
    x0    = RPLeftRight[11];
    xEnd  = RPLeftRight[12];
    tEnd  = RPLeftRight[13];

    primL = getPrim(qL);
    primR = getPrim(qR);

    // We sweep through the entire domain, and sample the solution at each point based on the number of 
    // cells defined as one of the function parameters, 'nSteps'.
    for (i = 2; i < (nSteps + 2); i++)
    {
        double error = 1;
        double x = x0+(i-1.5)*dx;
        double xOverT = (x-xJump)/tEnd;
        
        rhoL = primL[0];
        rhoR = primR[0];
        vL = primL[1];
        vR = primR[1];
        pL = primL[3];
        pR = primR[3];
        lambdaL = primL[4];
        lambdaR = primR[4];
        
        // Pre-calculate values/parameters
        AL = 2 / ((Gamma + 1) * rhoL);
        BL = pL * (Gamma - 1) / (Gamma + 1);
        AR = 2 / ((Gamma + 1) * rhoR);
        BR = pR * (Gamma - 1) / (Gamma + 1);
        aL = soundSpeed(pL, rhoL);
        aR = soundSpeed(pR, rhoR);

        // Guessed pressure value
        TOL = 1e-9;
        iter = 0;
        p0 = 0.5 * (pL + pR) - (1 / 8) * (vR - vL) * (rhoL + rhoR) * (aL + aR);
        p0 = std::min(TOL, p0);

        // Iterate on pressure
        while (error > TOL)
        {
            p1 = p0 - f(p0, pL, pR, aL, aR, vL, vR, AL, BL, AR, BR, Gamma) / fp(p0, pL, pR, aL, aR, rhoL, rhoR, AL, BL, AR, BR, Gamma);
            error = fabs(p1 - p0) / (0.5 * fabs(p1 + p0));
            p0 = std::max(0.0, p1);
            iter += 1;
            if (iter > 100)
            {
                std::cout << "no convergence" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        mfluxL = sqrt((p0 + BL) / AL);
        mfluxR = sqrt((p0 + BR) / AR);
        vStar = 0.5 * (vR + vL) + 0.5 * (fR(p0, pR, aR, AR, BR, Gamma) - fL(p0, pL, aL, AL, BL, Gamma));

        // Sample the solution
        // We have the star-state pressure and velocity, all we need is the star-state density 
        // to find the conservative flux at the interface

        // Calculate some constant values here
        double g1, g2, g4, g5, g7, g8;
        int pattern;

        g1 = p0 / pL;
        g2 = (Gamma - 1) / (Gamma + 1);
        g4 = (Gamma - 1) / (2 * Gamma);
        g5 = 2 / (Gamma + 1);
        g7 = p0 / pR;
        g8 = 2 / (Gamma - 1);

        // We have 4 possible wave patterns:

        // left shock, right shock
        if (p0 >= pL && p0 >= pR)
        {
            pattern = 1;
        }
        // left shock, right rarefaction
        else if (p0 >= pL && p0 < pR)
        {
            pattern = 2;
        }
        // left rarefaction, right rarefaction
        else if (p0 < pL && p0 < pR)
        {
            pattern = 3;
        }
        // left rarefaction, right shock
        else
        {
            pattern = 4;
        }

        // if we have a left shock wave,
        if (pattern == 1 || pattern == 2)
        {
            rhoLstar = rhoL * (g1 + g2) / (g1 * g2 + 1);
            lambdaLstar = lambdaL * (g1 + g2) / (g1 * g2 + 1);
            sL = vL - mfluxL / rhoL;
        }
        // if we have a left rarefaction,
        else
        {
            rhoLfan = rhoL * pow(g1, 1 / Gamma);
            lambdaLfan = lambdaL * pow(g1, 1 / Gamma);
            aLstar = aL * pow(g1, g4);
            sHL = vL - aL;
            sTL = vStar - aLstar;
        }
        // if we have a right shock wave,
        if (pattern == 1 || pattern == 4)
        {
            rhoRstar = rhoR * (g7 + g2) / (g7 * g2 + 1);
            lambdaRstar = lambdaR * (g7 + g2) / (g7 * g2 + 1);
            sR = vR + mfluxR / rhoR;
        }
        // if we have a right rarefaction wave,
        else
        {
            rhoRfan = rhoR * pow(g7, 1 / Gamma);
            lambdaRfan = lambdaR * pow(g7, 1 / Gamma);
            aRstar = aR * pow(g7, g4);
            sHR = vR + aR;
            sTR = vStar + aRstar;
        }        

        switch (pattern)
        {
            case 1: // shock-contact-shock
            {
                if (xOverT < sL)
                {
                    primBound = primL;
                }
                else if (sL <= xOverT && xOverT < vStar)
                {
                    primBound[0] = rhoLstar;
                    primBound[1] = vStar;
                    primBound[3] = p0;
                    primBound[4] = lambdaLstar;
                }
                else if (vStar <= xOverT && xOverT < sR)
                {
                    primBound[0] = rhoRstar;
                    primBound[1] = vStar;
                    primBound[3] = p0;
                    primBound[4] = lambdaRstar;
                }
                else if (xOverT >= sR)
                {
                    primBound = primR;
                }
                break;
            }
            case 2: // shock-contact-rarefaction
            {
                if (xOverT < sL)
                {
                    primBound = primL;
                }
                else if (sL <= xOverT && xOverT < vStar)
                {
                    primBound[0] = rhoLstar;
                    primBound[1] = vStar;
                    primBound[3] = p0;
                    primBound[4] = lambdaLstar;
                }
                else if (vStar <= xOverT && xOverT < sTR)
                {
                    primBound[0] = rhoRfan;
                    primBound[1] = vStar;
                    primBound[3] = p0;
                    primBound[4] = lambdaLfan;
                }
                else if (sTR <= xOverT && xOverT < sHR)
                {
                    primBound[0] = rhoR * pow(g5 - (g2 / aR) * (vR - xOverT), g8);
                    primBound[1] = g5 * (-aR + (1 / g8) * vR + xOverT);
                    primBound[3] = pR * pow(g5 - (g2 / aR) * (vR - xOverT), g8 * Gamma);
                    primBound[4] = lambdaL * pow(g5 - (g2 / aR) * (vR - xOverT), g8);
                }
                else if (xOverT >= sHR)
                {
                    primBound = primR;
                }
                break;
            }
            case 3: // rarefaction-contact-rarefaction
            {
                if (xOverT < sHL)
                {
                    primBound = primL;
                }
                else if (sHL <= xOverT && xOverT < sTL)
                {
                    primBound[0] = rhoL * pow(g5 + (g2 / aL) * (vL - xOverT), g8);
                    primBound[1] = g5 * (aL + (1 / g8) * vL + xOverT);
                    primBound[3] = pL * pow(g5 + (g2 / aL) * (vL - xOverT), g8 * Gamma);
                    primBound[4] = lambdaL * pow(g5 + (g2 / aL) * (vL - xOverT), g8);
                }
                else if (sTL <= xOverT && xOverT < vStar)
                {
                    primBound[0] = rhoLfan;
                    primBound[1] = vStar;
                    primBound[3] = p0;
                    primBound[4] = lambdaLfan;
                }
                else if (vStar <= xOverT && xOverT < sTR)
                {
                    primBound[0] = rhoRfan;
                    primBound[1] = vStar;
                    primBound[3] = p0;
                    primBound[4] = lambdaRfan;
                }
                else if (sTR <= xOverT && xOverT < sHR)
                {
                    primBound[0] = rhoR * pow(g5 - (g2 / aR) * (vR - xOverT), g8);
                    primBound[1] = g5 * (-aR + (1 / g8) * vR + xOverT);
                    primBound[3] = pR * pow(g5 - (g2 / aR) * (vR - xOverT), g8 * Gamma);
                    primBound[4] = lambdaR * pow(g5 - (g2 / aR) * (vR - xOverT), g8 * Gamma);

                }
                else if (xOverT >= sHR)
                {
                    primBound = primR;
                }
                break;
            }
            case 4: // rarefaction-contact-shock
            {
                if (xOverT < sHL)
                {
                    primBound = primL;
                }
                else if (sHL <= xOverT && xOverT < sTL)
                {
                    primBound[0] = rhoL * pow(g5 + (g2 / aL) * (vL - xOverT), g8);
                    primBound[1] = g5 * (aL + (1 / g8) * vL + xOverT);
                    primBound[3] = pL * pow(g5 + (g2 / aL) * (vL - xOverT), g8 * Gamma);
                    primBound[4] = lambdaL * pow(g5 + (g2 / aL) * (vL - xOverT), g8);
                }
                else if (sTL <= xOverT && xOverT < vStar)
                {
                    primBound[0] = rhoLfan;
                    primBound[1] = vStar;
                    primBound[3] = p0;
                    primBound[4] = lambdaLfan;
                }
                else if (vStar <= xOverT && xOverT < sR)
                {
                    primBound[0] = rhoRstar;
                    primBound[1] = vStar;
                    primBound[3] = p0;
                    primBound[4] = lambdaRstar;
                }
                else if (xOverT >= sR)
                {
                    primBound = primR;
                }
                break;
            }
        }

        double rhoExact = primBound[0];
        double vxExact  = primBound[1];
        double vyExact  = primBound[2];
        double pExact   = primBound[3];
        double lambdaExact = primBound[4];

        for (int h=0; h < NUM_STATE; h++){
            arrExact[i-2][h] = primBound[h];
        }

        // print out relevant quantities for plotting
        double energyExact = energy(rhoExact,vxExact,vyExact,pExact);
        double epsExact    = specIntEner(rhoExact,vxExact,vyExact,energyExact);
        exact << x << " " << rhoExact << " " << vxExact << " " << pExact << " " << epsExact << " " << lambdaExact << std::endl;
    }

}

// Below are the various function calls required when solving for the pressure iteratively.

double f(const double &p0, const double &pL, const double &pR, const double &aL, const double &aR, const double &vL, const double &vR,
         const double &AL, const double &BL, const double &AR, const double &BR, const double &Gamma)
{

    double val = fL(p0, pL, aL, AL, BL, Gamma) + fR(p0, pR, aR, AR, BR, Gamma) + vR - vL;
    return val;
}

double fL(const double &p0, const double &pL, const double &aL, const double &AL, const double &BL, const double &Gamma)
{

    double val;
    if (p0 > pL)
    {
        val = (p0 - pL) * sqrt(AL / (p0 + BL));
    }
    else
    {
        val = ((2 * aL) / (Gamma - 1)) * (pow((p0 / pL), (Gamma - 1) / (2 * Gamma)) - 1);
    }
    return val;
}

double fR(const double &p0, const double &pR, const double &aR, const double &AR, const double &BR, const double &Gamma)
{

    double val;
    if (p0 > pR)
    {
        val = (p0 - pR) * sqrt(AR / (p0 + BR));
    }
    else
    {
        val = ((2 * aR) / (Gamma - 1)) * (pow((p0 / pR), (Gamma - 1) / (2 * Gamma)) - 1);
    }
    return val;
}

double fp(const double &p0, const double &pL, const double &pR, const double &aL, const double &aR, const double &rhoL, const double &rhoR,
          const double &AL, const double &BL, const double &AR, const double &BR, const double &Gamma)
{

    double val, fLp, fRp;
    if (p0 > pL)
    {
        fLp = sqrt(AL / (BL + p0)) * (1 - (p0 - pL) / (2 * (BL + p0)));
    }
    else
    {
        fLp = (1 / (rhoL * aL)) * pow(p0 / pL, (-(Gamma + 1) / (2 * Gamma)));
    }
    if (p0 > pR)
    {
        fRp = sqrt(AR / (BR + p0)) * (1 - (p0 - pR) / (2 * (BR + p0)));
    }
    else
    {
        fRp = (1 / (rhoR * aR)) * pow(p0 / pR, (-(Gamma + 1) / (2 * Gamma)));
    }
    val = fLp + fRp;
    return val;
}
