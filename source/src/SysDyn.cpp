/*
 * SysDyn.cpp
 *
 *    VIABLAB : a numerical library for Mathematical Viability Computations
 *    Copyright (C) <2020>  <Anna DESILLES, LASTRE>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Created on: 9 december 2013
 *      Author: Anna DESILLES
 */

#include "../include/SysDyn.h"

SysDyn::SysDyn()
{

}

SysDyn::SysDyn(const systemParams &SP, int ds, const controlParams &cp, Grid *grRef)
{
    grid = grRef;
    dimS = ds;
    dimC = cp.DIMC;
    limInfC = cp.LIMINFC;
    limSupC = cp.LIMSUPC;
    nbPointsC = cp.NBPOINTSC;
    totalNbPointsC = 1;

    dimTy = cp.DIM_TY;
    limInfTy = cp.LIMINF_TY;
    limSupTy = cp.LIMSUP_TY;
    nbPointsTy = cp.NBPOINTS_TY;
    totalNbPointsTych = 1;

    image = new double[dimS];
    FXmoinsH = new double[dimS];
    xTemp = new double[dimS];
    FXplusH = new double[dimS];
    globalTimeStep = SP.globDeltat;

    localDynBounds = SP.LOCAL_DYN_BOUNDS;

    jacob = new double*[dimS];

    for (int i = 0; i < dimS; i++)
    {
        jacob[i] = new double[dimS];
    }

    MF = SP.MF;
    L = SP.LIP;

    lfunc_L = SP.L_LIP;  // constante de Lipschitz
    lfunc_MF = SP.L_MF;  // majoration de la norme de la dynamique

    if (L == 0)
    {
        L = 1.0;
    }
    if (MF == 0)
    {
        MF = 1.0;
    }

    dynSignFactor = 1.0; //forward by default
    unsigned long long int k;
    int dc;
    /*!
     *  attention! Pour l'instant, la convention par d�faut
     *  pour la grille de contr�les  est la suivante:
     *  les  points sont les centres de mailles
     *  avec cette convention on peut m�me avoit un point  d controle par axe
     *  sans problemes
     *  Donc le nombre de points=nombre d'intervalles !!!!!
     */

    if (dimC > 0)
    {
        stepC = new double[dimC];

        spdlog::debug("[System] : Control dimension {}", dimC);

        for (dc = 0; dc < dimC; dc++)
        {
            totalNbPointsC *= nbPointsC[dc];
            stepC[dc] = (limSupC[dc] - limInfC[dc]) / (nbPointsC[dc] - 1);
        }
        logVector("[System]: control discretization step :", stepC, dimC);
        logVector("[System]: control inf limits :", limInfC, dimC);
        logVector("[System]: control sup limits :", limSupC, dimC);
        controlCoords = new double*[totalNbPointsC];
        controlIntCoords = new unsigned long long int*[totalNbPointsC];
        unsigned long long int *coordsIntC = new unsigned long long int[dimC];

        for (k = 0; k < totalNbPointsC; k++)
        {
            controlCoords[k] = new double[dimC];
            controlIntCoords[k] = new unsigned long long int[dimC];
            numToIntCoords_gen(k, dimC, nbPointsC, coordsIntC);
            for (dc = 0; dc < dimC; dc++)
            {

                controlCoords[k][dc] = limInfC[dc] + stepC[dc] * coordsIntC[dc]; //+0.5*stepC[dc];
                controlIntCoords[k][dc] = coordsIntC[dc];
            }
        }
    }
    else
    {
        spdlog::warn("[System] : dynamic system without control");
        stepC = new double[1];
        for (dc = 0; dc < dimC; dc++)
        {
            totalNbPointsC *= nbPointsC[dc];
            stepC[dc] = (limSupC[dc] - limInfC[dc]) / (nbPointsC[dc]);
        }
        controlCoords = new double*[totalNbPointsC];
        controlIntCoords = new unsigned long long int*[totalNbPointsC];

        for (k = 0; k < totalNbPointsC; k++)
        {
            controlCoords[k] = new double[1];
            controlIntCoords[k] = new unsigned long long int[1];
        }
    }

    if (dimTy > 0)
    {

        spdlog::warn("[System] : dynamic system with tychastic control");
        spdlog::debug("[System] : Tychastic control dimansion {}", dimTy);
        stepTy = new double[dimTy];

        for (dc = 0; dc < dimTy; dc++)
        {
            totalNbPointsTych *= nbPointsTy[dc];
            stepTy[dc] = (limSupTy[dc] - limInfTy[dc]) / (nbPointsTy[dc] - 1);
        }
        logVector("[System]: tychastic control discretization step :", stepTy,
                  dimTy);
        logVector("[System]: tychastic control inf limits :", limInfTy, dimTy);
        logVector("[System]: tychastic control sup limits :", limSupTy, dimTy);

        tychCoords = new double*[totalNbPointsTych];
        tychIntCoords = new unsigned long long int*[totalNbPointsTych];
        unsigned long long int *coordsIntTy = new unsigned long long int[dimTy];
        cout << "  coords de tych \n";
        for (k = 0; k < totalNbPointsTych; k++)
        {
            tychCoords[k] = new double[dimTy];
            tychIntCoords[k] = new unsigned long long int[dimTy];
            numToIntCoords_gen(k, dimTy, nbPointsTy, coordsIntTy);
            for (dc = 0; dc < dimTy; dc++)
            {
                tychCoords[k][dc] = limInfTy[dc] + stepTy[dc] * coordsIntTy[dc];
                tychIntCoords[k][dc] = coordsIntTy[dc];
            }
        }

    }

    initializeMethods(SP);
    
    spdlog::info("[System] : Dynamic system initialization finished");
}

void SysDyn::initializeMethods(const systemParams &SP)
{
    dynamics = SP.DYNAMICS;
    dynamics_fd = SP.DYNAMICS_FD;
    dynamics_tych_fd = SP.DYNAMICS_TYCH_FD;
    dynamics_tych = SP.DYNAMICS_TYCH;

    dynType = SP.DYN_TYPE;

    if (SP.DYN_TYPE == DD)
    {
        cout << " DynSYS de type DD" << endl;
        if (SP.FD_DYN_TYPE == RETRO)
        {
            this->retroFileName = SP.RETRO_FILE_NAME;
        }
    }
    fd_dyn_type = SP.FD_DYN_TYPE;
    constraintsXU = SP.CONSTR_XU;
    constraintsXU_fd = SP.CONSTR_XU_fd;
    constraintsXV_tych = SP.CONSTR_XV_TYCH;
    constraintsX = SP.CONSTR_X;
    constraintsX_fd = SP.CONSTR_X_fd;
    controlEligibilityForTraj_fd = SP.CONTROL_ELIGIBILITY_FOR_TRAJ_fd;
    dynConstraintsForTraj = SP.DYN_CONSTR_FOR_TRAJ;
    target = SP.TARGET;
    target_fd = SP.TARGET_FD;
    
    jacobian = SP.JACOBIAN;
    jacobian_tych = SP.JACOBIAN_TYCH;

    lFunc = SP.L_FUNC;
    lFunc_tych = SP.L_FUNC_TYCH;
    lFunc_fd = SP.L_FUNC_FD;
    lFunc_tych_fd = SP.L_FUNC_TYCH_FD;
    muFunc_fd = SP.MU_FUNC_FD;
    mFunc = SP.M_FUNC;
    mFunc_tych = SP.M_FUNC_TYCH;

    computeMF = SP.COMPUTE_MF;
    computeLC = SP.COMPUTE_LC;
    isTychastic = (dimTy > 0);
    
    switch (computeMF)
    {
    case ANALYTICAL:
        SysDyn::calcul_M = &SysDyn::returnMF_local_ana;
        break;
    case ANALYTICAL_CALC:
        SysDyn::calcul_M = &SysDyn::calculMF_local_ana;
        break;
    case NUMERICAL_CALC:
        SysDyn::calcul_M = isTychastic ? &SysDyn::calculMF_local_num_tych : &SysDyn::calculMF_local_num;
        break;
    }

    switch (computeLC)
    {
    case ANALYTICAL:
        SysDyn::calcul_L =  &SysDyn::returnL_local_ana;
        break;
    case ANALYTICAL_CALC:
        SysDyn::calcul_L = isTychastic ? &SysDyn::calculL_local_ana_tych : &SysDyn::calculL_local_ana;
        break;
    case NUMERICAL_CALC:
        SysDyn::calcul_L = isTychastic ? &SysDyn::calculL_local_num_tych : &SysDyn::calculL_local_num;
        break;
    }

    discretisation = SP.SCHEME;

    if (dynType == DC)
    {
        discretisation = NO_DISCRETIZATION_SCHEME;
    }

    switch (discretisation)
    {
    case NO_DISCRETIZATION_SCHEME:
        SysDyn::discretDynamics = &SysDyn::FDiscret;
        SysDyn::discretDynamics_tych = &SysDyn::FDiscret_tych;
        break;
    case EL:
        SysDyn::discretDynamics = &SysDyn::FDiscretEuler;
        SysDyn::discretDynamics_tych = &SysDyn::FDiscretEuler_tych;
        break;
    case RK2:
        SysDyn::discretDynamics = &SysDyn::FDiscretRK2;
        SysDyn::discretDynamics_tych = &SysDyn::FDiscretRK2_tych;
        break;
    case RK4:
        SysDyn::discretDynamics = &SysDyn::FDiscretRK4;
        SysDyn::discretDynamics_tych = &SysDyn::FDiscretRK4_tych;
        break;
    }
}

double **SysDyn::getControlCoords() const
{
    return controlCoords;
}

unsigned long long int** SysDyn::getControlIntCoords()
{
    return controlIntCoords;
}

double* SysDyn::getLimSupC()
{
    return limSupC;
}
double* SysDyn::getLimInfC()
{
    return limInfC;
}
double* SysDyn::getStepC()
{
    return limInfC;
}
unsigned long long int SysDyn::getDimC() const
{
    return dimC;
}

double** SysDyn::getTychCoords() const
{
    return tychCoords;
}

unsigned long long int** SysDyn::getTychIntCoords()
{

    return tychIntCoords;
}

double* SysDyn::getLimSupTy()
{
    return limSupTy;
}
double* SysDyn::getLimInfTy()
{
    return limInfTy;
}
double* SysDyn::getStepTy()
{
    return stepTy;
}
unsigned long long int SysDyn::getDimTy() const
{
    return dimTy;
}
DynType SysDyn::getDynType()
{
    return dynType;
}

unsigned long long int* SysDyn::getNbPointsC()
{
    return nbPointsC;
}

unsigned long long int SysDyn::getTotalNbPointsC() const
{
    return totalNbPointsC;
}

unsigned long long int* SysDyn::getNbPointsTy()
{
    return nbPointsTy;
}

unsigned long long int SysDyn::getTotalNbPointsTy() const
{
    return totalNbPointsTych;
}

void SysDyn::FDiscretEuler(const double *x, const double *u, double *res, double rho) const
{

    int i;

    (*dynamics)(x, u, res);

    for (i = 0; i < dimS; i++)
    {
        res[i] = x[i] + rho * dynSignFactor * res[i];
    }

    grid->periodizePoint(res);

}

void SysDyn::FDiscretEuler_tych(const double *x, const double *u, const double *v, double *res, double rho) const
{

    int i;

    (*dynamics_tych)(x, u, v, res);

    for (i = 0; i < dimS; i++)
    {
        res[i] = x[i] + rho * dynSignFactor * res[i];
    }

    grid->periodizePoint(res);

}

void SysDyn::FDiscret(const double *x, const double *u, double *res, double rho) const
{

    (*dynamics)(x, u, res);
    grid->periodizePoint(res);
}

void SysDyn::FDiscret_tych(const double *x, const double *u, const double *v, double *res, double rho) const
{

    (*dynamics_tych)(x, u, v, res);
    grid->periodizePoint(res);
}


void SysDyn::FDiscretRK4(const double *x, const double *u, double *res, double rho) const
{
    int i;
    double *ki, *y;
    ki = new double[dimS];
    y = new double[dimS];
    (*dynamics)(x, u, ki);

    for (i = 0; i < dimS; i++)
    {
        res[i] = x[i] + rho * dynSignFactor * ki[i] / 6.0;
        y[i] = x[i] + 0.5 * rho * dynSignFactor * ki[i];
    }
    grid->periodizePoint(y);

    /*
     * k2=f(x+0.5*rho*k1,u)
     * res=res+rho*k2/3;
     *
     */
    (*dynamics)(y, u, ki);

    for (i = 0; i < dimS; i++)
    {
        y[i] = x[i] + 0.5 * rho * dynSignFactor * ki[i];
        res[i] = res[i] + rho * dynSignFactor * ki[i] / 3.0;
    }

    grid->periodizePoint(y);
    /*
     * k3=f(x+0.5*rho*k2,u)
     * res=res+rho*k3/3;
     *
     */
    (*dynamics)(y, u, ki);

    for (i = 0; i < dimS; i++)
    {
        y[i] = x[i] + rho * dynSignFactor * ki[i];
        res[i] = res[i] + rho * dynSignFactor * ki[i] / 3.0;
    }

    grid->periodizePoint(y);

    /*
     * k4=f(x+rho*k3,u)
     * res=res+rho*k4/6;
     *
     */
    (*dynamics)(y, u, ki);

    for (i = 0; i < dimS; i++)
    {
        res[i] = res[i] + rho * dynSignFactor * ki[i] / 6.0;
    }

    grid->periodizePoint(res);
    delete[] ki;
    delete[] y;
}

void SysDyn::FDiscretRK4_tych(const double *x, const double *u, const double *v, double *res, double rho) const
{
    int i;
    double *ki, *y;
    ki = new double[dimS];
    y = new double[dimS];
    (*dynamics_tych)(x, u, v, ki);

    for (i = 0; i < dimS; i++)
    {
        res[i] = x[i] + rho * dynSignFactor * ki[i] / 6.0;
        y[i] = x[i] + 0.5 * rho * dynSignFactor * ki[i];
    }
    grid->periodizePoint(y);

    /*
     * k2=f(x+0.5*rho*k1,u)
     * res=res+rho*k2/3;
     *
     */
    (*dynamics_tych)(y, u, v, ki);

    for (i = 0; i < dimS; i++)
    {
        y[i] = x[i] + 0.5 * rho * dynSignFactor * ki[i];
        res[i] = res[i] + rho * dynSignFactor * ki[i] / 3.0;
    }

    grid->periodizePoint(y);
    /*
     * k3=f(x+0.5*rho*k2,u)
     * res=res+rho*k3/3;
     *
     */
    (*dynamics_tych)(y, u, v, ki);

    for (i = 0; i < dimS; i++)
    {
        y[i] = x[i] + rho * dynSignFactor * ki[i];
        res[i] = res[i] + rho * dynSignFactor * ki[i] / 3.0;
    }

    grid->periodizePoint(y);

    /*
     * k4=f(x+rho*k3,u)
     * res=res+rho*k4/6;
     *
     */
    (*dynamics_tych)(y, u, v, ki);

    for (i = 0; i < dimS; i++)
    {
        res[i] = res[i] + rho * dynSignFactor * ki[i] / 6.0;
    }

    grid->periodizePoint(res);
    delete[] ki;
    delete[] y;
}


double SysDyn::calculRho_local(const double *x) const
{
    if (dynType == DC)
    {
        return 1.0;
    }
    double rho1;
    double h = grid->maxStep;
    double LL = ((this->*calcul_L))(x);
    double MFF = ((this->*calcul_M))(x);
    if (MFF * LL < 2.0 * h)
    {
        MFF = 1.0;
        LL = 1.0;
    }

    rho1 = sqrt((2.0 * h) / (LL * MFF));
    return rho1;
}

void SysDyn::FDiscretRK2(const double *x, const double *u, double *res, double rho) const
{

    int i;
    double *Fx, *Fres;
    Fx = new double[dimS];
    Fres = new double[dimS];
    (*dynamics)(x, u, Fx);

    //   calculRho_local(x  );

    for (i = 0; i < dimS; i++)
    {
        res[i] = x[i] + rho * dynSignFactor * Fx[i];
    }

    grid->periodizePoint(res);

    (*dynamics)(res, u, Fres);

    for (i = 0; i < dimS; i++)
    {
        res[i] = x[i] + 0.5 * rho * dynSignFactor * (Fx[i] + Fres[i]);
    }

    grid->periodizePoint(res);
    delete[] Fx;
    delete[] Fres;
}

void SysDyn::FDiscretRK2_tych(const double *x, const double *u, const double *v, double *res, double rho) const
{

    int i;
    double *Fx, *Fres;
    Fx = new double[dimS];
    Fres = new double[dimS];
    (*dynamics_tych)(x, u, v, Fx);

    //   calculRho_local(x  );

    for (i = 0; i < dimS; i++)
    {
        res[i] = x[i] + rho * dynSignFactor * Fx[i];
    }

    grid->periodizePoint(res);

    (*dynamics_tych)(res, u, v, Fres);

    for (i = 0; i < dimS; i++)
    {
        res[i] = x[i] + 0.5 * rho * dynSignFactor * (Fx[i] + Fres[i]);
    }

    grid->periodizePoint(res);
    delete[] Fx;
    delete[] Fres;
}

double SysDyn::calculL_local_num(const double *x) const
{
    double *xTempL = new double[dimS];
    double *FXmoinsHL = new double[dimS];
    double *FXplusHL = new double[dimS];

    double *infX = grid->limInf;
    double *pasX = grid->step;
    double *supX = grid->limSup;

    int i, j, k;
    for (i = 0; i < dimS; i++)
    {
        xTempL[i] = x[i];
    }
    double L1 = 0;

    bool test = false;

    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
    {
        for (j = 0; j < dimS; j++)
        {
            test = false;
            xTempL[j] = xTempL[j] - pasX[j];
            //on teste si l'indice courant de l'ensemble dilate n'est pas en dehors de l'espace
            if ((xTempL[j] <= supX[j]) && (xTempL[j] >= infX[j]))
            {
                // si l'indice est dans les limites de l'espace on  calcule le max de F(x,u) sur u

                (*dynamics)(xTempL, controlCoords[nu], FXmoinsHL);
                xTempL[j] = xTempL[j] + 2.0 * pasX[j];
            }
            else
            {
                xTempL[j] = x[j];
                (*dynamics)(xTempL, controlCoords[nu], FXmoinsHL);
                xTempL[j] = xTempL[j] + pasX[j];
                test = true;
            }

            if ((xTempL[j] <= supX[j]) && (xTempL[j] >= infX[j]))
            {
                // si l'indice est dans les limites de l'espace on  calcule le max de F(x,u) sur u

                (*dynamics)(xTempL, controlCoords[nu], FXplusHL);
                xTempL[j] = xTempL[j] - pasX[j];
            }
            else
            {
                xTempL[j] = xTempL[j] - pasX[j];
                (*dynamics)(xTempL, controlCoords[nu], FXplusHL);
                test = true;
            }

            for (k = 0; k < dimS; k++)
            {
                FXmoinsHL[k] = fabs(FXmoinsHL[k] - FXplusHL[k]);
                if (test)
                {
                    FXmoinsHL[k] /= pasX[k];
                }
                else
                {
                    FXmoinsHL[k] /= (2.0 * pasX[k]);
                }

                if (FXmoinsHL[k] > L1)
                {
                    L1 = FXmoinsHL[k];
                }
            }
        }
    }
    delete[] xTempL;
    delete[] FXplusHL;
    delete[] FXmoinsHL;
    return max(L1, lfunc_L);
}

double SysDyn::calculL_local_num_tych(const double *x) const
{
    double *xTempL = new double[dimS];
    double *FXmoinsHL = new double[dimS];
    double *FXplusHL = new double[dimS];

    double *infX = grid->limInf;
    double *pasX = grid->step;
    double *supX = grid->limSup;

    int i, j, k;
    for (i = 0; i < dimS; i++)
    {
        xTempL[i] = x[i];
    }
    double L1 = 0;

    bool test = false;

    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
    {
        for (unsigned long long int nv = 0; nv < totalNbPointsTych; nv++)
        {
            for (j = 0; j < dimS; j++)
            {
                test = false;
                xTempL[j] = xTempL[j] - pasX[j];
                //on teste si l'indice courant de l'ensemble dilate n'est pas en dehors de l'espace
                if ((xTempL[j] <= supX[j]) && (xTempL[j] >= infX[j]))
                {
                    // si l'indice est dans les limites de l'espace on  calcule le max de F(x,u) sur u

                    (*dynamics_tych)(xTempL, controlCoords[nu], tychCoords[nv], FXmoinsHL);
                    xTempL[j] = xTempL[j] + 2.0 * pasX[j];
                }
                else
                {
                    xTempL[j] = x[j];
                    (*dynamics_tych)(xTempL, controlCoords[nu], tychCoords[nv], FXmoinsHL);
                    xTempL[j] = xTempL[j] + pasX[j];
                    test = true;
                }

                if ((xTempL[j] <= supX[j]) && (xTempL[j] >= infX[j]))
                {
                    // si l'indice est dans les limites de l'espace on  calcule le max de F(x,u) sur u

                    (*dynamics_tych)(xTempL, controlCoords[nu], tychCoords[nv], FXplusHL);
                    xTempL[j] = xTempL[j] - pasX[j];
                }
                else
                {
                    xTempL[j] = xTempL[j] - pasX[j];
                    (*dynamics_tych)(xTempL, controlCoords[nu], tychCoords[nv], FXplusHL);
                    test = true;
                }

                for (k = 0; k < dimS; k++)
                {
                    FXmoinsHL[k] = fabs(FXmoinsHL[k] - FXplusHL[k]);
                    if (test)
                    {
                        FXmoinsHL[k] /= pasX[k];
                    }
                    else
                    {
                        FXmoinsHL[k] /= (2.0 * pasX[k]);
                    }

                    if (FXmoinsHL[k] > L1)
                    {
                        L1 = FXmoinsHL[k];
                    }
                }
            }
        }
    }
    delete[] xTempL;
    delete[] FXplusHL;
    delete[] FXmoinsHL;
    return max(L1, lfunc_L);
}

double SysDyn::calculL_local_ana(const double *x) const
{
    int j, k;
    double **jacob = new double*[dimS];

    for (int i = 0; i < dimS; i++)
    {
        jacob[i] = new double[dimS];
    }
    double L1 = 0;
    double norme;
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
    {

        (*jacobian)(x, controlCoords[nu], jacob);
        norme = 0.;
        for (k = 0; k < dimS; k++)
        {
            for (j = 0; j < dimS; j++)
            {
                norme = max(norme, abs(jacob[k][j]));
            }
        }
        L1 = max(L1, norme);
    }

    for (int i = 0; i < dimS; i++)
    {
        delete[] jacob[i];
    }
    delete[] jacob;
    return max(L1, lfunc_L);
}

double SysDyn::calculL_local_ana_tych(const double *x) const
{
    int j, k;
    double **jacob = new double*[dimS];

    for (int i = 0; i < dimS; i++)
    {
        jacob[i] = new double[dimS];
    }
    double L1 = 0;
    double norme;
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
    {
        for (unsigned long long int nv = 0; nv < totalNbPointsTych; nv++)
        {

            (*jacobian_tych)(x, controlCoords[nu], tychCoords[nv], jacob);
            norme = 0.;
            for (k = 0; k < dimS; k++)
            {
                for (j = 0; j < dimS; j++)
                {
                    norme = max(norme, abs(jacob[k][j]));
                }
            }
            L1 = max(L1, norme);
        }
    }

    for (int i = 0; i < dimS; i++)
    {
        delete[] jacob[i];
    }
    delete[] jacob;
    return max(L1, lfunc_L);
}

double SysDyn::returnL_local_ana(const double *x)  const
{
    return max(L, lfunc_L);
}
double SysDyn::returnMF_local_ana(const double *x) const
{
    return max(MF, lfunc_MF);
}

double SysDyn::calculMF_local_num(const double *x) const
{

    //calcul de la taille e prevoir pour les coordonnees des indices de debut de parcours
    //que la methode GPU va renvoyer

    double *image = new double[dimS];

    double MF1 = 0.0;
    double normeImage;
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
    {
        (*dynamics)(x, controlCoords[nu], image);
        normeImage = 0.0;
        for (int k = 0; k < dimS; k++)
        {
            normeImage = max(normeImage, abs(image[k]));
        }
        MF1 = max(MF1, normeImage);

    }
    delete[] image;
    return max(MF1, lfunc_MF);
}

double SysDyn::calculMF_local_num_tych(const double *x) const
{

    //calcul de la taille e prevoir pour les coordonnees des indices de debut de parcours
    //que la methode GPU va renvoyer

    double *image = new double[dimS];

    double MF1 = 0.0;
    double normeImage;
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
    {
        for (unsigned long long int nv = 0; nv < totalNbPointsTych; nv++)
        {
            (*dynamics_tych)(x, controlCoords[nu], tychCoords[nv], image);
            normeImage = 0.0;
            for (int k = 0; k < dimS; k++)
            {
                normeImage = max(normeImage, abs(image[k]));
            }
            MF1 = max(MF1, normeImage);
        }
    }
    delete[] image;
    return max(MF1, lfunc_MF);
}


double SysDyn::calculMF_local_ana(const double *x) const
{

    double *image = new double[dimS];
    double normeImage;
    (*localDynBounds)(x, image);

    normeImage = 0.0;
    for (int k = 0; k < dimS; k++)
    {
        normeImage = max(normeImage, abs(image[k]));
    }
    double MF1 = normeImage;
    delete[] image;
    return max(MF1, lfunc_MF);
}

SysDyn::~SysDyn()
{
    unsigned long long int k;
    for (k = 0; k < totalNbPointsC; k++)
    {
        delete[] controlCoords[k];
    }
    delete[] controlCoords;
    delete[] image;
    delete[] xTemp;
    delete[] FXmoinsH;
    delete[] FXplusH;
}

string SysDyn::getRetroFileName()
{
    return retroFileName;
}

int SysDyn::getFDDynType()
{
    return fd_dyn_type;
}
bool SysDyn::isTimeStepGlobal()
{
    return globalTimeStep;
}
void SysDyn::setDynamicsForward()
{
    dynSignFactor = 1.0;
}
void SysDyn::setDynamicsBackward()
{
    dynSignFactor = -1.0;
}

SysDyn::PointStatus SysDyn::checkKernelRelation(double *position) const
{
    int cptOK = 0;
    int cellNum;
    int posTemp;
    
    const long long int *indicesDecalCell = grid->getIndicesDecalCell();
    int pow2 = (int) pow(2.0, dimS);

    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */
    long long unsigned *testI = new long long unsigned[dimS];

    PointStatus res = VALID_TRAJECTORY_POINT;
    
    if (!grid->isPointInGrid(position))
    {
        res = OUTSIDE_GRID;
    }
    else if (constraintsX(position) >= PLUS_INF)
    {
        res = OUTSIDE_CONSTRAINTS;
    }
    else {
        cellNum = grid->localizePoint(position);

        int ii = 0;
        bool inSet = false;
        while (ii < pow2 && !inSet) {
            posTemp = cellNum + indicesDecalCell[ii];
            grid->numToIntCoords(posTemp, testI);
            inSet = (grid->isInSet(testI));
            ++ii;
        }        
        if (!inSet) {
            res = OUTSIDE_DOMAIN;
        }
    }
    
    delete[] testI;

    return res;
}

bool SysDyn::isViableControl(const double *currentPos, const double *controlCoord, double *imageVect, double rho) const {
    if (constraintsXU(currentPos, controlCoord) >= PLUS_INF) {
        return false;
    }
    else {
        std::invoke(discretDynamics, this, currentPos, controlCoord, imageVect, rho);
        return checkKernelRelation(imageVect) == VALID_TRAJECTORY_POINT;
    }
}

bool SysDyn::isViableControl_tych(const double *currentPos, const double *controlCoord, const double *tycheCoord, double *imageVect, double rho) const {
    if (constraintsXU(currentPos, controlCoord) >= PLUS_INF || constraintsXV_tych(currentPos, tycheCoord) >= PLUS_INF) {
        return false;
    }
    else {
        std::invoke(discretDynamics_tych, this, currentPos, controlCoord, tycheCoord, imageVect, rho);
        return checkKernelRelation(imageVect) == VALID_TRAJECTORY_POINT;
    }
}

bool SysDyn::isViableGuaranteedControl(const double *currentPos, const double *controlCoord, double rho) const {

    if (constraintsXU(currentPos, controlCoord) >= PLUS_INF) {
        return false;
    }
    
    bool allTychOk = true;
    unsigned long long int i = 0;
    double *imageVect = new double[getDim()];
    while (allTychOk && i < totalNbPointsTych) {
        const double *tychCoord = tychCoords[i];
        if (constraintsXV_tych(currentPos, tychCoord) < PLUS_INF) {
            std::invoke(discretDynamics_tych, this, currentPos, controlCoord, tychCoord, imageVect, rho);
            allTychOk = (checkKernelRelation(imageVect) == VALID_TRAJECTORY_POINT);
        }
        ++i;
    }
    delete [] imageVect;
    return allTychOk;
}

void SysDyn::getTychasticImage(const double *x, const double *u, const double *v, double *imageVect, double rho) const {
    std::invoke(discretDynamics_tych, this, x, u, v, imageVect, rho);
}

const Grid *SysDyn::getGrid() const {
    return grid;
}

int SysDyn::getDim() const {
    return dimS;
}
