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
    isTychastic = cp.DIM_TY > 0;
    isHybrid = false;

    spdlog::info("[System] : Looking for control parameter : control dim is {}", cp.DIMC);
    if (cp.DIMC > 0)
	{
	spdlog::info("[System] : Starting build of ControlGrid");
	controls = new ControlGrid(cp.DIMC, cp.LIMINFC, cp.LIMSUPC, cp.NBPOINTSC);
	}
    else
	{
	controls = new ControlGrid();
	}

    if (cp.DIM_TY > 0)
	{
	tyches = new ControlGrid(cp.DIM_TY, cp.LIMINF_TY, cp.LIMSUP_TY, cp.NBPOINTS_TY);
	}
    else
	{
	tyches = new ControlGrid();
	}

    if (cp.DIM_HT > 0)
	{
	hybridTransistionControls = new ControlGrid(cp.DIM_HT, cp.LIMINF_HT, cp.LIMSUP_HT, cp.NBPOINTS_HT);
	}
    else
	{
	hybridTransistionControls = new ControlGrid();
	}
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

    initializeMethods(SP);

    spdlog::info("[System] : Dynamic system initialization finished");
    }

SysDyn::SysDyn(const systemParams &SP, int continuousStateDim, int discretStateDim, const controlParams &cp, Grid *grRef)
    {
    grid = grRef;
    dimS = continuousStateDim + discretStateDim;
    dimS_hc = continuousStateDim;
    dimS_hd = discretStateDim;

    isTychastic = (cp.DIM_TY > 0);
    isHybrid = cp.DIM_HT > 0 || dimS_hd > 0;
    if (cp.DIMC > 0)
	{
	spdlog::info("[System] : Starting build of ControlGrid");
	controls = new ControlGrid(cp.DIMC, cp.LIMINFC, cp.LIMSUPC, cp.NBPOINTSC);
	}
    else
	{
	controls = new ControlGrid();
	}

    if (cp.DIM_TY > 0)
	{
	tyches = new ControlGrid(cp.DIM_TY, cp.LIMINF_TY, cp.LIMSUP_TY, cp.NBPOINTS_TY);
	}
    else
	{
	tyches = new ControlGrid();
	}

    if (cp.DIM_HT > 0)
	{
	hybridTransistionControls = new ControlGrid(cp.DIM_HT, cp.LIMINF_HT, cp.LIMSUP_HT, cp.NBPOINTS_HT);
	}
    else
	{
	hybridTransistionControls = new ControlGrid();
	}

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

    initializeMethods(SP);

    spdlog::info("[System] : Dynamic system initialization finished");
    }

void SysDyn::initializeMethods(const systemParams &SP)
    {
    dynamics = SP.DYNAMICS;
    dynamics_fd = SP.DYNAMICS_FD;
    dynamics_tych_fd = SP.DYNAMICS_TYCH_FD;
    dynamics_tych = SP.DYNAMICS_TYCH;
    dynamics_hybrid_c = SP.DYNAMICS_HYBRID_C;
    dynamics_hybrid_d = SP.DYNAMICS_HYBRID_D;
    resetmap_hybrid = SP.RESET_MAP_HYBRID;

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
    constraintsXU_hybrid = SP.CONSTR_XU_HYBRID;
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

    switch (computeMF)
	{
    case ANALYTICAL:
	SysDyn::calcul_M = &SysDyn::returnMF_local_ana;
	SysDyn::calcul_M_hybrid = &SysDyn::returnMF_local_ana_hybrid;
	break;
    case ANALYTICAL_CALC:

	    SysDyn::calcul_M_hybrid = &SysDyn::calculMF_local_ana_hybrid;
	    SysDyn::calcul_M = &SysDyn::calculMF_local_ana;
	break;
    case NUMERICAL_CALC:
	    SysDyn::calcul_M_hybrid = &SysDyn::calculMF_local_num_hybrid;
	    SysDyn::calcul_M = isTychastic ? &SysDyn::calculMF_local_num_tych : &SysDyn::calculMF_local_num;

	break;
	}

    switch (computeLC)
	{
    case ANALYTICAL:
	SysDyn::calcul_L = &SysDyn::returnL_local_ana;
	SysDyn::calcul_L_hybrid = &SysDyn::returnL_local_ana_hybrid;
	break;
    case ANALYTICAL_CALC:
	    SysDyn::calcul_L_hybrid = &SysDyn::calculL_local_ana_hybrid;
	    SysDyn::calcul_L = isTychastic ? &SysDyn::calculL_local_ana_tych : &SysDyn::calculL_local_ana;
	break;
    case NUMERICAL_CALC:
	    SysDyn::calcul_L_hybrid = &SysDyn::calculL_local_num_hybrid;
	    SysDyn::calcul_L = isTychastic ? &SysDyn::calculL_local_num_tych : &SysDyn::calculL_local_num;
	break;
	}

    discretisation = SP.SCHEME;

    if (dynType == DC || dynType == DH)
	{
	discretisation = NO_DISCRETIZATION_SCHEME;
	}

    switch (discretisation)
	{
    case NO_DISCRETIZATION_SCHEME:
	SysDyn::discretDynamics = &SysDyn::FDiscret;
	SysDyn::discretDynamics_tych = &SysDyn::FDiscret_tych;
	SysDyn::discretDynamics_hybrid = &SysDyn::FDiscret_hybrid;
	break;
    case EL:
	SysDyn::discretDynamics = &SysDyn::FDiscretEuler;
	SysDyn::discretDynamics_tych = &SysDyn::FDiscretEuler_tych;
	SysDyn::discretDynamics_hybrid = &SysDyn::FDiscretEuler_hybrid;
	break;
    case RK2:
	SysDyn::discretDynamics = &SysDyn::FDiscretRK2;
	SysDyn::discretDynamics_tych = &SysDyn::FDiscretRK2_tych;
	SysDyn::discretDynamics_hybrid = &SysDyn::FDiscretRK2_hybrid;
	break;
    case RK4:
	SysDyn::discretDynamics = &SysDyn::FDiscretRK4;
	SysDyn::discretDynamics_tych = &SysDyn::FDiscretRK4_tych;
	SysDyn::discretDynamics_hybrid = &SysDyn::FDiscretRK4_hybrid;
	break;
	}
    }

double** SysDyn::getControlCoords() const
    {
    return controls->GetControlCoords();
    }

unsigned long long int** SysDyn::getControlIntCoords()
    {
    return controls->GetControlIntCoords();
    }

double* SysDyn::getLimSupC()
    {
    return controls->GetLimSup();
    }
double* SysDyn::getLimInfC()
    {
    return controls->GetLimInf();
    }
double* SysDyn::getStepC()
    {
    return controls->GetStep();
    }
unsigned long long int SysDyn::getDimC() const
    {
    return controls->GetDim();
    }

double** SysDyn::getTychCoords() const
    {
    return tyches->GetControlCoords();
    }

unsigned long long int** SysDyn::getTychIntCoords()
    {

    return tyches->GetControlIntCoords();
    }

double* SysDyn::getLimSupTy()
    {
    return tyches->GetLimSup();
    }
double* SysDyn::getLimInfTy()
    {
    return tyches->GetLimInf();
    }
double* SysDyn::getStepTy()
    {
    return tyches->GetStep();
    }
unsigned long long int SysDyn::getDimTy() const
    {
    return tyches->GetDim();
    }
DynType SysDyn::getDynType()
    {
    return dynType;
    }

unsigned long long int* SysDyn::getNbPointsC()
    {
    return controls->GetNbPoints();
    }

unsigned long long int SysDyn::getTotalNbPointsC() const
    {
    return controls->GetTotalNbPoints();
    }

unsigned long long int* SysDyn::getNbPointsTy()
    {
    return tyches->GetNbPoints();
    }

unsigned long long int SysDyn::getTotalNbPointsTy() const
    {
    return tyches->GetTotalNbPoints();
    }

unsigned long long int* SysDyn::getNbPointsHybrid()
    {
    return hybridTransistionControls->GetNbPoints();
    }

unsigned long long int SysDyn::getTotalNbPointsHybrid() const
    {
    return hybridTransistionControls->GetTotalNbPoints();
    }

double** SysDyn::getHybridCoords() const
    {
    return hybridTransistionControls->GetControlCoords();
    }

unsigned long long int** SysDyn::getHybridIntCoords()
    {

    return hybridTransistionControls->GetControlIntCoords();
    }

double* SysDyn::getLimSupHybrid()
    {
    return hybridTransistionControls->GetLimSup();
    }
double* SysDyn::getLimInfHybrid()
    {
    return hybridTransistionControls->GetLimInf();
    }
double* SysDyn::getStepHybrid()
    {
    return hybridTransistionControls->GetStep();
    }
unsigned long long int SysDyn::getDimHybrid() const
    {
    return hybridTransistionControls->GetDim();
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

void SysDyn::FDiscretRK4_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int * resd, double *tempResC, unsigned long long int *tempResD, double rho) const
    {
    int i;
    double *ki, *y;
    ki = new double[dimS_hc];
    y = new double[dimS_hc];

    FDiscret_hybrid(xc, xd, uc, ud, tempResC, tempResD, ki, resd, rho);

    for (i = 0; i < dimS_hc; i++)
	{
	resc[i] = tempResC[i] + rho * dynSignFactor * ki[i] / 6.0;
	y[i] = tempResC[i] + 0.5 * rho * dynSignFactor * ki[i];
	}
    grid->periodizePoint(y);

    /*
     * k2=f(x+0.5*rho*k1,u)
     * res=res+rho*k2/3;
     *
     */

    (*dynamics_hybrid_c)(y, tempResD, uc, ki);

    for (i = 0; i < dimS_hc; i++)
	{
	y[i] = tempResC[i] + 0.5 * rho * dynSignFactor * ki[i];
	resc[i] = resc[i] + rho * dynSignFactor * ki[i] / 3.0;
	}

    grid->periodizePoint(y);
    /*
     * k3=f(x+0.5*rho*k2,u)
     * res=res+rho*k3/3;
     *
     */

    (*dynamics_hybrid_c)(y, tempResD, uc, ki);

    for (i = 0; i < dimS_hc; i++)
	{
	y[i] = tempResC[i] + rho * dynSignFactor * ki[i];
	resc[i] = resc[i] + rho * dynSignFactor * ki[i] / 3.0;
	}

    grid->periodizePoint(y);

    /*
     * k4=f(x+rho*k3,u)
     * res=res+rho*k4/6;
     *
     */

    (*dynamics_hybrid_c)(y, tempResD, uc, ki);
    for (i = 0; i < dimS_hc; i++)
	{
	resc[i] = resc[i] + rho * dynSignFactor * ki[i] / 6.0;
	}

    grid->periodizePoint(resc);
    delete[] ki;
    delete[] y;
    }


void SysDyn::FDiscretRK2_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int * resd, double *tempResC, unsigned long long int *tempResD, double rho) const
    {

    int i;
    double *Fx, *Fres;
    Fx = new double[dimS_hc];
    Fres = new double[dimS_hc];
    FDiscret_hybrid(xc, xd, uc, ud, tempResC, tempResD, Fx, resd, rho);

    //   calculRho_local(x  );

    for (i = 0; i < dimS_hc; i++)
	{
	resc[i] = tempResC[i] + rho * dynSignFactor * Fx[i];
	}

    grid->periodizePoint(resc);

    (*dynamics_hybrid_c)(resc, tempResD, uc, Fres);

    for (i = 0; i < dimS_hc; i++)
	{
	resc[i] = tempResC[i] + 0.5 * rho * dynSignFactor * (Fx[i] + Fres[i]);
	}

    grid->periodizePoint(resc);
    delete[] Fx;
    delete[] Fres;
    }


void SysDyn::FDiscretEuler_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int * resd, double *tempResC, unsigned long long int *tempResD,  double rho) const
    {
    FDiscret_hybrid(xc, xd, uc, ud, tempResC, tempResD, resc, resd, rho);
    int i;


    for (i = 0; i < dimS_hc; i++)
	{
	resc[i] = tempResC[i] + rho * dynSignFactor * resc[i];
	}

    grid->periodizePoint(resc);

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

void SysDyn::FDiscret_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int * resd, double *tempResC, unsigned long long int *tempResD, double rho) const
    {
    //Step 1 : impulse transistion
    (*resetmap_hybrid)(xc, xd, ud, tempResC, tempResD);
    //Step 2  dynamics

    // 2.A : dynamics of discrete state
    (*dynamics_hybrid_d)(tempResC, tempResD, ud, resd);

    //2.B Continuous evolution
    (*dynamics_hybrid_c)(tempResC, tempResD, uc, resc);
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

double SysDyn::calculRho_local_hybrid(const double *xc, const unsigned long long int * xd) const
    {
    if (dynType == DC)
	{
	return 1.0;
	}
    double rho1;
    double h = grid->maxStep;
    double LL = ((this->*calcul_L_hybrid))(xc, xd);
    double MFF = ((this->*calcul_M_hybrid))(xc, xd);
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
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();
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

double SysDyn::calculL_local_num_hybrid(const double *xc, const unsigned long long int * xd) const
    {
    double *xTempL = new double[dimS_hc];
    double *FXmoinsHL = new double[dimS_hc];
    double *FXplusHL = new double[dimS_hc];

    double *infX = grid->limInf;
    double *pasX = grid->step;
    double *supX = grid->limSup;

    int i, j, k;
    for (i = 0; i < dimS_hc; i++)
	{
	xTempL[i] = xc[i];
	}
    double L1 = 0;

    bool test = false;
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
	{
	for (j = 0; j < dimS_hc; j++)
	    {
	    test = false;
	    xTempL[j] = xTempL[j] - pasX[j];
	    //on teste si l'indice courant de l'ensemble dilate n'est pas en dehors de l'espace
	    if ((xTempL[j] <= supX[j]) && (xTempL[j] >= infX[j]))
		{
		// si l'indice est dans les limites de l'espace on  calcule le max de F(x,u) sur u

		(*dynamics_hybrid_c)(xTempL, xd, controlCoords[nu], FXmoinsHL);
		xTempL[j] = xTempL[j] + 2.0 * pasX[j];
		}
	    else
		{
		xTempL[j] = xc[j];
		(*dynamics_hybrid_c)(xTempL, xd, controlCoords[nu], FXmoinsHL);
		xTempL[j] = xTempL[j] + pasX[j];
		test = true;
		}

	    if ((xTempL[j] <= supX[j]) && (xTempL[j] >= infX[j]))
		{
		// si l'indice est dans les limites de l'espace on  calcule le max de F(x,u) sur u

		(*dynamics_hybrid_c)(xTempL, xd, controlCoords[nu], FXplusHL);
		xTempL[j] = xTempL[j] - pasX[j];
		}
	    else
		{
		xTempL[j] = xTempL[j] - pasX[j];
		(*dynamics_hybrid_c)(xTempL, xd, controlCoords[nu], FXplusHL);
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
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    unsigned long long int totalNbPointsTych = tyches->GetTotalNbPoints();

    double **controlCoords = controls->GetControlCoords();
    double **tychCoords = tyches->GetControlCoords();
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
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();
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

double SysDyn::calculL_local_ana_hybrid(const double *xc, const unsigned long long int * xd) const
    {
    int j, k;
    double **jacob = new double*[dimS_hc];

    for (int i = 0; i < dimS_hc; i++)
	{
	jacob[i] = new double[dimS_hc];
	}
    double L1 = 0;
    double norme;
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
	{

	(*jacobian_hybrid)(xc, xd, controlCoords[nu], jacob);
	norme = 0.;
	for (k = 0; k < dimS_hc; k++)
	    {
	    for (j = 0; j < dimS_hc; j++)
		{
		norme = max(norme, abs(jacob[k][j]));
		}
	    }
	L1 = max(L1, norme);
	}

    for (int i = 0; i < dimS_hc; i++)
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
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    unsigned long long int totalNbPointsTych = tyches->GetTotalNbPoints();

    double **controlCoords = controls->GetControlCoords();
    double **tychCoords = tyches->GetControlCoords();
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

double SysDyn::returnL_local_ana(const double *x) const
    {
    return max(L, lfunc_L);
    }
double SysDyn::returnMF_local_ana(const double *x) const
    {
    return max(MF, lfunc_MF);
    }

double SysDyn::returnL_local_ana_hybrid(const double *xc, const unsigned long long int * xd) const
    {
    return max(L, lfunc_L);
    }
double SysDyn::returnMF_local_ana_hybrid(const double *xc, const unsigned long long int * xd) const
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
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();

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

double SysDyn::calculMF_local_num_hybrid(const double *xc, const unsigned long long int * xd) const
    {

    //calcul de la taille e prevoir pour les coordonnees des indices de debut de parcours
    //que la methode GPU va renvoyer

    double *image = new double[dimS_hc];

    double MF1 = 0.0;
    double normeImage;
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();

    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
	{
	(*dynamics_hybrid_c)(xc, xd, controlCoords[nu], image);
	normeImage = 0.0;
	for (int k = 0; k < dimS_hc; k++)
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
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    unsigned long long int totalNbPointsTych = tyches->GetTotalNbPoints();

    double **controlCoords = controls->GetControlCoords();
    double **tychCoords = tyches->GetControlCoords();
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

double SysDyn::calculMF_local_ana_hybrid(const double *xc, const unsigned long long int * xd) const
    {

    double *image = new double[dimS_hc];
    double normeImage;
    (*localDynBounds_hybrid)(xc, xd, image);

    normeImage = 0.0;
    for (int k = 0; k < dimS_hc; k++)
	{
	normeImage = max(normeImage, abs(image[k]));
	}
    double MF1 = normeImage;
    delete[] image;
    return max(MF1, lfunc_MF);
    }

SysDyn::~SysDyn()
    {
    delete controls;
    delete tyches;
    delete hybridTransistionControls;
    delete[] image;
    delete[] xTemp;
    delete[] FXmoinsH;
    delete[] FXplusH;
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

void SysDyn::getTychasticImage(const double *x, const double *u, const double *v, double *imageVect, double rho) const
    {
    std::invoke(discretDynamics_tych, this, x, u, v, imageVect, rho);
    }

const Grid* SysDyn::getGrid() const
    {
    return grid;
    }

int SysDyn::getDim() const
    {
    return dimS;
    }
