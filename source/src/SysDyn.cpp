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
#include "../include/SimpleSysDyn.h"
#include "../include/TychasticSysDyn.h"
#include "../include/HybridSysDyn.h"
#include <stdexcept>

SysDyn::SysDyn()
    {

    }

SysDyn::SysDyn(const systemParams &SP, int ds, const controlParams &cp, Grid *grRef)
    {
    isTychastic = cp.DIM_TY > 0;
    isHybrid = false;

    spdlog::info("[System] : Looking for control parameter : control dim is {}", cp.DIMC);
    if (cp.DIMC > 0)
	{
	spdlog::info("[System] : Starting build of ControlGrid");
	controls = std::make_unique<ControlGrid>(cp.DIMC, cp.LIMINFC, cp.LIMSUPC, cp.NBPOINTSC);
	}
    else
	{
	controls = std::make_unique<ControlGrid>();
	}

    if (cp.DIM_TY > 0)
	{
	tyches = std::make_unique<ControlGrid>(cp.DIM_TY, cp.LIMINF_TY, cp.LIMSUP_TY, cp.NBPOINTS_TY);
	}
    else
	{
	tyches = std::make_unique<ControlGrid>();
	}

    if (cp.DIM_HT > 0)
	{
	hybridTransistionControls = std::make_unique<ControlGrid>(cp.DIM_HT, cp.LIMINF_HT, cp.LIMSUP_HT, cp.NBPOINTS_HT);
	}
    else
	{
	hybridTransistionControls = std::make_unique<ControlGrid>();
	}

    localDynBounds = SP.LOCAL_DYN_BOUNDS;


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

    initializeSubSystems(SP, ds, 0, cp, grRef);
    initializeMethods(SP);

    spdlog::info("[System] : Dynamic system initialization finished");
    }

SysDyn::SysDyn(const systemParams &SP, int continuousStateDim, int discretStateDim, const controlParams &cp, Grid *grRef)
    {

    isTychastic = (cp.DIM_TY > 0);
    isHybrid = cp.DIM_HT > 0 || discretStateDim > 0;
    if (cp.DIMC > 0)
	{
	spdlog::info("[System] : Starting build of ControlGrid");
	controls = std::make_unique<ControlGrid>(cp.DIMC, cp.LIMINFC, cp.LIMSUPC, cp.NBPOINTSC);
	}
    else
	{
	controls = std::make_unique<ControlGrid>();
	}

    if (cp.DIM_TY > 0)
	{
	tyches = std::make_unique<ControlGrid>(cp.DIM_TY, cp.LIMINF_TY, cp.LIMSUP_TY, cp.NBPOINTS_TY);
	}
    else
	{
	tyches = std::make_unique<ControlGrid>();
	}

    if (cp.DIM_HT > 0)
	{
	hybridTransistionControls = std::make_unique<ControlGrid>(cp.DIM_HT, cp.LIMINF_HT, cp.LIMSUP_HT, cp.NBPOINTS_HT);
	}
    else
	{
	hybridTransistionControls = std::make_unique<ControlGrid>();
	}

    localDynBounds = SP.LOCAL_DYN_BOUNDS;


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

    initializeSubSystems(SP, continuousStateDim, discretStateDim, cp, grRef);
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
    if (!simpleSystem)
	{
	throw std::runtime_error("[SysDyn] simpleSystem not initialized");
	}
    simpleSystem->FDiscretEuler(x, u, res, rho);
    }

void SysDyn::FDiscret(const double *x, const double *u, double *res, double rho) const
    {
    if (!simpleSystem)
	{
	throw std::runtime_error("[SysDyn] simpleSystem not initialized");
	}
    simpleSystem->FDiscret(x, u, res, rho);
    }

void SysDyn::FDiscretRK4(const double *x, const double *u, double *res, double rho) const
    {
    if (!simpleSystem)
	{
	throw std::runtime_error("[SysDyn] simpleSystem not initialized");
	}
    simpleSystem->FDiscretRK4(x, u, res, rho);
    }

void SysDyn::FDiscretRK2(const double *x, const double *u, double *res, double rho) const
    {
    if (!simpleSystem)
	{
	throw std::runtime_error("[SysDyn] simpleSystem not initialized");
	}
    simpleSystem->FDiscretRK2(x, u, res, rho);
    }


void SysDyn::FDiscretRK4_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int * resd, double *tempResC, unsigned long long int *tempResD, double rho) const
    {
    if (!hybridSystem)
	{
	throw std::runtime_error("[SysDyn] hybridSystem not initialized");
	}
    hybridSystem->FDiscretRK4_hybrid(xc, xd, uc, ud, resc, resd, tempResC, tempResD, rho);
    }

void SysDyn::FDiscretRK2_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int * resd, double *tempResC, unsigned long long int *tempResD, double rho) const
    {
    if (!hybridSystem)
	{
	throw std::runtime_error("[SysDyn] hybridSystem not initialized");
	}
    hybridSystem->FDiscretRK2_hybrid(xc, xd, uc, ud, resc, resd, tempResC, tempResD, rho);
    }

void SysDyn::FDiscretEuler_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int * resd, double *tempResC, unsigned long long int *tempResD,  double rho) const
    {
    if (!hybridSystem)
	{
	throw std::runtime_error("[SysDyn] hybridSystem not initialized");
	}
    hybridSystem->FDiscretEuler_hybrid(xc, xd, uc, ud, resc, resd, tempResC, tempResD, rho);
    }

void SysDyn::FDiscret_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int * resd, double *tempResC, unsigned long long int *tempResD, double rho) const
    {
    if (!hybridSystem)
	{
	throw std::runtime_error("[SysDyn] hybridSystem not initialized");
	}
    hybridSystem->FDiscret_hybrid(xc, xd, uc, ud, resc, resd, tempResC, tempResD, rho);
    }

void SysDyn::FDiscret_tych(const double *x, const double *u, const double *v, double *res, double rho) const
    {
    if (!tychasticSystem)
	{
	throw std::runtime_error("[SysDyn] tychasticSystem not initialized");
	}
    tychasticSystem->FDiscret_tych(x, u, v, res, rho);
    }

void SysDyn::FDiscretEuler_tych(const double *x, const double *u, const double *v, double *res, double rho) const
    {
    if (!tychasticSystem)
	{
	throw std::runtime_error("[SysDyn] tychasticSystem not initialized");
	}
    tychasticSystem->FDiscretEuler_tych(x, u, v, res, rho);
    }

void SysDyn::FDiscretRK2_tych(const double *x, const double *u, const double *v, double *res, double rho) const
    {
    if (!tychasticSystem)
	{
	throw std::runtime_error("[SysDyn] tychasticSystem not initialized");
	}
    tychasticSystem->FDiscretRK2_tych(x, u, v, res, rho);
    }

void SysDyn::FDiscretRK4_tych(const double *x, const double *u, const double *v, double *res, double rho) const
    {
    if (!tychasticSystem)
	{
	throw std::runtime_error("[SysDyn] tychasticSystem not initialized");
	}
    tychasticSystem->FDiscretRK4_tych(x, u, v, res, rho);
    }

double SysDyn::calculRho_local(const double *x) const
    {
    if (!simpleSystem)
	{
	throw std::runtime_error("[SysDyn] simpleSystem not initialized");
	}
    return simpleSystem->calculRho_local(x);
    }

double SysDyn::calculRho_local_hybrid(const double *xc, const unsigned long long int * xd) const
    {
    if (!hybridSystem)
	{
	throw std::runtime_error("[SysDyn] hybridSystem not initialized");
	}
    return hybridSystem->calculRho_local_hybrid(xc, xd);
    }

double SysDyn::calculL_local_num(const double *x) const
    {
    if (!simpleSystem)
	{
	throw std::runtime_error("[SysDyn] simpleSystem not initialized");
	}
    return simpleSystem->calculL_local_num(x);
    }

double SysDyn::calculL_local_num_hybrid(const double *xc, const unsigned long long int * xd) const
    {
    if (!hybridSystem)
	{
    	throw std::runtime_error("[SysDyn] hybridSystem not initialized");
	}
	return hybridSystem->calculL_local_num_hybrid(xc, xd);
    }

double SysDyn::calculL_local_num_tych(const double *x) const
    {
    if (!tychasticSystem)
	{
	throw std::runtime_error("[SysDyn] tychasticSystem not initialized");
	}
    return tychasticSystem->calculL_local_num_tych(x);
    }

double SysDyn::calculL_local_ana(const double *x) const
    {
    if (!simpleSystem)
	{
	throw std::runtime_error("[SysDyn] simpleSystem not initialized");
	}
    return simpleSystem->calculL_local_ana(x);
    }

double SysDyn::calculL_local_ana_hybrid(const double *xc, const unsigned long long int * xd) const
    {
    if (!hybridSystem)
	{
	throw std::runtime_error("[SysDyn] hybridSystem not initialized");
	}
    return hybridSystem->calculL_local_ana_hybrid(xc, xd);
    }

double SysDyn::calculL_local_ana_tych(const double *x) const
    {
    if (!tychasticSystem)
	{
	throw std::runtime_error("[SysDyn] tychasticSystem not initialized");
	}
    return tychasticSystem->calculL_local_ana_tych(x);
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
    if (!simpleSystem)
	{
	throw std::runtime_error("[SysDyn] simpleSystem not initialized");
	}
    return simpleSystem->calculMF_local_num(x);
    }

double SysDyn::calculMF_local_num_hybrid(const double *xc, const unsigned long long int * xd) const
    {
	if (!hybridSystem)
	{
		throw std::runtime_error("[SysDyn] hybridSystem not initialized");
	}
	return hybridSystem->calculMF_local_num_hybrid(xc, xd);
    }

double SysDyn::calculMF_local_num_tych(const double *x) const
    {
	if (!tychasticSystem)
	{
		throw std::runtime_error("[SysDyn] tychasticSystem not initialized");
	}
	return tychasticSystem->calculMF_local_num_tych(x);
    }

double SysDyn::calculMF_local_ana(const double *x) const
    {
	if (!simpleSystem)
	{
		throw std::runtime_error("[SysDyn] simpleSystem not initialized");
	}
	return simpleSystem->calculMF_local_ana(x);
    }

double SysDyn::calculMF_local_ana_hybrid(const double *xc, const unsigned long long int * xd) const
    {
	if (!hybridSystem)
	{
		throw std::runtime_error("[SysDyn] hybridSystem not initialized");
	}
	return hybridSystem->calculMF_local_ana_hybrid(xc, xd);
    }

SysDyn::~SysDyn()
    {
    simpleSystem.reset();
    tychasticSystem.reset();
    hybridSystem.reset();
    hybridTransistionControls.reset();
    tyches.reset();
    controls.reset();
    }

int SysDyn::getFDDynType()
    {
    return fd_dyn_type;
    }
bool SysDyn::isTimeStepGlobal()
    {
    return simpleSystem->isTimeStepGlobal();
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
    return simpleSystem->getGrid();
    }

int SysDyn::getDim() const
    {
    return simpleSystem->getDim();
    }

void SysDyn::initializeSubSystems(const systemParams &SP, int continuousStateDim, int discreteStateDim, const controlParams &cp, Grid *refGrid)
    {
    int totalDim = continuousStateDim + discreteStateDim;
    simpleSystem = std::make_unique<SimpleSysDyn>(SP, totalDim, cp, refGrid, controls.get());
    tychasticSystem = std::make_unique<TychasticSysDyn>(SP, totalDim, cp, refGrid, controls.get(), tyches.get());
    hybridSystem = std::make_unique<HybridSysDyn>(SP, continuousStateDim, discreteStateDim, cp, refGrid, controls.get(), hybridTransistionControls.get());
    }

