#include "../include/HybridSysDyn.h"

#include <algorithm>
#include <cmath>

#include <spdlog/spdlog.h>

using std::abs;
using std::max;

HybridSysDyn::HybridSysDyn()
    : SimpleSysDyn(), dimS_hc(0), dimS_hd(0), hybridTransistionControls(new ControlGrid()), dynamics_hybrid_d(nullptr), dynamics_hybrid_c(nullptr),
      resetmap_hybrid(nullptr), constraintsXU_hybrid(nullptr), constraintsX_hybrid(nullptr), localDynBounds_hybrid(nullptr), jacobian_hybrid(nullptr),
      calculLHybridFunc(), calculMHybridFunc(), discretDynamicsHybrid(), ownsHybridControls(true)
    {
    isHybrid = true;
    }

HybridSysDyn::HybridSysDyn(const systemParams &SP, int continuousStateDim, int discreteStateDim, const controlParams &cp, Grid *refGrid, ControlGrid *controlGrid, ControlGrid *hybridControls)
    : SimpleSysDyn(SP, continuousStateDim + discreteStateDim, cp, refGrid, controlGrid), dimS_hc(continuousStateDim), dimS_hd(discreteStateDim),
      hybridTransistionControls(hybridControls), dynamics_hybrid_d(nullptr), dynamics_hybrid_c(nullptr), resetmap_hybrid(nullptr),
      constraintsXU_hybrid(nullptr), constraintsX_hybrid(nullptr), localDynBounds_hybrid(nullptr), jacobian_hybrid(nullptr), calculLHybridFunc(),
      calculMHybridFunc(), discretDynamicsHybrid(), ownsHybridControls(false)
    {
    isHybrid = true;
    if (!hybridTransistionControls)
	{
	ownsHybridControls = true;
	if (cp.DIM_HT > 0)
	    {
	    spdlog::info("[HybridSysDyn] : Starting build of ControlGrid for hybrid controls");
	    hybridTransistionControls = new ControlGrid(cp.DIM_HT, cp.LIMINF_HT, cp.LIMSUP_HT, cp.NBPOINTS_HT);
	    }
	else
	    {
	    hybridTransistionControls = new ControlGrid();
	    }
	}

    HybridSysDyn::initializeMethods(SP);
    }

HybridSysDyn::~HybridSysDyn()
    {
    if (ownsHybridControls)
	{
	delete hybridTransistionControls;
	}
    hybridTransistionControls = nullptr;
    }

void HybridSysDyn::initializeMethods(const systemParams &SP)
    {
    SimpleSysDyn::initializeMethods(SP);

    dynamics_hybrid_c = SP.DYNAMICS_HYBRID_C;
    dynamics_hybrid_d = SP.DYNAMICS_HYBRID_D;
    resetmap_hybrid = SP.RESET_MAP_HYBRID;
    constraintsXU_hybrid = SP.CONSTR_XU_HYBRID;
    constraintsX_hybrid = SP.CONSTR_X_HYBRID;
    localDynBounds_hybrid = SP.LOCAL_DYN_BOUNDS_HYBRID;
    jacobian_hybrid = SP.JACOBIAN_HYBRID;

    switch (computeLC)
	{
    case ANALYTICAL:
	calculLHybridFunc = [this](const double *xc, const unsigned long long int *xd) { return returnL_local_ana_hybrid(xc, xd); };
	break;
    case ANALYTICAL_CALC:
	calculLHybridFunc = [this](const double *xc, const unsigned long long int *xd) { return calculL_local_ana_hybrid(xc, xd); };
	break;
    case NUMERICAL_CALC:
	calculLHybridFunc = [this](const double *xc, const unsigned long long int *xd) { return calculL_local_num_hybrid(xc, xd); };
	break;
	}

    switch (computeMF)
	{
    case ANALYTICAL:
	calculMHybridFunc = [this](const double *xc, const unsigned long long int *xd) { return returnMF_local_ana_hybrid(xc, xd); };
	break;
    case ANALYTICAL_CALC:
	calculMHybridFunc = [this](const double *xc, const unsigned long long int *xd) { return calculMF_local_ana_hybrid(xc, xd); };
	break;
    case NUMERICAL_CALC:
	calculMHybridFunc = [this](const double *xc, const unsigned long long int *xd) { return calculMF_local_num_hybrid(xc, xd); };
	break;
	}

    switch (discretisation)
	{
    case NO_DISCRETIZATION_SCHEME:
	discretDynamicsHybrid = [this](const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	    unsigned long long int *resd, double *tempReset, unsigned long long int *tempResD, double rho) { FDiscret_hybrid(xc, xd, uc, ud, resc, resd, tempReset, tempResD, rho); };
	break;
    case EL:
	discretDynamicsHybrid = [this](const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	    unsigned long long int *resd, double *tempReset, unsigned long long int *tempResD, double rho) { FDiscretEuler_hybrid(xc, xd, uc, ud, resc, resd, tempReset, tempResD, rho); };
	break;
    case RK2:
	discretDynamicsHybrid = [this](const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	    unsigned long long int *resd, double *tempReset, unsigned long long int *tempResD, double rho) { FDiscretRK2_hybrid(xc, xd, uc, ud, resc, resd, tempReset, tempResD, rho); };
	break;
    case RK4:
	discretDynamicsHybrid = [this](const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	    unsigned long long int *resd, double *tempReset, unsigned long long int *tempResD, double rho) { FDiscretRK4_hybrid(xc, xd, uc, ud, resc, resd, tempReset, tempResD, rho); };
	break;
    default:
	discretDynamicsHybrid = [this](const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	    unsigned long long int *resd, double *tempReset, unsigned long long int *tempResD, double rho) { FDiscret_hybrid(xc, xd, uc, ud, resc, resd, tempReset, tempResD, rho); };
	break;
	}
    }

double* HybridSysDyn::getLimInfHybrid()
    {
    return hybridTransistionControls->GetLimInf();
    }

double* HybridSysDyn::getLimSupHybrid()
    {
    return hybridTransistionControls->GetLimSup();
    }

double* HybridSysDyn::getStepHybrid()
    {
    return hybridTransistionControls->GetStep();
    }

unsigned long long int HybridSysDyn::getDimHybrid() const
    {
    return hybridTransistionControls->GetDim();
    }

unsigned long long int* HybridSysDyn::getNbPointsHybrid()
    {
    return hybridTransistionControls->GetNbPoints();
    }

unsigned long long int HybridSysDyn::getTotalNbPointsHybrid() const
    {
    return hybridTransistionControls->GetTotalNbPoints();
    }

double** HybridSysDyn::getHybridCoords() const
    {
    return hybridTransistionControls->GetControlCoords();
    }

unsigned long long int** HybridSysDyn::getHybridIntCoords()
    {
    return hybridTransistionControls->GetControlIntCoords();
    }

double HybridSysDyn::calculRho_local_hybrid(const double *xc, const unsigned long long int *xd) const
    {
    if (dynType == DC)
	{
	return 1.0;
	}
    double h = grid->maxStep;
    double LL = calculLHybridFunc ? calculLHybridFunc(xc, xd) : returnL_local_ana_hybrid(xc, xd);
    double MFF = calculMHybridFunc ? calculMHybridFunc(xc, xd) : returnMF_local_ana_hybrid(xc, xd);
    if (MFF * LL < 2.0 * h)
	{
	MFF = 1.0;
	LL = 1.0;
	}
    return sqrt((2.0 * h) / (LL * MFF));
    }

void HybridSysDyn::FDiscret_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
    unsigned long long int *resd, double *tempResC, unsigned long long int *tempResD, double rho) const
    {
    (void)rho;
    (*resetmap_hybrid)(xc, xd, ud, tempResC, tempResD);
    (*dynamics_hybrid_d)(tempResC, tempResD, ud, resd);
    (*dynamics_hybrid_c)(tempResC, tempResD, uc, resc);
    }

void HybridSysDyn::FDiscretEuler_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
    unsigned long long int *resd, double *tempResC, unsigned long long int *tempResD, double rho) const
    {
    FDiscret_hybrid(xc, xd, uc, ud, tempResC, tempResD, resc, resd, rho);
    for (int i = 0; i < dimS_hc; i++)
	{
	resc[i] = tempResC[i] + rho * dynSignFactor * resc[i];
	}
    grid->periodizePoint(resc);
    }

void HybridSysDyn::FDiscretRK2_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
    unsigned long long int *resd, double *tempResC, unsigned long long int *tempResD, double rho) const
    {
    double *Fx = new double[dimS_hc];
    double *Fres = new double[dimS_hc];
    FDiscret_hybrid(xc, xd, uc, ud, tempResC, tempResD, Fx, resd, rho);
    for (int i = 0; i < dimS_hc; i++)
	{
	resc[i] = tempResC[i] + rho * dynSignFactor * Fx[i];
	}
    grid->periodizePoint(resc);
    (*dynamics_hybrid_c)(resc, tempResD, uc, Fres);
    for (int i = 0; i < dimS_hc; i++)
	{
	resc[i] = tempResC[i] + 0.5 * rho * dynSignFactor * (Fx[i] + Fres[i]);
	}
    grid->periodizePoint(resc);
    delete[] Fx;
    delete[] Fres;
    }

void HybridSysDyn::FDiscretRK4_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
    unsigned long long int *resd, double *tempResC, unsigned long long int *tempResD, double rho) const
    {
    double *ki = new double[dimS_hc];
    double *y = new double[dimS_hc];
    FDiscret_hybrid(xc, xd, uc, ud, tempResC, tempResD, ki, resd, rho);
    for (int i = 0; i < dimS_hc; i++)
	{
	resc[i] = tempResC[i] + rho * dynSignFactor * ki[i] / 6.0;
	y[i] = tempResC[i] + 0.5 * rho * dynSignFactor * ki[i];
	}
    grid->periodizePoint(y);
    (*dynamics_hybrid_c)(y, tempResD, uc, ki);
    for (int i = 0; i < dimS_hc; i++)
	{
	y[i] = tempResC[i] + 0.5 * rho * dynSignFactor * ki[i];
	resc[i] = resc[i] + rho * dynSignFactor * ki[i] / 3.0;
	}
    grid->periodizePoint(y);
    (*dynamics_hybrid_c)(y, tempResD, uc, ki);
    for (int i = 0; i < dimS_hc; i++)
	{
	y[i] = tempResC[i] + rho * dynSignFactor * ki[i];
	resc[i] = resc[i] + rho * dynSignFactor * ki[i] / 3.0;
	}
    grid->periodizePoint(y);
    (*dynamics_hybrid_c)(y, tempResD, uc, ki);
    for (int i = 0; i < dimS_hc; i++)
	{
	resc[i] = resc[i] + rho * dynSignFactor * ki[i] / 6.0;
	}
    grid->periodizePoint(resc);
    delete[] ki;
    delete[] y;
    }

double HybridSysDyn::calculL_local_num_hybrid(const double *xc, const unsigned long long int *xd) const
    {
    double *xTempL = new double[dimS_hc];
    double *FXmoinsHL = new double[dimS_hc];
    double *FXplusHL = new double[dimS_hc];
    double *infX = grid->limInf;
    double *pasX = grid->step;
    double *supX = grid->limSup;
    double L1 = 0;
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();
    for (int i = 0; i < dimS_hc; i++)
	{
	xTempL[i] = xc[i];
	}
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
	{
	for (int j = 0; j < dimS_hc; j++)
	    {
	    bool test = false;
	    xTempL[j] = xTempL[j] - pasX[j];
	    if ((xTempL[j] <= supX[j]) && (xTempL[j] >= infX[j]))
		{
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
		(*dynamics_hybrid_c)(xTempL, xd, controlCoords[nu], FXplusHL);
		xTempL[j] = xTempL[j] - pasX[j];
		}
	    else
		{
		xTempL[j] = xTempL[j] - pasX[j];
		(*dynamics_hybrid_c)(xTempL, xd, controlCoords[nu], FXplusHL);
		test = true;
		}
	    for (int k = 0; k < dimS_hc; k++)
		{
		FXmoinsHL[k] = fabs(FXmoinsHL[k] - FXplusHL[k]);
		FXmoinsHL[k] /= test ? pasX[k] : (2.0 * pasX[k]);
		L1 = max(L1, FXmoinsHL[k]);
		}
	    }
	}
    delete[] xTempL;
    delete[] FXplusHL;
    delete[] FXmoinsHL;
    return max(L1, lfunc_L);
    }

double HybridSysDyn::calculL_local_ana_hybrid(const double *xc, const unsigned long long int *xd) const
    {
    double **localJac = new double*[dimS_hc];
    for (int i = 0; i < dimS_hc; i++)
	{
	localJac[i] = new double[dimS_hc];
	}
    double L1 = 0;
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
	{
	(*jacobian_hybrid)(xc, xd, controlCoords[nu], localJac);
	for (int k = 0; k < dimS_hc; k++)
	    {
	    for (int j = 0; j < dimS_hc; j++)
		{
		L1 = max<double>(L1, abs(localJac[k][j]));
		}
	    }
	}
    for (int i = 0; i < dimS_hc; i++)
	{
	delete[] localJac[i];
	}
    delete[] localJac;
    return max(L1, lfunc_L);
    }

double HybridSysDyn::returnL_local_ana_hybrid(const double *xc, const unsigned long long int *xd) const
    {
    (void)xc;
    (void)xd;
    return max(L, lfunc_L);
    }

double HybridSysDyn::calculMF_local_num_hybrid(const double *xc, const unsigned long long int *xd) const
    {
    double *localImage = new double[dimS_hc];
    double MF1 = 0.0;
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
	{
	(*dynamics_hybrid_c)(xc, xd, controlCoords[nu], localImage);
	double normeImage = 0.0;
	for (int k = 0; k < dimS_hc; k++)
	    {
	    normeImage = max(normeImage, abs(localImage[k]));
	    }
	MF1 = max(MF1, normeImage);
	}
    delete[] localImage;
    return max(MF1, lfunc_MF);
    }

double HybridSysDyn::calculMF_local_ana_hybrid(const double *xc, const unsigned long long int *xd) const
    {
    double *localImage = new double[dimS_hc];
    (*localDynBounds_hybrid)(xc, xd, localImage);
    double normeImage = 0.0;
    for (int k = 0; k < dimS_hc; k++)
	{
	normeImage = max(normeImage, abs(localImage[k]));
	}
    delete[] localImage;
    return max(normeImage, lfunc_MF);
    }

double HybridSysDyn::returnMF_local_ana_hybrid(const double *xc, const unsigned long long int *xd) const
    {
    (void)xc;
    (void)xd;
    return max(MF, lfunc_MF);
    }
