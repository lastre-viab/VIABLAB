#include "../include/SimpleSysDyn.h"

#include <algorithm>
#include <cmath>

#include <spdlog/spdlog.h>

using std::abs;
using std::max;

SimpleSysDyn::SimpleSysDyn()
    : SysdynBase(), controls(nullptr), dynamics(nullptr), dynamics_fd(nullptr), localDynBounds(nullptr), jacobian(nullptr), constraintsXU(nullptr),
      constraintsXU_fd(nullptr), constraintsX(nullptr), constraintsX_fd(nullptr), controlEligibilityForTraj_fd(nullptr), dynConstraintsForTraj(nullptr),
      target(nullptr), target_fd(nullptr), lFunc(nullptr), lFunc_fd(nullptr), muFunc_fd(nullptr), mFunc(nullptr), calculLFunc(), calculMFunc(),
      discretDynamics(), fd_dyn_type(0), timeHorizon(0.0), retroFileName(), discretisation(0), dynType(DD), lfunc_L(0.0), lfunc_MF(0.0),
      dynSignFactor(1.0), ownsControls(false)
    {
    }

SimpleSysDyn::SimpleSysDyn(const systemParams &SP, int stateDim, const controlParams &cp, Grid *refGrid, ControlGrid *control_grid)
    : SysdynBase(refGrid, stateDim), controls(control_grid), dynamics(nullptr), dynamics_fd(nullptr), localDynBounds(SP.LOCAL_DYN_BOUNDS),
      jacobian(nullptr), constraintsXU(nullptr), constraintsXU_fd(nullptr), constraintsX(nullptr), constraintsX_fd(nullptr),
      controlEligibilityForTraj_fd(nullptr), dynConstraintsForTraj(nullptr), target(nullptr), target_fd(nullptr), lFunc(nullptr), lFunc_fd(nullptr),
      muFunc_fd(nullptr), mFunc(nullptr), calculLFunc(), calculMFunc(), discretDynamics(), fd_dyn_type(0), timeHorizon(0.0), retroFileName(),
      discretisation(0), dynType(DD), lfunc_L(0.0), lfunc_MF(0.0), dynSignFactor(1.0), ownsControls(false)
    {
	timeStepFactor = SP.TIME_STEP_FACTOR;
    if (!controls)
	{
	ownsControls = true;
	if (cp.DIMC > 0)
	    {
	    spdlog::info("[SimpleSysDyn] : Starting build of ControlGrid");
	    controls = new ControlGrid(cp.DIMC, cp.LIMINFC, cp.LIMSUPC, cp.NBPOINTSC);
	    }
	else
	    {
	    controls = new ControlGrid();
	    }
	}

    isTychastic = false;
    isHybrid = false;

    globalTimeStep = SP.globDeltat;

    computeMF = SP.COMPUTE_MF;
    computeLC = SP.COMPUTE_LC;

    MF = SP.MF;
    L = SP.LIP;

    lfunc_L = SP.L_LIP;
    lfunc_MF = SP.L_MF;

    if (L == 0)
	{
	L = 1.0;
	}
    if (MF == 0)
	{
	MF = 1.0;
	}

    initializeMethods(SP);
    spdlog::info("[System] : Dynamic system initialization finished");
    }

SimpleSysDyn::~SimpleSysDyn()
    {
    if (ownsControls)
	{
	delete controls;
	}
    }

void SimpleSysDyn::initializeMethods(const systemParams &SP)
    {
    dynamics = SP.DYNAMICS;
    dynamics_fd = SP.DYNAMICS_FD;

    dynType = SP.DYN_TYPE;

    if (SP.DYN_TYPE == DD)
	{
	spdlog::debug(" DynSYS de type DD");
	if (SP.FD_DYN_TYPE == RETRO)
	    {
	    retroFileName = SP.RETRO_FILE_NAME;
	    }
	}
    fd_dyn_type = SP.FD_DYN_TYPE;
    constraintsXU = SP.CONSTR_XU;
    constraintsXU_fd = SP.CONSTR_XU_fd;
    constraintsX = SP.CONSTR_X;
    constraintsX_fd = SP.CONSTR_X_fd;
    controlEligibilityForTraj_fd = SP.CONTROL_ELIGIBILITY_FOR_TRAJ_fd;
    dynConstraintsForTraj = SP.DYN_CONSTR_FOR_TRAJ;
    target = SP.TARGET;
    target_fd = SP.TARGET_FD;

    jacobian = SP.JACOBIAN;

    lFunc = SP.L_FUNC;
    lFunc_fd = SP.L_FUNC_FD;
    muFunc_fd = SP.MU_FUNC_FD;
    mFunc = SP.M_FUNC;

    switch (computeMF)
	{
    case ANALYTICAL:
	calculMFunc = [this](const double *x) { return returnMF_local_ana(x); };
	break;
    case ANALYTICAL_CALC:
	calculMFunc = [this](const double *x) { return calculMF_local_ana(x); };
	break;
    case NUMERICAL_CALC:
	calculMFunc = [this](const double *x) { return calculMF_local_num(x); };
	break;
	}

    switch (computeLC)
	{
    case ANALYTICAL:
	calculLFunc = [this](const double *x) { return returnL_local_ana(x); };
	break;
    case ANALYTICAL_CALC:
	calculLFunc = [this](const double *x) { return calculL_local_ana(x); };
	break;
    case NUMERICAL_CALC:
	calculLFunc = [this](const double *x) { return calculL_local_num(x); };
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
	discretDynamics = [this](const double *x, const double *u, double *res, double rho) { FDiscret(x, u, res, rho); };
	break;
    case EL:
	discretDynamics = [this](const double *x, const double *u, double *res, double rho) { FDiscretEuler(x, u, res, rho); };
	break;
    case RK2:
	discretDynamics = [this](const double *x, const double *u, double *res, double rho) { FDiscretRK2(x, u, res, rho); };
	break;
    case RK4:
	discretDynamics = [this](const double *x, const double *u, double *res, double rho) { FDiscretRK4(x, u, res, rho); };
	break;
	}
    }

unsigned long long int** SimpleSysDyn::getControlIntCoords()
    {
    return controls->GetControlIntCoords();
    }

double* SimpleSysDyn::getLimSupC()
    {
    return controls->GetLimSup();
    }

double* SimpleSysDyn::getLimInfC()
    {
    return controls->GetLimInf();
    }

double* SimpleSysDyn::getStepC()
    {
    return controls->GetStep();
    }

unsigned long long int SimpleSysDyn::getDimC() const
    {
    return controls->GetDim();
    }

unsigned long long int* SimpleSysDyn::getNbPointsC()
    {
    return controls->GetNbPoints();
    }

unsigned long long int SimpleSysDyn::getTotalNbPointsC() const
    {
    return controls->GetTotalNbPoints();
    }

double** SimpleSysDyn::getControlCoords() const
    {
    return controls->GetControlCoords();
    }

double SimpleSysDyn::calculRho_local(const double *x) const
    {
    if (dynType == DC)
	{
	return 1.0;
	}

    double h = grid->maxStep;
    double LL = calculLFunc ? calculLFunc(x) : returnL_local_ana(x);
    double MFF = calculMFunc ? calculMFunc(x) : returnMF_local_ana(x);
  if (MFF * LL < 2.0 * h)
  {
    MFF = 1.0;
    LL = 1.0;
  }

  double rho1 = sqrt((2.0 * h) / (LL * MFF));
  return timeStepFactor * rho1;
    //return timeStepFactor * sqrt((h) / max(h *h, LL * MFF));
    }

int SimpleSysDyn::getFDDynType() const
    {
    return fd_dyn_type;
    }

std::string SimpleSysDyn::getRetroFileName() const
    {
    return retroFileName;
    }

DynType SimpleSysDyn::getDynType() const
    {
    return dynType;
    }

void SimpleSysDyn::setDynamicsForward()
    {
    dynSignFactor = 1.0;
    }

void SimpleSysDyn::setDynamicsBackward()
    {
    dynSignFactor = -1.0;
    }

void SimpleSysDyn::FDiscretEuler(const double *x, const double *u, double *res, double rho) const
    {
    int i;

    (*dynamics)(x, u, res);

    for (i = 0; i < dimS; i++)
	{
	res[i] = x[i] + rho * dynSignFactor * res[i];
	}

    grid->periodizePoint(res);
    }

void SimpleSysDyn::FDiscret(const double *x, const double *u, double *res, double rho) const
    {
    (*dynamics)(x, u, res);
    grid->periodizePoint(res);
    }

void SimpleSysDyn::FDiscretRK4(const double *x, const double *u, double *res, double rho) const
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

    (*dynamics)(y, u, ki);

    for (i = 0; i < dimS; i++)
	{
	y[i] = x[i] + 0.5 * rho * dynSignFactor * ki[i];
	res[i] = res[i] + rho * dynSignFactor * ki[i] / 3.0;
	}

    grid->periodizePoint(y);

    (*dynamics)(y, u, ki);

    for (i = 0; i < dimS; i++)
	{
	y[i] = x[i] + rho * dynSignFactor * ki[i];
	res[i] = res[i] + rho * dynSignFactor * ki[i] / 3.0;
	}

    grid->periodizePoint(y);

    (*dynamics)(y, u, ki);

    for (i = 0; i < dimS; i++)
	{
	res[i] = res[i] + rho * dynSignFactor * ki[i] / 6.0;
	}

    grid->periodizePoint(res);
    delete[] ki;
    delete[] y;
    }

void SimpleSysDyn::FDiscretRK2(const double *x, const double *u, double *res, double rho) const
    {
    int i;
    double *Fx, *Fres;
    Fx = new double[dimS];
    Fres = new double[dimS];
    (*dynamics)(x, u, Fx);

    for (i = 0; i < dimS; i++)
	{
	res[i] = x[i] + rho * dynSignFactor * Fx[i];
	}

    grid->periodizePoint(res);

    (*dynamics)(res, u, Fres);

    for (i = 0; i < dimS; i++)
	{
	res[i] = x[i] + 0.5 * rho * dynSignFactor *  (Fx[i] + Fres[i]);
	}

    grid->periodizePoint(res);
    delete[] Fx;
    delete[] Fres;
    }

double SimpleSysDyn::calculL_local_num(const double *x) const
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
	    if ((xTempL[j] <= supX[j]) && (xTempL[j] >= infX[j]))
		{
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

double SimpleSysDyn::calculL_local_ana(const double *x) const
    {
    int j, k;
    double **localJac = new double*[dimS];

    for (int i = 0; i < dimS; i++)
	{
	localJac[i] = new double[dimS];
	}
    double L1 = 0;
    double norme;
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
	{

	(*jacobian)(x, controlCoords[nu], localJac);
	norme = 0.;
	for (k = 0; k < dimS; k++)
	    {
	    for (j = 0; j < dimS; j++)
		{
		norme = max(norme, abs(localJac[k][j]));
		}
	    }
	L1 = max(L1, norme);
	}

    for (int i = 0; i < dimS; i++)
	{
	delete[] localJac[i];
	}
    delete[] localJac;
    return max(L1, lfunc_L);
    }

double SimpleSysDyn::returnL_local_ana(const double *x) const
    {
    return max(L, lfunc_L);
    }


double SimpleSysDyn::calculMF_local_num(const double *x) const
    {
    double *localImage = new double[dimS];

    double MF1 = 0.0;
    double normeImage;
    unsigned long long int totalNbPointsC = controls->GetTotalNbPoints();
    double **controlCoords = controls->GetControlCoords();

    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
	{
	(*dynamics)(x, controlCoords[nu], localImage);
	normeImage = 0.0;
	for (int k = 0; k < dimS; k++)
	    {
	    normeImage = max(normeImage, abs(localImage[k]));
	    }
	MF1 = max(MF1, normeImage);
	}
    delete[] localImage;
    return max(MF1, lfunc_MF);
    }

double SimpleSysDyn::calculMF_local_ana(const double *x) const
    {
    double *localImage = new double[dimS];
    double normeImage;
    (*localDynBounds)(x, localImage);

    normeImage = 0.0;
    for (int k = 0; k < dimS; k++)
	{
	normeImage = max(normeImage, abs(localImage[k]));
	}
    double MF1 = normeImage;
    delete[] localImage;
    return max(MF1, lfunc_MF);
    }

double SimpleSysDyn::returnMF_local_ana(const double *x) const
    {
    return max(MF, lfunc_MF);
    }

