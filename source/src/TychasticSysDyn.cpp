#include "../include/TychasticSysDyn.h"

#include <algorithm>
#include <cmath>

#include <spdlog/spdlog.h>

using std::abs;
using std::max;

TychasticSysDyn::TychasticSysDyn()
    : SimpleSysDyn(), tyches(new ControlGrid()), dynamics_tych(nullptr), dynamics_tych_fd(nullptr), constraintsXV_tych(nullptr), lFunc_tych(nullptr),
      lFunc_tych_fd(nullptr), mFunc_tych(nullptr), jacobian_tych(nullptr), discretDynamicsTych(), ownsTyches(true)
    {
    isTychastic = true;
    }

TychasticSysDyn::TychasticSysDyn(const systemParams &SP, int stateDim, const controlParams &cp, Grid *refGrid, ControlGrid *controlGrid, ControlGrid *tychesGrid)
    : SimpleSysDyn(SP, stateDim, cp, refGrid, controlGrid), tyches(tychesGrid), dynamics_tych(nullptr), dynamics_tych_fd(nullptr), constraintsXV_tych(nullptr),
      lFunc_tych(nullptr), lFunc_tych_fd(nullptr), mFunc_tych(nullptr), jacobian_tych(nullptr), discretDynamicsTych(), ownsTyches(false)
    {
    isTychastic = (cp.DIM_TY > 0);
    if (!tyches)
	{
	ownsTyches = true;
	if (cp.DIM_TY > 0)
	    {
	    spdlog::info("[TychasticSysDyn] : Starting build of ControlGrid for tyches");
	    tyches = new ControlGrid(cp.DIM_TY, cp.LIMINF_TY, cp.LIMSUP_TY, cp.NBPOINTS_TY);
	    }
	else
	    {
	    tyches = new ControlGrid();
	    }
	}

    TychasticSysDyn::initializeMethods(SP);
    }

TychasticSysDyn::~TychasticSysDyn()
    {
    if (ownsTyches)
	{
	delete tyches;
	}
    tyches = nullptr;
    }

void TychasticSysDyn::initializeMethods(const systemParams &SP)
    {
    SimpleSysDyn::initializeMethods(SP);

    dynamics_tych = SP.DYNAMICS_TYCH;
    dynamics_tych_fd = SP.DYNAMICS_TYCH_FD;
    constraintsXV_tych = SP.CONSTR_XV_TYCH;
    lFunc_tych = SP.L_FUNC_TYCH;
    lFunc_tych_fd = SP.L_FUNC_TYCH_FD;
    mFunc_tych = SP.M_FUNC_TYCH;
    jacobian_tych = SP.JACOBIAN_TYCH;

    if (computeLC == ANALYTICAL_CALC)
	{
	calculLFunc = [this](const double *x) { return calculL_local_ana_tych(x); };
	}
    else if (computeLC == NUMERICAL_CALC)
	{
	calculLFunc = [this](const double *x) { return calculL_local_num_tych(x); };
	}

    if (computeMF == NUMERICAL_CALC)
	{
	calculMFunc = [this](const double *x) { return calculMF_local_num_tych(x); };
	}

    switch (discretisation)
	{
    case NO_DISCRETIZATION_SCHEME:
	discretDynamicsTych = [this](const double *x, const double *u, const double *v, double *res, double rho) { FDiscret_tych(x, u, v, res, rho); };
	break;
    case EL:
	discretDynamicsTych = [this](const double *x, const double *u, const double *v, double *res, double rho) { FDiscretEuler_tych(x, u, v, res, rho); };
	break;
    case RK2:
	discretDynamicsTych = [this](const double *x, const double *u, const double *v, double *res, double rho) { FDiscretRK2_tych(x, u, v, res, rho); };
	break;
    case RK4:
	discretDynamicsTych = [this](const double *x, const double *u, const double *v, double *res, double rho) { FDiscretRK4_tych(x, u, v, res, rho); };
	break;
    default:
	discretDynamicsTych = [this](const double *x, const double *u, const double *v, double *res, double rho) { FDiscret_tych(x, u, v, res, rho); };
	break;
	}
    }

void TychasticSysDyn::getTychasticImage(const double *x, const double *u, const double *v, double *imageVect, double rho) const
    {
    if (discretDynamicsTych)
	{
	discretDynamicsTych(x, u, v, imageVect, rho);
	}
    }

double* TychasticSysDyn::getLimInfTy()
    {
    return tyches->GetLimInf();
    }

double* TychasticSysDyn::getLimSupTy()
    {
    return tyches->GetLimSup();
    }

double* TychasticSysDyn::getStepTy()
    {
    return tyches->GetStep();
    }

unsigned long long int TychasticSysDyn::getDimTy() const
    {
    return tyches->GetDim();
    }

unsigned long long int* TychasticSysDyn::getNbPointsTy()
    {
    return tyches->GetNbPoints();
    }

unsigned long long int TychasticSysDyn::getTotalNbPointsTy() const
    {
    return tyches->GetTotalNbPoints();
    }

double** TychasticSysDyn::getTychCoords() const
    {
    return tyches->GetControlCoords();
    }

unsigned long long int** TychasticSysDyn::getTychIntCoords()
    {
    return tyches->GetControlIntCoords();
    }

void TychasticSysDyn::FDiscret_tych(const double *x, const double *u, const double *v, double *res, double rho) const
    {
    (*dynamics_tych)(x, u, v, res);
    grid->periodizePoint(res);
    }

void TychasticSysDyn::FDiscretEuler_tych(const double *x, const double *u, const double *v, double *res, double rho) const
    {
    int i;

    (*dynamics_tych)(x, u, v, res);

    for (i = 0; i < dimS; i++)
	{
	res[i] = x[i] + rho * dynSignFactor * res[i];
	}

    grid->periodizePoint(res);
    }

void TychasticSysDyn::FDiscretRK2_tych(const double *x, const double *u, const double *v, double *res, double rho) const
    {
    int i;
    double *Fx, *Fres;
    Fx = new double[dimS];
    Fres = new double[dimS];
    (*dynamics_tych)(x, u, v, Fx);

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

void TychasticSysDyn::FDiscretRK4_tych(const double *x, const double *u, const double *v, double *res, double rho) const
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

    (*dynamics_tych)(y, u, v, ki);

    for (i = 0; i < dimS; i++)
	{
	y[i] = x[i] + 0.5 * rho * dynSignFactor * ki[i];
	res[i] = res[i] + rho * dynSignFactor * ki[i] / 3.0;
	}

    grid->periodizePoint(y);

    (*dynamics_tych)(y, u, v, ki);

    for (i = 0; i < dimS; i++)
	{
	y[i] = x[i] + rho * dynSignFactor * ki[i];
	res[i] = res[i] + rho * dynSignFactor * ki[i] / 3.0;
	}

    grid->periodizePoint(y);

    (*dynamics_tych)(y, u, v, ki);

    for (i = 0; i < dimS; i++)
	{
	res[i] = res[i] + rho * dynSignFactor * ki[i] / 6.0;
	}

    grid->periodizePoint(res);
    delete[] ki;
    delete[] y;
    }

double TychasticSysDyn::calculL_local_num_tych(const double *x) const
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
		if ((xTempL[j] <= supX[j]) && (xTempL[j] >= infX[j]))
		    {
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

double TychasticSysDyn::calculL_local_ana_tych(const double *x) const
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
    unsigned long long int totalNbPointsTych = tyches->GetTotalNbPoints();

    double **controlCoords = controls->GetControlCoords();
    double **tychCoords = tyches->GetControlCoords();
    for (unsigned long long int nu = 0; nu < totalNbPointsC; nu++)
	{
	for (unsigned long long int nv = 0; nv < totalNbPointsTych; nv++)
	    {

	    (*jacobian_tych)(x, controlCoords[nu], tychCoords[nv], localJac);
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
	}

    for (int i = 0; i < dimS; i++)
	{
	delete[] localJac[i];
	}
    delete[] localJac;
    return max(L1, lfunc_L);
    }

double TychasticSysDyn::calculMF_local_num_tych(const double *x) const
    {
    double *localImage = new double[dimS];

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
	    (*dynamics_tych)(x, controlCoords[nu], tychCoords[nv], localImage);
	    normeImage = 0.0;
	    for (int k = 0; k < dimS; k++)
		{
		normeImage = max(normeImage, abs(localImage[k]));
		}
	    MF1 = max(MF1, normeImage);
	    }
	}
    delete[] localImage;
    return max(MF1, lfunc_MF);
    }

