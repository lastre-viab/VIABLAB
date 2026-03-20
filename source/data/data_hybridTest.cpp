//
// Created by adesi on 12/11/2025.
//
#include "../include/utilities.h"
#include "../include/ParametersManager.h"
#include "../include/defs.h"
#include <cstdlib>


extern "C" {
std::string paramsFile = "hybridTest.json";

//------------------------------------------------------------------------------------------------------
//  Model description source for the test problem
//------------------------------------------------------------------------------------------------------
using std::abs;
using std::max;


double a = 0.7;
double c = 1.5;
double eps = 0.25;
double tol = 0.001;
double reset = 0.5;
    double resetAlpha1 = 0.5;
    double resetAlpha2 = -0.5;
    double resetBeta1 = 1.0;
    double resetBeta2 = 1.0;
    int resetNbLines = 1;
double * resetAlphas;
    double * resetBetas;
bool withReset = false;

void loadModelData(const ParametersManager* PM)
{
    spdlog::info("[LoadModelData] : STARt");
    const modelParams* modelParams = PM->getModelParameters();
    a = modelParams->getDouble("CIRCLE_PARAM");
    c = modelParams->getDouble("CONSTR_PARAM");
    eps = modelParams->getDouble("EPS_PARAM");
    tol = modelParams->getDouble("TOL_PARAM");
    reset = modelParams->getDouble("RESET_PARAM");
    resetAlpha1 = modelParams->getDouble("RESET_ALPHA1_PARAM");
    resetAlpha2 = modelParams->getDouble("RESET_ALPHA2_PARAM");
    resetBeta1 = modelParams->getDouble("RESET_BETA1_PARAM");
    resetBeta2 = modelParams->getDouble("RESET_BETA2_PARAM");
    resetAlpha1 = modelParams->getDouble("RESET_ALPHA1_PARAM");
    resetNbLines = modelParams->getInt("RESET_NBLINES");
    withReset = modelParams->getBool("WITH_RESET_PARAM");
    resetAlphas = new double[resetNbLines];
    resetBetas = new double[resetNbLines];
    if (resetNbLines >= 1)
    {
        resetAlphas[0] = resetAlpha1;
        resetBetas[0] = resetBeta1;
    }
    if (resetNbLines >= 2)
    {
        resetAlphas[1] = resetAlpha2;
        resetBetas[1] = resetBeta2;
    }

    spdlog::info("[LoadModelData] : model params ok");
    spdlog::info("[LoadModelData] : Finished");
}


double norm(const double* x)
{
    double res = 0;
    for (int i = 0; i < 2; i++)
    {
        res += x[i] * x[i];
    }
    return sqrt(res);
}

void dynamics_hybrid_c(const double* x, const unsigned long long int* xd, const double* u, double* image)
{
    double nx = norm(x);
    double alpha = (nx - a) / nx;
    image[0] = alpha * x[0] + x[1];
    image[1] = - x[0] +  alpha * x[1];
}

void jacobian_hybrid(const double* x, const unsigned long long int* xd, const double* u, double** image)
{
    double nx = norm(x);
    double alpha = (nx - a) / nx;
    image[0][0] = a * x[0] * x[0] / (nx * nx * nx) + alpha;
    image[0][1] = 1 + x[0] * a * x[1] / (nx * nx * nx);
    image[1][0] = -1 + x[0] * a * x[1] / (nx * nx * nx);
    image[1][1] = a * x[1] * x[1] / (nx * nx * nx) + alpha;
}

void localDynBounds_hybrid(const double* x, const unsigned long long int* xd, double* res)
{
    double nx = norm(x);
    double alpha = (nx - a) / nx;
    res[0] = abs(x[1] + alpha * x[0]);
    res[1] = abs(-x[0] + alpha * x[1]);
}

void dynamics_hybrid_d(const double* x, const unsigned long long int* xd, const unsigned long long int* u,
                       unsigned long long int* image)
{
    image[0] = xd[0];
}

bool isInResetSet(const double* x, const unsigned long long int* currentDiscreteState)
{
    if (!withReset)
    {
        return false;
    }
    bool isInReset = false;
    for (int i = 0; i < resetNbLines; i++)
    {
        isInReset |= abs(x[0] + resetAlphas[i] * x[1] - resetBetas[i]) < tol;
    }
    return isInReset;
}

void resetMap_hybrid(const double* x, const unsigned long long int* currentDiscreteState,
                     const unsigned long long int* discreteControl, double* image,
                     unsigned long long int* nextDiscreteState)
{
    if (isInResetSet(x, currentDiscreteState))
    {
        image[0] = x[0] - reset;
        image[1] = x[1];
    }
    else
    {
        image[0] = x[0];
        image[1] = x[1];
    }
    nextDiscreteState[0] = currentDiscreteState[0];
}


double constraintsX(const double* x)
{
    if (norm(x) < eps)
    {
        return PLUS_INF;
    }
    if ((x[0] < -c) || (x[0] > c))
    {
        return PLUS_INF;
    }
    if ((x[1] < -c) || (x[1] > c))
    {
        return PLUS_INF;
    }
    return 1.0;
}
}
