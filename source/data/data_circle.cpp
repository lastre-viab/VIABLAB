//
// Created by adesi on 05/03/2026.
//
#include "../include/utilities.h"
#include "../include/ParametersManager.h"
#include "../include/defs.h"
#include <cstdlib>


extern "C" {
std::string paramsFile = "CircleViab.json";

//------------------------------------------------------------------------------------------------------
//  Model description source for the test problem
//------------------------------------------------------------------------------------------------------
using std::abs;
using std::max;


double a = 0.7;
double c = 1.5;
double eps = 0.25;
double tol = 0.001;

void loadModelData(const ParametersManager* PM)
{
    spdlog::info("[LoadModelData] : STARt");
    const modelParams* modelParams = PM->getModelParameters();
    a = modelParams->getDouble("CIRCLE_PARAM");
    c = modelParams->getDouble("CONSTR_PARAM");
    eps = modelParams->getDouble("EPS_PARAM");
    tol = modelParams->getDouble("TOL_PARAM");

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

void dynamics(const double* x, const double* u, double* image)
{
    double nx = norm(x);
    double alpha = (nx - a) / nx;
    image[0] = alpha * x[0] + x[1];
    image[1] = -x[0] + alpha * x[1];
}

/*!
 * jacobian matrix of the dynamics
 */
void jacobian(const double* x, const double* u, double** image)
{
    double nx = norm(x);
    double alpha = (nx - a) / nx;
    image[0][0] = a * x[0] * x[0] / (nx * nx * nx) + alpha;
    image[0][1] = 1 + x[0] * a * x[1] / (nx * nx * nx);
    image[1][0] = -1 + x[0] * a * x[1] / (nx * nx * nx);
    image[1][1] = a * x[1] * x[1] / (nx * nx * nx) + alpha;
}

/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
void localDynBounds(const double* x, double* res)
{
    double nx = norm(x);
    double alpha = (nx - a) / nx;
    res[0] = abs(x[1] + alpha * x[0]);
    res[1] = abs(-x[0] + alpha * x[1]);
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
