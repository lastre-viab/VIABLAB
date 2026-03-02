//
// Created by adesi on 12/11/2025.
//
#include "../include/utilities.h"
#include "../include/ParametersManager.h"
#include "../include/defs.h"
#include <cstdlib>

#include "Exemple_multiDim_data.h"

extern "C" {
std::string paramsFile = "hybridTest.json";

//------------------------------------------------------------------------------------------------------
//  Model description source for the test problem
//------------------------------------------------------------------------------------------------------


void loadModelData(const ParametersManager* PM)
{
}

double a = 1.0;
double c = 2.0;
double eps = 0.1;

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
    double n = norm(x);
    image[0] = x[1] + (n - a) * x[0] / n;
    image[1] = -x[0] + (n - a) * x[1] / n;
}

void dynamics_hybrid_d(const unsigned long long int* x, const unsigned long long int* u, unsigned long long int* image)
{

    image[0] = x[0]; // nouveau niveau
    image[1] = x[1]; // n
}


void resetMap_hybrid(const double* x, const unsigned long long int* currentDiscreteState,
                     const unsigned long long int* nextDiscreteState,
                     const unsigned long long int* discreteControl, double* image)
{
    image[0] = x[0] - 1;
    image[1] = x[1];
}

// Constraints

// Check M[](t) constraint: M[](t) in {0,E(t),I+(t)}
bool CheckMConstraints(const unsigned long long int M, const unsigned long long int E,
                       const unsigned long long int Iplus)
{
    if (E == 3)
    {
        return true;
    }

    if (M != 0 && M != E && M != Iplus)
    {
        return false;
    }
    return true;
}

// Check I+(t) constraint
bool CheckIConstraint(const unsigned long long int E, const unsigned long long int Iplus)
{
    if (Iplus < 0)
    {
        return false;
    }
    // si E(t)=0，I+(t)=0,1,2
    if (E == 0)
    {
        return (Iplus == 0 || Iplus == 1 || Iplus == 2);
    }
    // si E(t)=1, I+(t)=0,2
    else if (E == 1)
    {
        return (Iplus == 0 || Iplus == 2);
    }
    // si E(t)=2, I+(t)=0,1
    else if (E == 2)
    {
        return (Iplus == 0 || Iplus == 1);
    }
    // si E(t)=3, I+(t)=0
    else if (E == 3)
    {
        return (Iplus == 0);
    }
    else
    {
        return false;
    }
}


double constraintsXU_hybrid(const double* x, const unsigned long long int* xd, const double* u,
                            const unsigned long long int* ud)
{
    // discrete control : ud[0] = surcreusement
    // discrete control : ud[1] = investissement
    // x[0] = capital
    // x[1] = niveau de la nappe
    // discreteState xd[0] = niveau puit
    // DiscreteState xd[1] = equipement
    // cont control M : mode irrigation = u[1]
    // cont control S : spec = u[0]

    getIntControlCoords(u, continuousControlIntCoords);

    if (!CheckMConstraints(continuousControlIntCoords[1], xd[1], ud[1]))
    {
        //spdlog::info("[contraints XU ] : M pas OK");
        return PLUS_INF;
    }

    // Check S(t) constraint: S(t) in [P(t)/5, 8]
    if (xd[0] > ud[0])
    {
        return PLUS_INF;
    }

    // Check I+(t) constraint
    if (!CheckIConstraint(xd[1], ud[1]))
    {
        return PLUS_INF;
    }

    // Check the constraint: h(t)+BES_f(t)*a+BES_c(t)*a<=S(t)*5
    double BES_f = BesoinEnEau(continuousControlIntCoords);
    // double BES_c = BesoinEnEau(M_c, Sp_c, BES_PARAM);

    if (x[1] + BES_f * BES_para > ud[0] * 5)
    {
        // spdlog::info("[contraints XU ] :niveau de la nappe {}", x[1]);
        // spdlog::info("[contraints XU ] :besoin eau {}", BES_f * BES_para);
        // spdlog::info("[contraints XU ] :niveau puit {}", ud[0] * 5);
        return PLUS_INF;
    }

    // Check the constraint: h>=30 and M_k=0, then Sp_k=0

    if (x[1] >= 30 && continuousControlIntCoords[1] == 0 && continuousControlIntCoords[0] != 0)
    {
        // spdlog::info("[contraints XU ] : Pas droit d'arroser à la main");
        return PLUS_INF;
    }


    // Check the constraint: C(t)

    double CS = CoutSurcreuse(xd[0], ud[0]);
    double CI = CoutInvest(ud[1]);
    double PV_f = PrixVente(continuousControlIntCoords);

    double CP_f = CoutProduction(continuousControlIntCoords);

    // Constraint 1: γ*C(t)-CS(t)-CI(t)-CP_f(t)-CP_c(t)-FC>=0
    if (gamma_para * (x[0] - CS - CI) - CP_f - FIXED_COSTSR < 0)
    {
        //	 spdlog::info("[contraints XU ] : pas assez de capital");
        //	 spdlog::info("[contraints XU ] : capital {}", x[0]);
        //	 spdlog::info("[contraints XU ] : cout creuser {}", CS);
        //	 spdlog::info("[contraints XU ] : cout investir {}", CI);
        //	 spdlog::info("[contraints XU ] : cout produire {}", CP_f);
        //	 spdlog::info("[contraints XU ] : cout fixe {}", FIXED_COSTSR);
        return PLUS_INF;
    }
    if (gamma_para* x[0]
    <
    0
    )
    {
        return PLUS_INF;
    }

    // If all constraints are satisfied, return 1.0 to indicate feasibility
    // spdlog::info("[contraints XU ] : OK");
    return 1.0;
}
}
