/*
 * initParams.h
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
 *  Created on: reated on: 24 aug. 2023
 *      Author: Anna DESILLES
 */

#ifndef SRC_INITPARAMS_H_
#define SRC_INITPARAMS_H_

#include <dlfcn.h>

#include "../include/defs.h"
#include "../include/Params.h"
#include "../include/TrajectoryHelpers.h"
#include "../include/ParametersManager.h"

template<typename T>
T getUserSymbolOrDefaultTemplate(void *modelHandle, const char *symbolName, T defaultSymbol) {
    T symbol = (T) dlsym(modelHandle, symbolName);
    if (symbol == nullptr) {
        symbol = defaultSymbol;
        spdlog::info("Using default symbol for {}", symbolName);
    }
    else {
        spdlog::info("Succesfully loaded symbol : {}", symbolName);
    }
    return symbol;
}

template<typename T>
T getUserSymbol(void *modelHandle, const char *symbolName) {
    T symbol = (T) dlsym(modelHandle, symbolName);
    if (symbol == nullptr) {
        spdlog::info("No {} symbol found", symbolName);
    }
    else {
        spdlog::info("Succesfully loaded symbol : {}", symbolName);
    }
    return symbol;
}

#define getUserSymbolOrDefault(handle, f) getUserSymbolOrDefaultTemplate(handle, #f, &f)

ParametersManager* initParams(void *modelHandle, const string &paramsFile, gridParams &gp, algoViabiParams &avp, controlParams &cp, systemParams &sp, int nbOmpThreads)
    {

        ParametersManager *pm = new ParametersManager(&gp, &avp, &cp, &sp, nbOmpThreads, paramsFile, modelHandle);
    
	sp.DYNAMICS = getUserSymbolOrDefault(modelHandle, dynamics);

	sp.DYNAMICS_TYCH = getUserSymbolOrDefault(modelHandle, dynamics_tych);

    switch (sp.DYN_TYPE)
	{
    case CC:
        break;
    case DC:
    case DD:	
        sp.COMPUTE_LC = ANALYTICAL;
        sp.COMPUTE_MF = ANALYTICAL;
        
        sp.LIP = 1.0;
        sp.MF = 1.0;
        break;	
    case DH:
        sp.DYNAMICS = getUserSymbolOrDefault(modelHandle, dynamics_hybrid);
        break;
	}
	sp.CONSTR_X = getUserSymbolOrDefault(modelHandle, constraintsX);

	sp.CONSTR_XU = getUserSymbolOrDefault(modelHandle, constraintsXU);
    sp.CONSTR_XV_TYCH = getUserSymbolOrDefault(modelHandle, constraintsXV_tych);
	sp.JACOBIAN = getUserSymbolOrDefault(modelHandle, jacobian);
	sp.JACOBIAN_TYCH = getUserSymbolOrDefault(modelHandle, jacobian_tych);
	sp.LOCAL_DYN_BOUNDS = getUserSymbolOrDefault(modelHandle, localDynBounds);
    
	sp.M_FUNC = getUserSymbolOrDefault(modelHandle, m);
    sp.L_FUNC = getUserSymbolOrDefault(modelHandle, l);

	sp.M_FUNC_TYCH = getUserSymbolOrDefault(modelHandle, m_tych);
    sp.L_FUNC_TYCH = getUserSymbolOrDefault(modelHandle, l_tych);

	sp.TARGET = getUserSymbolOrDefault(modelHandle, target);

	sp.L_FUNC_FD = getUserSymbolOrDefault(modelHandle, l_fd);

	sp.L_FUNC_TYCH_FD = getUserSymbolOrDefault(modelHandle, l_tych_fd);

	sp.TARGET_FD = getUserSymbolOrDefault(modelHandle, target_fd);

	sp.CONSTR_XU_fd = getUserSymbolOrDefault(modelHandle, constraintsXU_fd);

	sp.CONTROL_ELIGIBILITY_FOR_TRAJ_fd = getUserSymbolOrDefault(modelHandle, controlEligibilityForTraj_fd);

	sp.DYNAMICS_TYCH_FD = getUserSymbolOrDefault(modelHandle, dynamics_tych_fd);
    
	sp.CONSTR_X_fd = getUserSymbolOrDefault(modelHandle, constraintsX_fd);
    
	sp.DYN_CONSTR_FOR_TRAJ = getUserSymbolOrDefault(modelHandle, dynConstraintsForTraj);

	sp.MU_FUNC_FD = getUserSymbolOrDefault(modelHandle, constraintsXUY_fd);

    std::vector<trajectoryParams> &trajParamsList = pm->getTrajectoryParametersList();

    // Il ne peut pas y avoir de CONTROL_WEIGHT par défaut
    controlWeight_t controlWeight = getUserSymbol<controlWeight_t>(modelHandle, "controlWeight");
    neighborValidator_t validator = getUserSymbol<neighborValidator_t>(modelHandle, "isValidNeighbor");
    temporalControl_t temporalControl = getUserSymbol<temporalControl_t>(modelHandle, "temporalControl");
    userTyche_t userTyche = getUserSymbol<userTyche_t>(modelHandle, "tycheValue");
    cumulativeDistribution_t cumulativeDistribution = getUserSymbol<cumulativeDistribution_t>(modelHandle, "cumulativeDistribution");
    probabilityDensity_t probabilityDensity = getUserSymbol<probabilityDensity_t>(modelHandle, "probabilityDensity");

    for (trajectoryParams &tp : trajParamsList) {
        tp.CONTROL_WEIGHT = controlWeight;
        switch (tp.TRAJECTORY_TYPE) {
        case STOCHASTIC:
            tp.SORT_INDEXES = &shuffleControlIndexes;
            break;
        case WEIGHTED_CONTROLS_CAUTIOUS_HEAVY:
        case WEIGHTED_CONTROLS_CAUTIOUS:
        case WEIGHTED_CONTROLS:
            tp.SORT_INDEXES = &sortControlIndexesByWeight;            
            if (tp.CONTROL_WEIGHT == nullptr) {
                spdlog::warn("No control weight function defined even though WEIGHTED_CONTROLS was requested. Defaulting to arbitrary control order");
                tp.SORT_INDEXES = &noSort;
            }
            else {
                spdlog::info("Succesfully loaded symbol: controlWeight");
            }
            break;
        case STRATEGY_LIST:
            for (unsigned long long int i = 0; i < cp.DIM_TY; ++i) {
                tp.TYCHE_PARAMS[i].USER_TYCHE = userTyche;
                tp.TYCHE_PARAMS[i].PROBABILITY_DENSITY = probabilityDensity;
                tp.TYCHE_PARAMS[i].CUMULATIVE_DISTRIBUTION = cumulativeDistribution;
            }
            break;
        default:
            tp.SORT_INDEXES = &noSort;
        }

        tp.IS_VALID_NEIGHBOR = validator;
        switch (tp.BUBBLE_INTERPRETATION) {
        case MOORE:
        case MOORE_PX:
            tp.IS_VALID_NEIGHBOR = &belowInfiniteDistance;
            break;
        case EUCLIDEAN:
        case EUCLIDEAN_PX:
            tp.IS_VALID_NEIGHBOR = &belowEuclideanDistance;
            break;
        case ELLIPTIC:
        case ELLIPTIC_PX:
            tp.IS_VALID_NEIGHBOR = &belowEllipticDistance;
            break;
        case CUSTOM:
            if (tp.IS_VALID_NEIGHBOR == nullptr) {
                spdlog::warn("Custom bubble requested but no custom isValidNeighbor function given. Defaulting to MOORE");
                tp.BUBBLE_INTERPRETATION = MOORE;
                tp.IS_VALID_NEIGHBOR = &belowInfiniteDistance;
            }
            break;
        case CUSTOM_PX:
            if (tp.IS_VALID_NEIGHBOR == nullptr) {
                spdlog::warn("Custom bubble requested but no custom isValidNeighbor function given. Defaulting to MOORE_PX");
                tp.BUBBLE_INTERPRETATION = MOORE_PX;
                tp.IS_VALID_NEIGHBOR = &belowInfiniteDistance;
            }
            break;
        }
        tp.TEMPORAL_CONTROL = temporalControl;
    }
    
    /*
     * Initialisation des paramètres de systèmes dynamique
     * Ici toutes les valeurs sont par defaut, non utilisés
     */

    return pm;

    }

#endif /* SRC_INITPARAMS_H_ */
