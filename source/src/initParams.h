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

#include "../include/defs.h"
#include "../include/ParametersManager.h"

ParametersManager* initParams(gridParams &gp, algoViabiParams &avp,
	controlParams &cp, systemParams &sp, int nbOmpThreads);
ParametersManager* initParams(gridParams &gp, algoViabiParams &avp,
	controlParams &cp, systemParams &sp, int nbOmpThreads)
    {

    ParametersManager *pm = new ParametersManager(&gp, &avp, &cp, &sp,
	    nbOmpThreads, paramsFile);

    if (dynamics)
	{
	sp.DYNAMICS = &dynamics;
	}
    else
	{
	sp.DYNAMICS = &dynamics_default;
	}

    if (dynamics_tych)
    	{
    	sp.DYNAMICS_TYCH = &dynamics_tych;
    	}
        else
    	{
    	sp.DYNAMICS_TYCH = &dynamics_tych_default;
    	}
    if (sp.DYN_TYPE == 4)
	{
	sp.DYNAMICS = &dynamics_hybrid;
	}

    switch (sp.DYN_TYPE)
	{
    case 2:
	{
	sp.COMPUTE_LC = 0.0;
	sp.COMPUTE_MF = 0.0;

	sp.LIP = 1.0;
	sp.MF = 1.0;
	break;
	}
    case 3:
	{
	sp.COMPUTE_LC = 0.0;
	sp.COMPUTE_MF = 0.0;
	sp.LIP = 1.0;
	sp.MF = 1.0;
	break;
	}
	}
    if (constraintsX)
	{
	sp.CONSTR_X = &constraintsX;
	}
    else
	{
	sp.CONSTR_X = &constraintsX_default;
	}

    if (constraintsXU)
	{
	sp.CONSTR_XU = &constraintsXU;
	}
    else
	{
	sp.CONSTR_XU = &constraintsXU_default;
	}
    if (jacobian)
	{
	sp.JACOBIAN = &jacobian;
	}
    else
	{
	sp.JACOBIAN = &jacobian_default;
	}
    if (jacobian_tych)
    	{
    	sp.JACOBIAN_TYCH = &jacobian_tych;
    	}
        else
    	{
    	sp.JACOBIAN_TYCH = &jacobian_tych_default;
    	}
    if (localDynBounds)
	{
	sp.LOCAL_DYN_BOUNDS = &localDynBounds;
	}
    else
	{
	sp.LOCAL_DYN_BOUNDS = &localDynBounds_default;
	}
    if (m)
	{
	sp.M_FUNC = &m;
	}
    else
	{
	sp.M_FUNC = &m_default;
	}
    if (l)
	{
	sp.L_FUNC = &l;
	}
    else
	{
	sp.L_FUNC = &l_default;
	}

    if (m_tych)
    	{
    	sp.M_FUNC_TYCH = &m_tych;
    	}
        else
    	{
    	sp.M_FUNC_TYCH = &m_tych_default;
    	}
        if (l_tych)
    	{
    	sp.L_FUNC_TYCH = &l_tych;
    	}
        else
    	{
    	sp.L_FUNC_TYCH = &l_tych_default;
    	}

    if (target)
	{
	sp.TARGET = &target;
	}
    else
	{
	sp.TARGET = &target_default;
	}

    if (l_fd)
	{
	sp.L_FUNC_FD = &l_fd;
	}
    else
	{
	sp.L_FUNC_FD = &l_fd_default;
	}

    if (l_tych_fd)
	{
	sp.L_FUNC_TYCH_FD = &l_tych_fd;
	}
    else
	{
	sp.L_FUNC_TYCH_FD = &l_fd_tych_default;
	}

    if (target_fd)
	{
	sp.TARGET_FD = &target_fd;
	}
    else
	{
	sp.TARGET_FD = &target_fd_default;
	}

    if (constraintsXU_fd)
	sp.CONSTR_XU_fd = &constraintsXU_fd;
    else
	sp.CONSTR_XU_fd = &constraintsXU_fd_default;

    if (controlEligibilityForTraj_fd)
	sp.CONTROL_ELIGIBILITY_FOR_TRAJ_fd = &controlEligibilityForTraj_fd;
    else
	sp.CONTROL_ELIGIBILITY_FOR_TRAJ_fd =
		&controlEligibilityForTraj_fd_default;

    if (dynamics_fd)
	sp.DYNAMICS_FD = &dynamics_fd;
    else
	sp.DYNAMICS_FD = NULL;
    if (dynamics_tych_fd)
	sp.DYNAMICS_TYCH_FD = &dynamics_tych_fd;
    else
	sp.DYNAMICS_TYCH_FD = NULL;
    if (constraintsX_fd)
	sp.CONSTR_X_fd = &constraintsX_fd;
    else
	sp.CONSTR_X_fd = NULL;
    if (dynConstraintsForTraj)
	sp.DYN_CONSTR_FOR_TRAJ = &dynConstraintsForTraj;
    else
	sp.DYN_CONSTR_FOR_TRAJ = &dynConstraintsForTraj_default;

    if (constraintsXUY_fd)
	sp.MU_FUNC_FD = &constraintsXUY_fd;
    else
	sp.MU_FUNC_FD = &constraintsXUY_fd_default;
    /*
     * Initialisation des paramètres de systèmes dynamique
     * Ici toutes les valeurs sont par defaut, non utilisés
     */

    return pm;

    }

#endif /* SRC_INITPARAMS_H_ */
