/*
 * ViabiMicroMacroTrajectoryHelper.h
 *
 * VIABLAB : a numerical library for Mathematical Viability Computations
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
 *  Created on: 7 juil. 2024
 *      Author: adesi
 */

#ifndef SRC_VIABIMICROMACROTRAJECTORYHELPER_H_
#define SRC_VIABIMICROMACROTRAJECTORYHELPER_H_

#include "GridMicroMacro.h"
#include "SysDyn.h"
#include "ParametersManager.h"
#include "ViabiTrajectoryHelper.h"

class ViabiMicroMacroTrajectoryHelper final :public ViabiTrajectoryHelper
{
public:
    ViabiMicroMacroTrajectoryHelper() = default;

    ViabiMicroMacroTrajectoryHelper(GridMicroMacro *gr, SysDyn *ds, TrajectoryParametersManager *params);

    virtual ~ViabiMicroMacroTrajectoryHelper();

    double computeOptimalTrajectory(double *initPosition, double timeHorizon, string fileName, bool &succes);

    int findOptiControl(double *currentPos, double budget, unsigned long long int optimDiscreteSuccessor, double &dt, int nbStepIter,
	    double stepCoeff, double *resPos, bool &succes, double &nextBudget);

    unsigned long long int findOptimalDiscreteSuccessor(unsigned long long int pos, double normalizedTime, double dt);

    double computeOptimalTrajectory_Lmin(double *initPosition, string fileName, bool &succes);

    int findOptiControl_Lmin(double budget, double *currentPos, double &dt, int nbStepIter, double stepCoeff, double *resPos, double &newBudget,
	    bool &succes);

    double computeViableTrajectory(double *initPosition, double initValue, string fileName, bool &succes);
    double computeViableTrajectory_DD(unsigned long long int *initPosition, double initValue, string fileName, bool &succes);
    double computeViableTrajectory_tych_DD(unsigned long long int *initPosition, double initValue, string fileName, bool &succes);

    unsigned long long int (ViabiMicroMacroTrajectoryHelper::*findfViableControl_DD)(double budget, unsigned long long int *currentPos,
	    unsigned long long int currentControl, unsigned long long int *resPos, double &newBudget, bool &succes);
    unsigned long long int findViabControlDiffControl_DD(double budget, unsigned long long int *currentPos, unsigned long long int currentControl,
	    unsigned long long int *resPos, double &newBudget, bool &succes);
    unsigned long long int findViabControlDefault_DD(double budget, unsigned long long int *currentPos, unsigned long long int currentControl,
	    unsigned long long int *resPos, double &newBudget, bool &succes);
    unsigned long long int findViabControlDefault_tych_DD(double budget, unsigned long long int *currentPos, unsigned long long int currentControl,
	    unsigned long long int currentTych, unsigned long long int *resPos, double &newBudget, bool &succes);
    unsigned long long int findViabControlMinValue_DD(double budget, unsigned long long int *currentPos, unsigned long long int currentControl,
	    unsigned long long int *resPos, double &newBudget, bool &succes);

    int *getPreferedControlIndexes();
    int *sortPreferedControlIndexes(const double *x, double t, int strategyIndex = 0);
    virtual PointStatus checkKernelRelation(double *position) const;
    
private:    
    GridMicroMacro *grid;

    /*!
     *  \brief Pointeur sur la base de données servant à enregister la rétroaction optimale
     */
    double *vTab;
    int typeTraj;
    string filePrefix;
};

#endif /* SRC_VIABIMICROMACROTRAJECTORYHELPER_H_ */
