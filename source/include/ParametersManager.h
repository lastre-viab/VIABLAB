/*
 * ParametersManager.h
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
 *  Created on: 27 juil. 2021
 *      Author: Anna DESILLES
 */

#ifndef SRC_PARAMETERSMANAGER_H_
#define SRC_PARAMETERSMANAGER_H_

#include "../include/defs.h"

class ParametersManager
    {
public:
    ParametersManager();
    ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams *cp,
	    systemParams *sp);
    ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams *cp,
	    systemParams *sp, int nbOmpThreads, string paramsFile);
    gridParams* getGridParameters();
    algoViabiParams* getAlgoParameters();
    controlParams* getControlParameters();
    systemParams* getSystemParameters();
    void readControlParametersFromJson();
    void readGridParametersFromJson();
    void readSystemParametersFromJson();
    void readAlgoParametersFromJson();
    virtual ~ParametersManager();

private:
    gridParams *gridParameters;
    algoViabiParams *algoParameters;
    controlParams *controlParameters;
    systemParams *systemParameters;
    int nbOmpThreads;
    string parametersFileName;
    void readTabData(ptree *dataRoot, double *target, string label,
	    int nbElements);
    void readTabData(ptree *dataRoot, unsigned long long int *target,
	    string label, int nbElements);
    void readTabData(ptree *dataRoot, int *target, string label,
	    int nbElements);
    };

#endif /* SRC_PARAMETERSMANAGER_H_ */
