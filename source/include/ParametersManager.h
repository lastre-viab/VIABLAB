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

#include "Params.h"
#include "ModelParameters.h"
#include "viablab_export.h"

using namespace boost::property_tree;

class VIABLAB_LIBRARY_EXPORT ParametersManager
{       
public:
    ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams *cp, systemParams *sp, int nbOmpThreads, const string &paramsFile, void *modelHandle);
    const gridParams* getGridParameters() const;
    const algoViabiParams* getAlgoParameters() const;
    const controlParams* getControlParameters() const;
    const systemParams* getSystemParameters() const;
    const std::vector<trajectoryParams> &getTrajectoryParametersList() const;
    std::vector<trajectoryParams> &getTrajectoryParametersList();
    int getNbTrajectories() const;
    
    void *getModelHandle();
    const modelParams *getModelParameters() const;
    
    void readControlParametersFromJson(ptree &allParamsRoot);
    void readGridParametersFromJson(ptree &allParamsRoot);
    void readSystemParametersFromJson(ptree &allParamsRoot);
    void readAlgoParametersFromJson(ptree &allParamsRoot);
    void readTrajectoryParametersListFromJson(ptree &allParamsRoot);
    void readModelParametersFromJson(ptree &allParamsRoot);
    virtual ~ParametersManager();

private:
    gridParams *gridParameters;
    algoViabiParams *algoParameters;
    controlParams *controlParameters;
    systemParams *systemParameters;
    std::vector<trajectoryParams> trajectoryParametersList;
    
    void *modelHandle;
    int nbOmpThreads;
    string parametersFileName;

    modelParams modelParameters;

    template<typename T>
    void readTabDataSkipInvalid(ptree *dataRoot, T *target, const string &label, int &nbElements);

    template<typename T>
    void readTabData(ptree *dataRoot, T *target, const string &label, int nbElements, const T &defaultValue);
        

    void readTrajectoryParameters(ptree &trajParamsRoot, trajectoryParams &params);
    void readTycheParameters(ptree &tycheParamsRoot, tycheParams &params);
    
    void toPixels(double *bubbleRadius);
    void initBubble(ptree &dataRoot, trajectoryParams &params);
    void initSeed(ptree &dataRoot, trajectoryParams &params);
};

class TrajectoryParametersManager {
public:
    TrajectoryParametersManager(ParametersManager *pm, int trajIndex);
    const trajectoryParams *getTrajectoryParameters() const;
    int getNbTrajectories() const;
    int getTrajectoryIndex() const;
    const gridParams* getGridParameters() const;    
    const algoViabiParams* getAlgoParameters() const;    
    const controlParams* getControlParameters() const;    
    const systemParams* getSystemParameters() const;        
    void *getModelHandle(void);    
    const modelParams *getModelParameters() const;   
private:
    ParametersManager *parametersManager;
    int trajIndex;
};

void mergeJSONPtreeInto(const ptree &from, ptree &to);

#endif /* SRC_PARAMETERSMANAGER_H_ */
