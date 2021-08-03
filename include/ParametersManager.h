/*
 * ParametersManager.h
 *
 *  Created on: 27 juil. 2021
 *      Author: adesi
 */

#ifndef SRC_PARAMETERSMANAGER_H_
#define SRC_PARAMETERSMANAGER_H_

#include "../include/defs.h"

class ParametersManager {
public:
	ParametersManager();
	ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams *cp,systemParams *sp);
	gridParams * getGridParameters();
	algoViabiParams * getAlgoParameters();
	controlParams * getControlParameters();
	systemParams * getSystemParameters();
	virtual ~ParametersManager();

private:
	gridParams *gridParameters;
	algoViabiParams *algoParameters;
	controlParams *controlParameters;
	systemParams *systemParameters;
};

#endif /* SRC_PARAMETERSMANAGER_H_ */
