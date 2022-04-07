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
	ParametersManager(gridParams *gp, algoViabiParams *avp,  controlParams * cp, string cpFile,systemParams *sp);
		gridParams * getGridParameters();
	algoViabiParams * getAlgoParameters();
	controlParams * getControlParameters();
	systemParams * getSystemParameters();
	void readControlParametersFromJson( string file);

	virtual ~ParametersManager();

private:
	gridParams *gridParameters;
	algoViabiParams *algoParameters;
	controlParams *controlParameters;
	systemParams *systemParameters;
	void readTabData (ptree *dataRoot, double * target, string label, int nbElements );
	void readTabData (ptree *dataRoot, unsigned long long int * target, string label, int nbElements );
};

#endif /* SRC_PARAMETERSMANAGER_H_ */
