/*
 * ParametersManager.cpp
 *
 *  Created on: 27 juil. 2021
 *      Author: adesi
 */

#include "../include/ParametersManager.h"


ParametersManager::ParametersManager() {
	// TODO Auto-generated constructor stub

}

ParametersManager::ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams *cp,systemParams *sp)
{
	gridParameters =gp;
	algoParameters=avp;
	controlParameters=cp;
	systemParameters=sp;
}

gridParams * ParametersManager::getGridParameters(){
	return gridParameters;
}
algoViabiParams * ParametersManager::getAlgoParameters(){
	return algoParameters;
}
controlParams * ParametersManager::getControlParameters(){
	return controlParameters;
}
systemParams * ParametersManager::getSystemParameters(){
	return systemParameters;
}
ParametersManager::~ParametersManager() {
	// TODO Auto-generated destructor stub
}

