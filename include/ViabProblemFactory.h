/*
 * ViabProblemFabric.h
 *
 *  Created on: 1 août 2021
 *      Author: adesi
 */

#ifndef SRC_VIABPROBLEMFABRIC_H_
#define SRC_VIABPROBLEMFABRIC_H_

#include "ParametersManager.h"
#include "Viabi.h"
#include "ViabiBitSet.h"
#include "ViabiMicroMacro.h"


class ViabProblemFactory {
public:
	ViabProblemFactory();
	ViabProblemFactory(ParametersManager * pm);
	Viabi * constructViabilityProblem( int gridMethod);
	virtual ~ViabProblemFactory();

private :
	ParametersManager *problemParameters;

};

#endif /* SRC_VIABPROBLEMFABRIC_H_ */
