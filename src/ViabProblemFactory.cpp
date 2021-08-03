/*
 * ViabProblemFabric.cpp
 *
 *  Created on: 1 août 2021
 *      Author: adesi
 */

#include "../include/ViabProblemFactory.h"


ViabProblemFactory::ViabProblemFactory() {
	// TODO Auto-generated constructor stub

}

ViabProblemFactory::ViabProblemFactory(ParametersManager * pm)
{
	problemParameters=pm;
}

ViabProblemFactory::~ViabProblemFactory() {
	// TODO Auto-generated destructor stub
}

Viabi * ViabProblemFactory::constructViabilityProblem( int gridMethod){
	Viabi * vp;
	if(gridMethod==BS)
	{
		vp= new ViabiBitSet(problemParameters);
	}
	else
	{
		vp= new ViabiMicroMacro(problemParameters);
	}
	return vp;
}

