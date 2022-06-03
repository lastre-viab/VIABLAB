/*
 * Viabi.h
 *  *
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
 *  Created on: 10 sept. 2013
 *      Author: ANYA
 */

#ifndef VIABI_H_
#define VIABI_H_
//#include "defs.h"
#include "ParametersManager.h"
#include "SysDyn.h"

using namespace std;



class Viabi {
public:
	Viabi();
	Viabi( ParametersManager * pm);
	virtual ~Viabi();
	virtual void printViabiInfo() =0;
	virtual void initialiseTarget() =0;
	virtual void initialiseConstraints() =0;
	virtual void  ViabilityKernel( bool sortieOK,int nbArret)  =0;
	virtual void  CaptureBasin( )  =0;
	virtual void  GarantedViabilityKernel( bool sortieOK,int nbArret)  =0;
	virtual void computeTrajectories() =0;
	virtual void loadViableSets() =0;
	virtual void saveViableSets() =0;

protected:
	ParametersManager* modelParams;
	SysDyn* dynsys;
	int nbOMPThreads;
};



#endif /* VIABI_H_ */
