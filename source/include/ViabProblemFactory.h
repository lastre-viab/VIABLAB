/*
 * ViabProblemFabric.h
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
 *  Created on: 1 août 2021
 *      Author: adesi
 */

#ifndef SRC_VIABPROBLEMFABRIC_H_
#define SRC_VIABPROBLEMFABRIC_H_

#include "ParametersManager.h"
#include "Viabi.h"
#include "ViabiBitSet.h"
#include "ViabiMicroMacro.h"
#include "ViabiMicroMacroDiscrete.h"

class ViabProblemFactory
    {
public:
    ViabProblemFactory();
    ViabProblemFactory(ParametersManager *pm);
    Viabi* constructViabilityProblem(int gridMethod);
    virtual ~ViabProblemFactory();

private:
    ParametersManager *problemParameters;

    };

#endif /* SRC_VIABPROBLEMFABRIC_H_ */
