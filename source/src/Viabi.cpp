/*
 * Viabi.cpp
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
 *  Created on: 7 juil. 2013
 *      Author: Anna DESILLES
 */

#include "../include/Viabi.h"

Viabi::Viabi(ParametersManager *pm)
    {
    modelParams = pm;
    nbOMPThreads = 1;
    }

Viabi::Viabi()
    {
    // TODO Auto-generated constructor stub
//cout<< "  coucou constructeur par defaut\n";
    }

Viabi::~Viabi()
    {
    // TODO Auto-generated destructor stub
    }
