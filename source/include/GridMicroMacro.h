/*
 * GridMicroMacro.h
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
 *  Created on: 26 nov. 2013
 *      Author: ANYA
 */

#ifndef GRIDHJB_H_
#define GRIDHJB_H_

#include "Grid.h"
#include "SysDyn.h"
#include "Params.h"

#include <iostream>

using namespace std;
class GridMicroMacro final: public Grid
    {
public:
    GridMicroMacro();
    GridMicroMacro(const gridParams &gp);
    virtual ~GridMicroMacro();
    virtual void printGrid(void) const;
    void addPointToSet(const unsigned long long int *coords, double value);
    void addPointToSet(unsigned long long int pos, double value);

    void loadSet(const string &fileName);
    virtual bool isInSet(const unsigned long long int *coords) const;
    virtual unsigned long long int getNearestPointInSet(const double *coords) const;
    unsigned long long int getBestNearPointInSet(const double *coords) const;
    double getOptimalValue(const double *coords) const;
    virtual void savePointsList(const string &fileName) const;
    virtual void saveValOnGrid(const string &fileName) const;
    virtual void saveValOnGridLight(const string &fileName) const;
    void saveValOnGrid_DD(const string &fileName) const;

    void saveProjection(const string &fileName, const unsigned long long int *projection) const;
    void computeMinMaxValues(double &minV, double &maxV) const;
    const double* getGridPtr() const;
    double* getGridPtr();
    const double* getGridPtr_tmp() const;
    void copyGrid(const double *grIn, double *grOut) const;

    void saveSubLevelset(double level, const string &fileName) const;
    void saveSubLevelset_DD(double level, const string &fileName) const;
    void saveCoupeBoundary(const string &nomFichier) const;
    void saveCoupeBoundary_DD(const string &nomFichier) const;
    int nbOMPThreads;
private:

    double *gridPtr;
    double *gridPtr_tmp;

    };

#endif /* GRIDHJB_H_ */
