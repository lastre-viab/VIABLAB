/*
 * GridHJB.h
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

#include <iostream>




using namespace std;
class GridMicroMacro: public Grid {
public:
  GridMicroMacro();
  GridMicroMacro(  gridParams gp);
  virtual ~GridMicroMacro();
  virtual void printGrid(void);
  void addPointToSet(unsigned long long int * coords, double value)  ;
  void addPointToSet(unsigned long long int pos, double value)  ;

  void loadSet(string fileName);
  virtual bool isInSet(unsigned long long int * coords );
  virtual void savePointsList(string fileName);
  virtual void saveValOnGrid(string fileName);
  void saveValOnGrid_DD(string fileName);

  void printFoncValVox_fd(string fileName, double scaleVal, int epi) ;
  void  addPointToSet_withData(unsigned long long int *pos, double val, void * dataBuff, int dataLength) ;
   void saveProjetion(string fileName, unsigned long long int * projection);
  void computeMinMaxValues(double &minV, double & maxV);
  double * getGridPtr();
  double * getGridPtr_tmp();
  void copyGrid( double *grIn, double *grOut);

  void saveSubLevelset(double level, string fileName );
  void saveSubLevelset_DD(double level, string fileName );
  void saveCoupeBoundary(string nomFichier);
  int nbOMPThreads;
private:


  double *gridPtr;
  double *gridPtr_tmp;



};

#endif /* GRIDHJB_H_ */
