+
/*! \file  testZermelo.h
 *
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
 *  \author: A. DESILLES, LATSRE
 *  \brief  This is a model file for test example : computing Julia sets using viability algorithms
 *
 *
 */


#ifndef JULIA2DDATA_H_
#define JULIA2DDATA_H_

/***********************************************************
 * Some model specific parameters can be defined here
 *************************************************************/

double R=2.0;      // Ball radius for computations
double a=-0.7998;  // Dynamics parameters = u=(a,b)
double b= 0.15517;


void dynamics(double * x, double *u, double * image)
{
  image[0]=     x[0]*x[0] - x[1]*x[1]+a;
  image[1]=   2*x[0]*x[1]+ b;
}

inline double constraintsX( double * x )
{
  double res1=  x[0]*x[0] + x[1]*x[1]-R*R;
  double res=(res1<=0)?1.0:PLUS_INF;
  return res;
}

#endif /* JULIA2DDATA_H_ */
