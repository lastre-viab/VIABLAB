/*! \file  testZermelo.h
 *
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
 *  \author: A. D�silles, LATSRE
 *  \brief  Fichier contant la d�finition d'un mod�le de viabilit� pour les tests
 *
 *
 *
 *
 */

#include "../include/utilities.h"

extern "C" {

//	    A = autre cas de la forme f(x,u,v)
//		    ou u=(u[1],u[2],u[3],u[4],u[5]) in [alphamin,alphamax]x[betamin,betamax]x[gammamin,gammamax] }

/*
 * Some model specific parameters can be defined here
 */
double b = 0.4;
double  r = 1.0;
double Pmax = 0.5;
double Lmin = 0.1;
double Lmax = 1.0;
double c = 0.1;

void dynamics(const double *x, const double *u, double *image)
{

	image[0]=  c*u[0]; //k'
	double x1p8 = x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1];
	image[1]=  -b*x[1] + x[0]+ r * x1p8 /(1+x1p8); //m'
	//cout<< " dynamique renvoie "<<image[0]<< " "<<image[1]<<endl;
}

/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
void localDynBounds(const double *x, double * res)
{

	res[0]=  c;
	double x1p8 = x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1];

	res[1] = abs(-b*x[1] + x[0]+ r * x1p8 /(1+x1p8));
}

void jacobian(const double *x, const double *u , double ** jacob)
{

	jacob[0][0]=0.0;
	jacob[0][1]=0.0;
	double x1p7 = x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1];

	double x1p8 = x[1]*x1p7;

	jacob[1][0]= 1.0;
	jacob[1][1]=-b + r*8*x1p7/((1+x1p8)*(1+x1p8) );

}

 double l(const double *x, const double *u )
{
   double res = PLUS_INF;
   if(x[0]>=Lmin)
   {
	   res = ( (x[0] <= Lmax) && (x[1] <= Pmax)) ? 0.0 : 1.0;
   }
   return res;
}

 double constraintsX( const double *x )
{
 double kSign = max( (x[0] - Lmin)*(x[0] - Lmax), (x[1] - Pmax));
 return  kSign <= 0 ? 0.0 : (x[0] >= Lmin)? 0: PLUS_INF;
	//return 0.0;
}

}
