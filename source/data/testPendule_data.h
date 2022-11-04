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


#ifndef TESTDATA_H_
#define TESTDATA_H_





//	    A = autre cas de la forme f(x,u,v)
//		    ou u=(u[1],u[2],u[3],u[4],u[5]) in [alphamin,alphamax]x[betamin,betamax]x[gammamin,gammamax] }

/*
 * Some model specific parameters can be defined here
 */
double a = 0.01;
double b = 0.0;
double c = 9.81;



void loadModelData()
{


}

void dynamics(double * x, double *u, double * image)
{
	image[0] = x[1];
	image[1]=  -a*x[1] - c*sin(x[0])+b +c*u[0]; //m'
}

/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
inline void localDynBounds(double * x, double * res)
{
	res[0] = abs(x[1]);
	res[1]= abs(-a*x[1] - c*sin(x[0])+b) +c;
}

inline void jacobian(double *x, double *u , double ** jacob)
{

	jacob[0][0]=0.0;
	jacob[0][1]=1.0;

	jacob[1][0]=-c*cos(x[0]);
	jacob[1][1]=-a;
}

double eps = 0.1;
inline  double target (double * x)
{
	double res;
	if(sqrt((x[0])*(x[0])+x[1]*x[1])<eps)
	{
		res=0.0;
	}
	else
	{
		res=PLUS_INF;
	}
	return res;
}
#endif /* TESTDATA_H_ */
