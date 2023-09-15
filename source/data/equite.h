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
double alpha=0.75;
double beta = 0.25;

void dynamics(double * x, double *u, double * image)
{
	double p = x[0], q= x[1], r = x[2], sig = u[0], v = u[1];
	double k = -log(1 - p), s = -log(1 - q), c = -log(1 - r);
	double kprime = pow(k, alpha)*pow( sig * s, beta) - c;
	double sprime = -sig * s;
	double cprime = v * c;
	if(p >=1.0)
	{
		image[0]=   0.0; //k'
	}
	else
	{
		image[0]=  (1.0 - p) * kprime; //k'
	}
	if(q >=1.0)
	{
		image[1]=   0.0; //k'
	}
	else
	{
		image[1]=  (1-q) * sprime;
	}
	//m'
	if(r >=1.0)
	{
		image[2]=   0.0; //k'
	}
	else
	{
		image[2]=   (1-r) * cprime;
	}
	//cout<< " dynamique renvoie "<<image[0]<< " "<<image[1]<<endl;
}

/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */

inline void localDynBounds(double * x, double * res)
{
	double p = x[0], q= x[1], r = x[2];

	if(p>=1.0)
	{
		res[0] = 0.0;
	}
	else
	{
		double k = -log(1 - p), s = -log(1 - q), c = -log(1 - r);
				double kprime1 = abs( pow(k, alpha)*pow( s, beta) - c);
				double kprime0 = abs( -c);
		res[0]=  (1.0 - p) * max(kprime0, kprime1);
	}

	res[1]=  abs((1-q) * log( 1 - q));
	res[2] = abs((1-r) * log( 1 - r));
}

inline void jacobian(double *x, double *u , double ** jacob)
{
	double p = x[0], q= x[1], r = x[2], sig = u[0], v = u[1];
	if(p>=1.0)
	{
		jacob[0][0]= 0.0;
		jacob[0][1]=  0.0;
		jacob[0][2]= 0.0;
	}
	else
	{
		double fact = pow(-log(1-p), alpha) * pow(-sig * log(1 - q), beta) + log(1 - r);

		jacob[0][0]= - fact + alpha *  pow(sig * log(1-q) /  log(1.0 - p), beta);
		jacob[0][1]=  sig * beta * (( 1.0 - p) / (1-q) ) * pow(sig * log(1-p) /  log(1.0 - q), alpha);
		jacob[0][2]= (1.0 - p) / (1.0 - r );
	}

	jacob[1][0]=0.0;
	jacob[1][1]=-u[0]*(1.0 + log(1.0 - q));
	jacob[1][2]=0.0;

	jacob[2][0]=0.0;
	jacob[2][1]=0.0;
	jacob[2][2]= u[1]*(1.0 + log(1.0 - r));

}

inline double l(double * x, double * u )
{
	return 0.0;
}

inline double constraintsX( double * x )
{
	return 1.0;
}

#endif /* TESTDATA_H_ */
