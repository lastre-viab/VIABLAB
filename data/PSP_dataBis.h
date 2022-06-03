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



string paramsFile = "testPSPBis.json";

//	    A = autre cas de la forme f(x,u,v)
//		    ou u=(u[1],u[2],u[3],u[4],u[5]) in [alphamin,alphamax]x[betamin,betamax]x[gammamin,gammamax] }

/*
 * Some model specific parameters can be defined here
 */
double alpha=0.75;
double beta = 0.25;

double c = 1.0;
double r = 5.0;

void loadModelData()
{


}

void dynamics(double * x, double *u, double * image)
{
	if(u[0]<=0.0)
	{
		image[0]=   - x[2]; //k'
	}
	else
	{
		if(x[0]<=0.0)
		{
			image[0]=   - x[2]; //k'
		}
		else
		{
			image[0]=  pow(x[0], alpha)*pow(r*u[0], beta) - x[2]; //k'
		}
	}
	image[1]=  -r*u[0]; //m'
	image[2] = c*u[1]; //c'
	//cout<< " dynamique renvoie "<<image[0]<< " "<<image[1]<<endl;
}

/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
inline void localDynBounds(double * x, double * res)
{
	if(x[0]<=0.0)
	{
		res[0] = x[2];
	}
	else
	{
		res[0]= max(x[2], abs( pow(r,beta)*pow(x[0], alpha) - x[2])); //k'
	}
	res[1]=  r;
	res[2] = c;
}

inline void jacobian(double *x, double *u , double ** jacob)
{

	if(u[0]<=0.0)
		{
			jacob[0][0]=0.0;
				jacob[0][1]=0.0;
				jacob[0][2]=-1.0;
		}
		else
		{
			if(x[0]<=0.0)
			{
				jacob[0][0]=0.0;
				jacob[0][1]=0.0;
				jacob[0][2]=-1.0;
			}
			else
			{
				jacob[0][0]=alpha/pow(x[0], 1.0-alpha)*pow(r*u[0], beta);
				jacob[0][1]=0.0;
				jacob[0][2]=-1.0;
			}
		}

	jacob[1][0]=0.0;
	jacob[1][1]=0.0;
	jacob[1][2]=0.0;

	jacob[2][0]=0.0;
	jacob[2][1]=0.0;
	jacob[2][2]=0.0;

}

inline double l(double * x, double * u )
{
  return 0.0;
}

inline double constraintsX( double * x )
{
  return -x[2];
}

#endif /* TESTDATA_H_ */
