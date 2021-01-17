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
 *  \author: A. D�silles, LATSRE
 *  \brief  Fichier contant la d�finition d'un mod�le de viabilit� pour les tests
 *
 *
 */


#ifndef TESTDATA_H_
#define TESTDATA_H_


/***********************************************************
 * Some model specific parameters can be defined here
 *************************************************************/

double R=2.0;      //rayon de la boule pour le calcul
double a=-0.7998;  // Paramètres = u=(a,b)
double b= 0.15517;

/*! \var dim
 *  \brief State  dimension
 */
const int dim=2;
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=0;
/*! \var dicret_type
 *  \brief the  type of discretization scheme used  for the dynamics
 *  EE or 1 = Euler Explicit  scheme
 *  EI or 2 = Euler Implicit
 *  RK2I or 3 = RK2 Implicit (RK2 for -F)
 *  RK2E or 4 = RK2 Explicit (RK2 for F)
 */
const int discret_type=0;

/*!
 * \var dynType  defines the type of dynamics in the model
 *      1 or CC : continuous in time and space
 *      2 or DC : discrete time continuous space
 *      2 or DD : discrete time discrete space
 *      4 or HD : hybrid TODO
 */
const int dynType=DC;


/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */
double STATE_MIN[dim]={-R, -R};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={R, R};

/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int nbPointsState[dim]    =  {401,401};
/*
 * Direction le long de laquelle l'ensemble sera représenté par des suites de bits
 */
unsigned long long int dirTramage =1;
/*
 * Sélection de la méthode de représentation de l'ensemble
 * Ceparamètre détermine quelle classe sera utilisée pour les calculs
 *
 *    - BS = BitSet, représentation par fonction caractéristique , pour
 *                   tout type de calculs
 *    - MM = MicroMacro, représentation par valeurs réelles, pour les calculs
 *          d'ensembles épigraphiques , associés aux système micro-macro
 */
int gridMethod=BS;
/*
 * Sélection de l'ensemble à calculer
 */
int setType=VIAB;
/*
 * Nombre de raffinements successifs dans le calcul
 */
int refine=2;

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0};

/*!
 * \var dbFileName
 * the  name for the  main data base file for the roject
 */

string prefix="JuliaTest-";

 /*
 * Definition of the dynamics  and associated functions and constants
 */

/*
 * Paramètre qui indique si l'ensemble doit être calculé
 * dans le cas contrainte il devra pouvoir être chargé depuis un fichier pour une étude de
 * trajectoires
 */
int computeSet=1;
int saveBoundary=1;

const int nbTrajs=0;
double initPoints[dim*nbTrajs]={};


void dynamics(double * x, double *u, double * image);
void dynamics(double * x, double *u, double * image)
{
  image[0]=     x[0]*x[0] - x[1]*x[1]+a;
  image[1]=   2*x[0]*x[1]+ b;
}


/*      *****************************************
 *  Definition of constraints and target
 *************************************************** */

/*!
 * \brief Function  defining the state   constraints, corresonds  to k(x)
 *
 * This function defines the set K for admissible  states
 * @param x state variable
 * @return  value that caraterise the constraints set
 */


inline double constraintsX( double * x )
{
  double res1=  x[0]*x[0] + x[1]*x[1]-R*R;
  double res=(res1<=0)?1.0:PLUS_INF;
  return res;
}



#include "Julia2D_unused.h"
#endif /* TESTDATA_H_ */
