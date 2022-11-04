/*! \file  Equilibres4D_data.h
 *
 * *
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
 *  \author: A. D�silles, LASTRE
 *  \brief  Fichier contant la d�finition d'un mod�le de viabilit� pour les tests
 *
 *  Ce fichier  contient les d�finitions  de tous les �l�ments qui d�crivent un
 *  un probl�me de viabilit� concret et  qui sont consid�r�s comme
 *  param�tres par le code; Ce fichier repr�sente ainsi une interface d'entr�e primitive pour le code.
 *
 *  Les �l�ments principaux d�finis ici  sont :
 *    -la  dynamique f(t,x,u)
 *    -les contraintes:
 *      -# sur l'�tat k(t,x)
 *      -# sur le controles U(x)
 *    - la cible c(t,x)
 *    -les fonctions d�finissant l'objectif  d'un probl�me d'optimisation
 *    	-# fonction l(t,x,u)
 *    	-# fonction m(t,x,u)
 *    - diff�rentes options d�finissant plus pr�cis�ment la nature du probl�me � r�soudre
 *    - diff�rentes options d�finissance la m�thode num�rique � utiliser
 *
 *    Voir plus loins dans les commentaires  du fichier la d�finition de chaque fonction et de chaque param�tre
 *
 *
 *
 */


#ifndef EQUILIBRES4D_DATA_H_
#define EQUILIBRES4D_DATA_H_
/*! \var dim
 *  \brief State  dimension
 */
const int dim=4;
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=0;

/*!
 * \var dynType  defines the type of dynamics in the model
 *      1 or CC : continuous in time and space
 *      2 or DC : discrete time continuous space
 *      2 or DD : discrete time discrete space
 *      4 or HD : hybrid \todo
 */
 const int dynType=CC;

/*! \var dicret_type
 *  \brief the  type of discretization scheme used  for the dynamics
 *  EE or 1 = Euler Explicit  scheme
 *  EI or 2 = Euler Implicit
 *  RK2I or 3 = RK2 Implicit (RK2 for -F)
 *  RK2E or 4 = RK2 Explicit (RK2 for F)
 *  RK4I or 5 = RK2 Implicit (RK4 for -F)
 *  RK4E or 6 = RK2 Explicit (RK4 for F)
 */
const int discret_type=EE;

double c=0.01;

/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */

double STATE_MIN[dim]={-1., -1., -2.0, -c};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={1.,  1., 2.0, c};

/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int nbPointsState[dim]    =  {101,101,201,5};

unsigned long long int dirTramage =2;

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0,0,0};

int saveProjection=1;
unsigned long long int projection[dim]={0,0,0,1};
int intermediate_savings = 1;
string prefix="equi4D-testOMP-";

/*!
 * \var globalDeltaT
 *  bool�en indique si le pas de temps  doit �re choisi globalement
 *  ou localement pour les algorithmes de viabilit�
 */
bool globalDeltaT=false;
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
 * Definition of the dynamics  and associated functions and constants
 */

void dynamics(double * x, double *u, double * image)
{
  image[0]=  0.0;
  image[1]=  0.0;
  image[2]=  0.0;
  image[3]= ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) );
}

/*!
 * jacobian matrix of the dynamics
 */

inline void jacobian(double *x, double *u , double ** jacob)
{
  jacob[0][0]=0.0;
  jacob[0][1]=0.0;
  jacob[0][2]=0.0;
  jacob[0][3]=0.0;

  jacob[1][0]=0.0;
  jacob[1][1]=0.0;
  jacob[1][2]=0.0;
  jacob[1][3]=0.0;

  jacob[2][0]=0.0;
  jacob[2][1]=0.0;
  jacob[2][2]=0.0;
  jacob[2][3]=0.0;

  jacob[3][0]= sign(x[0]-sin(x[2]) )*( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) )+
             ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( sign(x[0]+sin(x[2]) )  );

  jacob[3][1]=( sign( x[1]-cos(x[2]) ) )* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) )+
              ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( sign( x[1]+cos(x[2]) ) );

  jacob[3][2]=( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* (sign(x[0]+sin(x[2]) )*cos(x[2])- sign( x[1]+cos(x[2]) )*sin(x[2]) )+
              ( -sign(x[0]-sin(x[2]) )*cos(x[2]) + sign( x[1]-cos(x[2]) ) * sin(x[2]))* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) );

  jacob[3][3]=0.0;

}
/*!
 * \var computeLC : indicates which method  to use to copute Lipschitz constant
 * 0= analytic  global value used
 * 1= local calculation using the jacobian matrix
 * 2= local calculation using the finite differences
 */
const int computeLC=1;

/*!
 * \var LC Lipschitz constant if known analytically
 */
double LC= 20.0;

/*!
 * \var M : global bound  for the dynamics if defined
 */
double M= 20.0;
/*!
 * \var computeM : indicates which method  to use to compute the bound for dynamics
 * 0= analytic  global value used
 * 1= local calculation using the localDynBounds function
 * 2= local calculation using  explicit maximization on controls
 */
const int computeM=1;


/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
inline void localDynBounds(double * x, double * res)
{
  res[0]=0.;
  res[1]=0.;
  res[2]=0.;
  res[3]= ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) );

}



/*      *****************************************
 *  Definition of constraints and target
 *************************************************** */

/*!
 * \brief Function  defining the state   constraints, corresponds  to k(x)
 *
 * This function defines the set K for admissible  states
 * @param x state variable
 * @return  value that characterize the constraints set
 */

inline double constraintsX( double * x )
{
  return 1.0;
}

#include "equilibres4D_unused.h"

#endif /* TESTDATA_H_ */
