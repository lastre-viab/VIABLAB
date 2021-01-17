/*! \file  Equilibres4D_data.h
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
#ifndef THIBAUT_DATA_H_
#define THIBAUT_DATA_H_

/*! \var dim
 *  \brief State  dimension
 */
const int dim=3;
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=0;


/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MIN[dimc]={};
/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX[dimc]={};

/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {};

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
const int discret_type=RK4E;

// Work area
//double localXMin,localXMax,localYMin,localYMax,localZMin,localZMax;
double localXMin = -1;
double localXMax = 1;
double localYMin = -1;
double localYMax = 1;
double localZMin = -1;
double localZMax = 1;

// Real K
//double localXMinK,localXMaxK,localYMinK,localYMaxK,localZMinK,localZMaxK;
double localXMinK = localXMin;
double localXMaxK = localXMax;
double localYMinK = localYMin;
double localYMaxK = localYMax;
double localZMinK = localZMin;
double localZMaxK = localZMax;

/*
 * Pavé de contraintes
 */
double STATE_K_MIN[dim]={localXMinK, localYMinK, localZMinK};

double STATE_K_MAX[dim]={localXMaxK,  localYMaxK, localZMaxK};


int res=100;
int resX=res;
int resY=res;
int resZ=res;

/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */

//double *CONTROL_MIN, *CONTROL_MAX;
double STATE_MIN[dim]={localXMin, localYMin, localZMin};
/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={localXMax,  localYMax, localZMax};
/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int nbPointsState[dim] = {resX,resY,10};

unsigned long long int dirTramage =2;
/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0,0};

/*
 * sélection de type de reconstruction de trajectoire
 *  VD= 1, viable par défaut, on sélectionne le premier contrôle viable trouvé
 *  VL= 2, viable lourd: le contrôle reste constant tant qu'il viable;
 *  cette méthode nécessite une initialisation de contrôle
 */
int typeTraj=VL;


string prefix="ThibautCircleRes"+to_string(res)+"NumericalScheme"+to_string(discret_type);

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
  image[0]= -x[1];
  image[1] = x[0];
  image[2]=0;
}

/*!
 * jacobian matrix of the dynamics
 */

inline void jacobian(double *x, double *u , double ** jacob)
{
  jacob[0][0] = 0;
  jacob[0][1] = -1;
  jacob[0][2] = 0;

  jacob[1][0] = 1;
  jacob[1][1] = 0;
  jacob[1][2] = 0;

  jacob[2][0] = 0;
  jacob[2][1] = 0;
  jacob[2][2] = 0;

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
double LC= 200;

/*!
 * \var M : global bound  for the dynamics if defined
 */
double M= 200;
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

  res[0]=abs(x[1]) ;
  res[1]= abs(x[0]);
  res[2]=0.0;


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
  /* bool b1min,b2min,b3min, b1max,b2max,b3max;
  b1min = x[0] >= localXMinK;
  b2min = x[1] >= localYMinK;
  b3min = x[2] >= localZMinK;

  b1max = x[0] <= localXMaxK;
  b2max = x[1] <= localYMaxK;
  b3max = x[2] <= localZMaxK;

  if (b1min && b2min && b3min && b1max && b2max && b3max)
    return 1.0;
  else
    return 0;*/
/*
 * Pour les contraintes de type pavé on calcule max {i=1:dim}   max(STATE_K_MIN[i]-x[i], x[i]-STATE_K_MAX[i])
 * Si ce max  est < 0  c'es équivalent à :
 *  pour tout i  STATE_K_MIN[i] < x[i]< STATE_K_MAX[i]
 *  et donc les contraintes sont vérifiées
 */
  double res=-PLUS_INF;
  for(int  i=0;i<dim;i++){
    res = max(res, max(STATE_K_MIN[i]-x[i], x[i]-STATE_K_MAX[i]));
  }
  double resFinal=(res>0)?PLUS_INF:0.0;
  return resFinal;
}

#include "Cylinder_unused.h"

#endif /* TESTDATA_H_ */
