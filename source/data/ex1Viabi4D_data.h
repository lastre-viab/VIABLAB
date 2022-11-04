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


#ifndef TESTDATA_H_
#define TESTDATA_H_

#include "ex1Viabi4D_unused.h"


void loadModelData();

/*! \var dim
 *  \brief State  dimension
 */
const int dim=4;
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=1;
/*! \var dicret_type
 *  \brief the  type of discretization scheme used  for the dynamics
 *  EE or 1 = Euler Explicit  scheme
 *  EI or 2 = Euler Implicit
 *  RK2I or 3 = RK2 Implicit (RK2 for -F)
 *  RK2E or 4 = RK2 Explicit (RK2 for F)
 */
const int discret_type=RK4E;

/*!
 * \var T maximum time horizon for the study
 */
double T=4.5;

/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MIN[dimc]={-0.9};
/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX[dimc]={0.9};

/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */
double STATE_MIN[dim]={0.0,0.0,0.0, 0.0};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={2.0,2.0,2.0, 2.5};

int sortieOK[dim]={0,0,0,0};


/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {21};
/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int nbPointsState[dim]    =  {41,41,41,41};
unsigned long long int dirTramage =1;
int refine=3;
int saveCoupeBound=1;
int saveCoupe=0;
int sliceDirs[dim]={0, 0, 1, 0};
double sliceVals[dim]={0., 0, 1.0, 0.0};

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0,0,0};
unsigned long long int projection[dim]={0,0,0,0};

unsigned long long int trajProjection[dim]={0,0,0,0};
double initPoint[dim]={0.5,0,0,0};

int scaling[dim]={0,0,0,0};


/*
 * target = 1
 * departure =0;
 * Ce parametre determine le sens des trajectoires
 */
int target_or_departure_problem=1;



/*!
 * \var dbFileName
 * the  name for the  main data base file for the roject
 */

string prefix="ex1viab4D--";

/*!
 * \var globalDeltaT
 *  bool�en indique si le pas de temps  doit �re choisi globalement
 *  ou localement pour les algorithmes de viabilit�
 */
bool globalDeltaT=false;

//	    A = autre cas de la forme f(x,u,v)
//		    ou u=(u[1],u[2],u[3],u[4],u[5]) in [alphamin,alphamax]x[betamin,betamax]x[gammamin,gammamax] }

/*
 * Some model specific parameters can be defined here
 */
double alpha=0.2;

/*
 * Definition of the dynamics  and associated functions and constants
 */


int compteOnlyOptimalRetro=1;

int computeSet=1;


void loadModelData()
{

}


void dynamics_continuous(double * x, double *u, double * image);
void dynamics_discrete(double * x, double *u, double * image);


void dynamics_continuous(double * x, double *u, double * image)
{

  double m=0.0;
  for(int i=1;i<dim-1;i++)
    {
      m+=x[i];
    }
  m=alpha*m*m/((dim-2)*(dim-2));
  image[0]=   x[0]+m-x[3];
  image[1]=  0.0;
  image[2]=  0.0;
  image[3]=  u[0];
  //cout<< " dynamique renvoie "<<image[0]<< " "<<image[1]<<endl;
}
void dynamics_discrete(double * x, double *u, double * image)
{
  image[0]=  0.0;
  image[1]=  0.0;
}


/*!
 * jacobian matrix of the dynamics
 */
inline void jacobian(double *x, double *u , double ** jacob)
{
double m=0.0;
  for(int i=1;i<dim-1;i++)
    {
      m+=x[i];
    }

  jacob[0][0]=1.0;
  jacob[0][1]=2.0*alpha*m/((dim-2)*(dim-2));
  jacob[0][2]=2.0*alpha*m/((dim-2)*(dim-2));
  jacob[0][3]=-1.0;

  jacob[1][0]=0.0;
  jacob[1][1]=0.0;
  jacob[1][2]=0.0;
  jacob[1][3]=0.0;

  jacob[2][0]=0.0;
  jacob[1][1]=0.0;
  jacob[2][2]=0.0;
  jacob[2][3]=0.0;

  jacob[3][0]=0.0;
  jacob[3][1]=0.0;
  jacob[3][2]=0.0;
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
double LC= max(1.0,  alpha*STATE_MAX[1]);
/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
inline void localDynBounds(double * x, double * res)
{
  double m=0.0;
    for(int i=1;i<dim-1;i++)
      {
        m+=x[i];
      }
    m=alpha*m*m/((dim-2)*(dim-2));

  res[0]=  abs( x[0]+ m-x[3]);
  res[1]=0.0;
  res[2]=0.0;
  res[3]=  max(abs(CONTROL_MAX[0]),abs( CONTROL_MIN[0]));
}
/*!
 * \var M : global bound  for the dynamics if defined
 */
double M= max( alpha*STATE_MAX[1]*STATE_MAX[1]/2.0+STATE_MAX[3]+STATE_MAX[0],max(abs(CONTROL_MAX[0]),abs( CONTROL_MIN[0])));
/*!
 * \var computeM : indicates which method  to use to compute the bound for dynamics
 * 0= analytic  global value used
 * 1= local calculation using the localDynBounds function
 * 2= local calculation using  explicit maximization on controls
 */
const int computeM=1;



/*      *****************************************
 *  Definition of constraints and target
 *************************************************** */



/*!
 * \brief Function  defining the mixed  constraints
 *
 * This function defines the set U(x) for admissible controls as function of the state
 * @param x state variable
 * @param u control variable
 * @return  value that caraterise the constraints set
 */
inline double constraintsXU( double * x, double * u )
{

  return 1.0;
}


/*!
 * \brief Function  defining the state   constraints, corresonds  to k(x)
 *
 * This function defines the set K for admissible  states
 * @param x state variable
 * @return  value that caraterise the constraints set
 */


inline double constraintsX( double * x )
{
  double res=-PLUS_INF;
  res= max(res, max(0.2-x[0], x[0]-2.0));
  res= max(res, max(0.2-x[1], x[1]-2.0));
  res= max(res, max(0.2-x[2], x[2]-2.0));
  res= max(res,   -x[3]);

  double val=(res<=0)?1.0:PLUS_INF;

  return val;
}
/*!
 * Function that  characterise  the target set C, corresponds to c(x)
 * @param x state variable
 * @return value  to charaterise the target set
 */
inline  double target (double * x)
{
  double res =PLUS_INF;

  return res;
}

/* *   *************************************************************************
 *      Definition of value function
 ***************************************************************************/


/*!
 * Function  for optimisation criterion, corresponds to l(x,u)
 * @param x state variable
 * @param u control
 * @return value of l
 */
inline double l(double * x, double * u )
{
  return 1.0;
}



/*!
 * Function  for optimisation criterion, corresponds to m(x, u)
 * @param x state variable
 * @param u control
 * @return value of m
 */
inline double m(double * x, double * u )
{
  return 0.0;
}

void postProcess();
void postProcess()
{

}
#endif /* TESTDATA_H_ */
