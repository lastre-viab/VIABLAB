/*! \file  ex3bis_Viabi2D_data.h
 *
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
 *  \author: A. Désilles, LATSRE
 *  \brief  Fichier contant la définition d'un modèle de viabilité pour les tests
 *
 *  Ce fichier  contient les définitions  de tous les éléments qui décrivent un
 *  un problème de viabilité concret et  qui sont considérés comme
 *  paramètres par le code; Ce fichier représente ainsi une interface d'entrée primitive pour le code.
 *
 *  Les éléments principaux définis ici  sont :
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
 */
#ifndef EX3BIS_VIABI2D_DATA_H_
#define EX3BIS_VIABI2D_DATA_H_
/*
 * paramètres du modèle
 */
double r=1.5;      // rayon de la contrainte
double alpha=1.0;  // rayon de contrôles
int d=4;
/*! \var dim
 *  \brief State  dimension
 */
const int dim=d;
/*! \var dicret_type
 *  \brief the  type of discretization scheme used  for the dynamics
 *  EE or 1 = Euler Explicit  scheme
 *  EI or 2 = Euler Implicit
 *  RK2I or 3 = RK2 Implicit (RK2 for -F)
 *  RK2E or 4 = RK2 Explicit (RK2 for F)
 */
const int discret_type=RK2E;
/*!
 * \var T maximum time horizon for the study
 */
double T=10.0; // non utilisé ici
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=1;
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
/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {181};
/*! \var STATE_MIN[i]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */
double *STATE_MIN;
/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double *STATE_MAX;
/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int *nbPointsState;
string prefix="ex1_viab-d4-";
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
 * Choix de l'axe qui sera utilisé pour la définition des segments
 * de maillage ; Dans ce modèle, parfaitement symétrique, ce choix ne
 * joue pas de rôle sugnificatif pour la vitesse de convergence
 */
int dirTramage=2;
/*
 * Nombre d'étapes de rafinement
 */
int refine=2;
/*
 * Paramètre qui définit si l'ensemble doit être recalculé.
 * Si computeSet=1 alors l'ensemble sea calculé, si
 */
int computeSet=1;
int ompThreads=1;
int saveBoundary=1;
int saveCoupeBound=1;
int saveCoupe=1;
int *sliceDirs;
double *sliceVals;

int *sortieOK;
int *scaling;

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int *periodic;
unsigned long long int *projection;

unsigned long long int *trajProjection;
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

/*!
 * \var globalDeltaT
 *  bool�en indique si le pas de temps  doit �re choisi globalement
 *  ou localement pour les algorithmes de viabilit�
 */
bool globalDeltaT=false;
/*
 * Definition of the dynamics  and associated functions and constants
 */


unsigned long long int NBS=21;

void loadModelData();
void loadModelData(){


  STATE_MIN =new double [dim];
  STATE_MAX=new double [dim];
  nbPointsState=new unsigned long long int[dim];

  periodic=new int[dim];
  projection=new unsigned long long int[dim];
  trajProjection=new unsigned long long int[dim];
  sliceDirs=new int[dim];
  sliceVals =new double [dim];
  sortieOK=new int[dim];
  scaling=new int[dim];
  for(int i=0;i<d;i++)
    {
   STATE_MIN[i]=0.0;
    STATE_MAX[i]=1.0;
    nbPointsState[i]=NBS;

    periodic[i]=0;
    projection[i]=0;
    trajProjection[i]=0;
    scaling[i]=0;
    sortieOK[i]=0;
    }
  STATE_MAX[0]=2.0;
  STATE_MAX[1]=3.0;

  nbPointsState[0]=41;
  nbPointsState[1]=61;



  sliceDirs[0]=0;
  sliceDirs[1]=0;
  sliceVals[0]=0.0;
  sliceVals[1]=0.0;
  for(int i=2;i<d;i++)
    {
    sliceDirs[i]=1;
    sliceVals[i]=0.5;
    }
}
void dynamics(double * x, double *u, double * image)
{
  double m=0.0;
  for(int i=2;i<dim;i++)
    {
    m+=x[i];
    }
  m=m*m/((dim-2)*(dim-2));
  image[0]=   x[0]+m-x[1];
  image[1]=  u[0];
  for(int k=2;k<dim;k++)
    {
    image[k]=  0.0;
    }


}
/*!
 * jacobian matrix of the dynamics
 */
inline void jacobian(double *x, double *u , double ** jacob)
{
  double m=0.0;
  for(int i=2;i<dim;i++)
    {
    m+=x[i];
    }
  jacob[0][0]=1.0;
  jacob[0][1]=-1.0;

  for(int j=2;j<d;j++)
    {
    jacob[0][j]=2.0*m/((dim-2)*(dim-2));
    }

  for(int i=1;i<d;i++)
    {
    for(int j=0;j<d;j++)
      {
      jacob[i][j]=0.0;
      }
    }
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
double LC= 1.0;
/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
inline void localDynBounds(double * x, double * res)
{
  double m=0.0;
  for(int i=2;i<dim;i++)
    {
    m+=x[i];
    }
  m=m*m/((dim-2)*(dim-2));
  res[0]=  abs(  m+x[0]-x[1]);
  res[1]=  max(abs(CONTROL_MAX[0]),abs( CONTROL_MIN[0]));
  for(int i=2;i<d;i++)
    res[i]=  0.0;
}
/*!
 * \var M : global bound  for the dynamics if defined
 */
double M=2.5;// max( alpha*STATE_MAX[0]+STATE_MAX[1],max(abs(CONTROL_MAX[0]),abs( CONTROL_MIN[0])));
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
  /*double res1= max(0.0-x[0], x[0]-2.0);
        double res2=  -x[1];

        double res=(max(res1,res2)<=0)?1.0:PLUS_INF;
   */
  return 1.0;
}

#include "ex1Viabi_MultiDim_unused.h"
#endif /* EX3BIS_VIABI2D_DATA_H_ */
