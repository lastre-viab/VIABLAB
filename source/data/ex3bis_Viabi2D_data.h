/*! \file  ex3bis_Viabi2D_data.h
 *
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

 /*! \var dim
 *  \brief State  dimension
 */
const int dim=2;
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
const int dimc=2;
/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MIN[dimc]={-1.0, -1.0};
/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX[dimc]={1.0, 1.0};
/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {101, 101};
/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */
double STATE_MIN[dim]={-r, -r};
/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={r, r};
/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int nbPointsState[dim]    =  {151, 151};
string prefix="ex3bis_viab2D-refine5-";
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
int dirTramage=0;
/*
 * Nombre d'étapes de rafinement
 */
int refine=0;
/*
 * Paramètre qui définit si l'ensemble doit être recalculé.
 * Si computeSet=1 alors l'ensemble sea calculé, si
 */
int computeSet=1;
int ompThreads=1;
int saveBoundary=1;

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0};
unsigned long long int projection[dim]={0,0};

unsigned long long int trajProjection[dim]={0,0};
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
void dynamics(double * x, double *u, double * image)
{
  image[0]=   x[0]+alpha*u[0];
  image[1]=  x[1]+alpha*u[1];
}
/*!
 * jacobian matrix of the dynamics
 */
inline void jacobian(double *x, double *u , double ** jacob)
{
  jacob[0][0]=1.0;
  jacob[0][1]=0.0;
  jacob[1][0]=1.0;
  jacob[1][1]=0.0;
}
/*!
 * \var computeLC : indicates which method  to use to copute Lipschitz constant
 * 0= analytic  global value used
 * 1= local calculation using the jacobian matrix
 * 2= local calculation using the finite differences
 */
const int computeLC=0;
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
  res[0]=  abs(x[0])+alpha;
  res[1]= abs(x[1])+alpha;
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
inline double constraintsXU( double * x, double * u ){
  double res=(u[0]*u[0]+u[1]*u[1]-1.0<=0.0)?1.0:PLUS_INF;
  return res;
}
/*!
 * \brief Function  defining the state   constraints, corresonds  to k(x)
 * This function defines the set K for admissible  states
 * @param x state variable
 * @return  value that caraterise the constraints set
 */
inline double constraintsX( double * x ){
  double res1=  x[0]*x[0]+x[1]*x[1]-r*r;
  double res=( res1<=0)?1.0:PLUS_INF;
  return res;
}

#include "ex3bis_Viabi2D_unused.h"
#endif /* EX3BIS_VIABI2D_DATA_H_ */
