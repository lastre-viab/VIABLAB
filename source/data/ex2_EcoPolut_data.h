/*! \file  testZermelo.h
 *
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


void loadModelData();

/*! \var dim
 *  \brief State  dimension
 */
const int dim=4;
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=2;
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
double T=10.0;

/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MIN[dimc]={0.0, 0.0};
/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX[dimc]={1.0, 1.0};

/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */
double STATE_MIN[dim]={0.0, 0.0, 0.0, 0.0};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={1., 1.0, 1.0, 1.0};

/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {5, 5};
/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int nbPointsState[dim]    =  {11,11,11, 11};
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
int dirTramage=3;
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

/*
 * sélection de type de reconstruction de trajectoire
 *  VD= 1, viable par défaut, on sélectionne le premier contrôle viable trouvé
 *  VL= 2, viable lourd: le contrôle reste constant tant qu'il viable;
 *  cette méthode nécessite une initialisation de contrôle
 */
int typeTraj=VL;

/*

const int nbTrajs=5;
double initPoints[dim*nbTrajs]={1.2, 1.0, 1.0, 1.5, 1.4, 1.8, 0.4, 0.2, 0.5, 0.6};
double initControls[dimc*nbTrajs]={0.0, 0.0, 0.0, 0.0, 0.0};

 */
const int nbTrajs=0;
double initPoints[dim*nbTrajs]={};
double initControls[dimc*nbTrajs]={};


// const int nbTrajs=1;
// double initPoints[dim*nbTrajs]={ 0.5, 0.6};


/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0};
unsigned long long int projection[dim]={0};

unsigned long long int trajProjection[dim]={0};

int scaling[dim]={0, 0, 0, 0};


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

string prefix="ex1_eco-fin";

/*!
 * \var globalDeltaT
 *  bool�en indique si le pas de temps  doit �re choisi globalement
 *  ou localement pour les algorithmes de viabilit�
 */
bool globalDeltaT=false;

//          A = autre cas de la forme f(x,u,v)
//                  ou u=(u[1],u[2],u[3],u[4],u[5]) in [alphamin,alphamax]x[betamin,betamax]x[gammamin,gammamax] }

/*
 * Definition of the dynamics  and associated functions and constants
 */


int compteOnlyOptimalRetro=1;


void loadModelData()
{

}

double a=0.06;
double rho=0.05;
double d=4.0;
double b=0.08;

inline double prodFunc(double K)
{
  double x=max(K,0.0);
  return  sqrt(K);
}
inline double prodFuncPrime(double K)
{
  double x=max(K,0.000000001);

  return  1.0/(sqrt(x)*2.0);
}
inline double utility(double C)
{
  return   (sqrt(C) );
  // return  log(C);
}
inline double utilityPrime(double C)
{
  return   (0.5/sqrt(max(C, 0.000000001)) );
  // return  log(C);
}




void dynamics(double * x, double *u, double * image)
{
  double alpha=u[0];
  double beta=u[1];

  double K=x[0];
  double P=x[1];
  double z1=x[2];
  double z2=x[3];

  image[0]=  (1.0-alpha-beta)*prodFunc(K)-a*K;
  image[1]=  (1-beta*d)*prodFunc(K)-b*P;
  image[2]=  0.0;//rho*z1+utility(alpha*prodFunc(K));
  image[3]= rho*z2-pow(P, 1.1)/1.1;
  //cout<< " dynamique renvoie "<<image[0]<< " "<<image[1]<<endl;
}



/*!
 * jacobian matrix of the dynamics
 */
inline void jacobian(double *x, double *u , double ** jacob)
{
  double alpha=u[0];
  double beta=u[1];

  double K=x[0];
  double P=x[1];
  double z1=x[2];
  double z2=x[3];

  jacob[0][0]=prodFuncPrime(K)-a;
  jacob[0][1]=0.0;
  jacob[0][2]=0.0;
  jacob[0][3]=0.0;

  jacob[1][0]=(1-beta*d)*prodFuncPrime(K);
  jacob[1][1]=-b;
  jacob[1][2]=0.0;
  jacob[1][3]=0.0;

  jacob[2][0]=0.0;
  jacob[2][1]=0.0;
  jacob[2][2]=0.0;
  jacob[2][2]=0.0;

  jacob[3][0]=0.0;
  jacob[3][1]=-pow(P, 0.1);
  jacob[3][2]=0.0;
  jacob[3][2]=rho;

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
inline void localDynBounds(double * x, double * image)
{
  double K=x[0];
  double P=x[1];
  double z1=x[2];
  double z2=x[3];

  image[0]=  abs(prodFunc(K)-a*K);
  image[1]= max(abs( (1- d)*prodFunc(K)-b*P),abs( prodFunc(K)-b*P));
  image[2]= abs( rho*z1+utility(prodFunc(K)));
  image[3]= abs(rho*z2-pow(P, 1.1));
}


/*!
 * \var M : global bound  for the dynamics if defined
 */
double M=3.0;// max( alpha*STATE_MAX[0]+STATE_MAX[1],max(abs(CONTROL_MAX[0]),abs( CONTROL_MIN[0])));
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
  double alpha=u[0];
  double K=x[0];
  double beta=u[1];
double res=max( -utility(alpha*prodFunc(K))+x[2],u[0]+u[1]-1.0);
  return (res<=0)?1.0:PLUS_INF;
}


/*!
 * \brief Function  defining the state   constraints, corresonds  to k(x)
 *
 * This function defines the set K for admissible  states
 * @param x state variable
 * @return  value that caraterise the constraints set
 */

double Wmin=0.0;
inline double constraintsX( double * x )
{
  /*double res1= -x[2];


  double res=(res1<=0)?1.0:PLUS_INF;
*/
  return 1.0;
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

#include "ex2_EcoPolut_unused.h"

#endif /* TESTDATA_H_ */
