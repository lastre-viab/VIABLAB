/*! \file  testZermelo.h
 *
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

/*
// params general
double L = 3 ;   // >0 
double s = 1 ;   // >0
double delta = 0.5 ; // in ]0,1[ 
double h = 0.3 ;   // in ]0,1[


// params indivs
double W_i = 5 ;   // in [1,10]
double W = 8*W_i ;  // in [W_i,10*W_i]
double C_i = L * W_i / (2* (1+ exp(s*h)) ) ;    // in [ L * W_i / (2* (1+e^(s*h)) ) ;  L * W_i / (2* (1+e^(s*(h-1))) ) ]
*/




// params generaux
double L = 1.0 ;   // >0 
double s = 0.2 ;   // >0
double delta = 0.5 ; // in ]0,1[ 
double h = 0.2 ;   // in ]0,1[


// params indivs
double W_i = 5.0 ;   // in [1,10]
double W = 8.0*W_i ;  // in [W_i,10*W_i]
double C_i = 2.0 ;    // in [ L * W_i / (2* (1+e^(s*h)) ) ;  L * W_i / (2* (1+e^(s*(h-1))) ) ]







double e_inf = double(L)  / double(1+exp(s*(h)))  ;  
double e_sup = double(L)  / double(1+exp(s*(h-1))) ;

/*

// functions indivs
// benefits
void b_i(double * e, double * image);
void b_i(double * e, double * image)
{
  image[0]=  W_i * e[0];
}

// costs
void c_i(double * u_i, double * image);
void c_i(double * u_i, double * image)
{
  image[0]=  C_i / (1- u_i[0]/2);
}

*/


// durée (horison temporel)
//double Tfinal = 10;  
double Tfinal = 10000;  




////////////////////////////////////


/*! \var dim
 *  \brief State  dimension
 */
const int dim=2;

// controles 
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=1;

/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MIN[dimc]={0};

/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX[dimc]={1};

/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {2};

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
    // doit renvoyer + \infty  si u \not \in U(x).

    double res = 1.0;
    if ( u[0] < 0.0  ) { res = PLUS_INF;
                  } else if ( u[0] > 1 ) {res = PLUS_INF;
                  } else {res = 1.0;};
   return res;
   }




// tychastic (voir unused)


/*!
 * \var T maximum time horizon for the study
 */
double T= Tfinal;


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
double STATE_MIN[dim]={e_inf, 0};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={e_sup, 10*e_sup};

/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
//unsigned long long int nbPointsState[dim]    =  {401,401, 11};
unsigned long long int nbPointsState[dim]    =  {51,51};

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
int refine= 2;

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0};

/*!
 * \var dbFileName
 * the  name for the  main data base file for the roject
 */

string prefix="F_IPGG-";

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

int sortieOKinf[dim]={0, 0};
int sortieOKsup[dim]={1,1};


void dynamics(double * x, double *u, double * image);
void dynamics(double * x, double *u, double * image)
{
  // x[0] : e, x[1]:k_i
  image[0]= x[0]+1.0;//delta * x[0] + (1-delta) * double(L) / double(1+exp(s* (h- W_i * double(u[0]) / double(W) ) ) ) ;

  image[1]= x[1];//x[1] + W_i * image[0] - double(C_i) / double(1- double(u[0])/ double(2) ) ;
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

  double res = 1.0;
  if(x[1] < 0) {res = PLUS_INF;} else {res = 1.0;};
  return res;

}




// pour essayer de résoudre le pb de non mise a jour des params
void static updateParams(double L_, double s_, double h_){

    e_inf = double(L_)  / double(1+exp(s_*(h_)))  ;  
    e_sup = double(L_)  / double(1+exp(s_*(h_-1))) ;
    STATE_MIN[0]=e_inf;
    STATE_MIN[1]=0;
    STATE_MAX[0]=e_sup;
    STATE_MAX[1]=10*e_sup;
printf("merde");
printf("sort[0] %d",sortieOKsup[0]);
printf("sort[1] %d",sortieOKsup[1]);

getchar();

}



#include "IPGG_unused.h"
#endif /* TESTDATA_H_ */
