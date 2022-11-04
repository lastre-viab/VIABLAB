/*! \file  Equilibres4D_data.h
 *
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
const int discret_type=RK2E;

double localMu_Max=0.02;
double localMu_Min=-0.006213;


double localAlpha=1;
double localBeta=1;
double localArgMuMax = 0.2;

/*
 * Toujours la même valeur, il vaut mieux la déclarer comme globale
 * et ne pas la recalculer à chaque appel de la dynamqiue et des autres fonctions
 */
double k = (localMu_Max - localMu_Min)*1.0 / localArgMuMax;


double localA = 1; 
double localK = 1;
double localR = 1;

double localP = 0.7772;
double localDR = 0.45;
double localDC = 0.112;

// Work area
//double localXMin,localXMax,localYMin,localYMax,localZMin,localZMax;
double localXMin = 1.9;
double localXMax = 2.2;
double localYMin = 0.1;
double localYMax = 0.4;
double localZMin = 0.3;
double localZMax = 0.5;

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



int res=600;
int resX=res;
int resY=res;
int resZ=400;

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
unsigned long long int nbPointsState[dim] = {resX,resY,resZ};

unsigned long long int dirTramage =0;
/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0,0};

unsigned long long int projection[dim]={0,0,0};

string prefix="ThibautNewVersionRes"+to_string(res);

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



inline double muFonction(double w){

  double muRate = 0;
  if (w  <= localArgMuMax){
    muRate = localMu_Min + k * w;
  }else{
    muRate = localMu_Max * exp(localArgMuMax-w);
  }
  return muRate;
}


void dynamics(double * x, double *u, double * image)
{
  double naturalCapitalConsumption = localDR * localP*x [0] * x[1];

  // biens manufactures
  double prodZ,consZ;
  prodZ = localA * exp(localAlpha * log ((1-localP)*x[0]) ) * exp(localBeta*log(naturalCapitalConsumption));
  consZ = localDC*x[0]*x[2];

  //population humaine
  double w = consZ / x[0];

  image[0] = x[0] * muFonction(w);
  // capital naturel
  image[1] = localR*x[1] * ( localK - x[1]) - naturalCapitalConsumption ;

  image[2] = prodZ - consZ;



}


inline double muFonctionPrime(double w){
  double k = (localMu_Max - localMu_Min)*1.0 / localArgMuMax;

  if ( w <= localArgMuMax){
    return k;
  } else {
    return -muFonction(w);
  }
}

/*!
 * jacobian matrix of the dynamics
 */

inline void jacobian(double *x, double *u , double ** jacob)
{

  double  w = localDC * x[2];


  jacob[0][0] = muFonction(w) ;
  jacob[0][1] = 0;
  jacob[0][2] = muFonctionPrime(w) * localDC*x[0];

  jacob[1][0] = -localDR*localP*x[1];
  jacob[1][1] = -2*localR*x[1]+localR*localK-localDR*localP*x[0];
  jacob[1][2] = 0;

  jacob[2][0] = localA*(localAlpha+localBeta)*exp(localAlpha*log(1-localP))*exp(localBeta*log(localDR*localP))*exp((localAlpha+localBeta-1)*log(x[0]))*exp(localBeta*log(x[1]))-w;
  jacob[2][1] = localA*localBeta*exp(localAlpha*log(1-localP))*exp(localBeta*log(localDR*localP))*exp((localAlpha+localBeta)*log(x[0]))*exp((localBeta-1)*log(x[1]));
  jacob[2][2] = -localDC*x[0]; 

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
double LC= 20;

/*!
 * \var M : global bound  for the dynamics if defined
 */
double M= 20;
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

   /*
    * pour améliorer les performances ici
    * il vaut mieux recopier le code de la dynamquue
    * cela évite d'allouer un tableau à chaque appel
    */
  double naturalCapitalConsumption = localDR * localP*x [0] * x[1];

   // biens manufactures
   double prodZ,consZ;
   prodZ = localA * exp(localAlpha * log ((1-localP)*x[0]) ) * exp(localBeta*log(naturalCapitalConsumption));
   consZ = localDC*x[0]*x[2];

   //population humaine
   double w = consZ / x[0];

   res[0] = x[0] * muFonction(w);
   // capital naturel
   res[1] = localR*x[1] * ( localK - x[1]) - naturalCapitalConsumption ;

   res[2] = prodZ - consZ;

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

  double res=-PLUS_INF;
  for(int  i=0;i<dim;i++){
    res = max(res, max(STATE_K_MIN[i]-x[i], x[i]-STATE_K_MAX[i]));
  }
  double resFinal=(res>0)?PLUS_INF:0.0;
  return resFinal;


}

#include "ChasseurCueilleur_unused.h"

#endif /* TESTDATA_H_ */
