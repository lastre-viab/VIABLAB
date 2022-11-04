/*! \file  zermelo_tmin_data.h
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


#ifndef EQUILIBRES4D_DATA_H_
#define EQUILIBRES4D_DATA_H_
/*! \var dim
 *  \brief State  dimension
 */
const int dim=2;
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=1;
/*!
 * \var T maximum time horizon for the study
 */
double T=10.5;


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
const int discret_type=RK2I;


/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */

double STATE_MIN[dim]={-6.0, -2.5};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={1.0, 2.5};

/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int nbPointsState[dim]    =  {301, 301};

unsigned long long int dirTramage =0;

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0};

int saveProjection=0;
unsigned long long int projection[dim]={0,0};

int sortieOKinf[dim]={0, 0};
int sortieOKsup[dim]={0,0};
int intermediate_savings = 0;
int computeSet=1;
int saveBoundary=1;
int saveSubLevel=1;


double level=T;

/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MIN[dimc]={-pi};
/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX[dimc]={pi};

/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {501};




string prefix="zer-";

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
int gridMethod=MM;
/*
 * Sélection de l'ensemble à calculer
 */
int setType=CAPT;
int compute_tmin= 1;

const int nbTrajs=1;
double initPoints[dim*nbTrajs]={-5.0, 0.0};

double l_Lip = 1.0;
double l_max=1.0;
/*
 * Definition of the dynamics  and associated functions and constants
 */
double c=0.0, a=0.25;

void dynamics(double * x, double *u, double * image)
{
  image[0]= cos(u[0]) + c-a*tanh(x[1]);
  image[1]= sin(u[0]);
}

/*!
 * jacobian matrix of the dynamics
 */

inline void jacobian(double *x, double *u , double ** jacob)
{
  jacob[0][0]=0.0;
  jacob[0][1]=-a*(1.0-tanh(x[1])*tanh(x[1]));


  jacob[1][0]=0.0;
  jacob[1][1]=0.0;


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
  res[0]=1.0+c+a;
  res[1]=1.0;

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
int nbobs = 2; double obstacle[2*4]={-3.5, 1.0, 0.5, 1.0, -1.5, 0.0, 0.5, 0.75};

inline double constraintsX( double * x )
{
  double x1=x[0], x2=x[1];
   double res=-PLUS_INF, current;
   int i;
   for( i=0;i<nbobs;i++){
     current = -max(abs(x1-(obstacle[i*4]))-obstacle[i*4+2], abs(x2-(obstacle[i*4+1]))-obstacle[i*4+3]);
     res=max(res, current);
   }
   double resFinal=(res>0)?PLUS_INF:0.0;
   return resFinal;
}
/*!
 * Function that  characterise  the target set C, corresponds to c(x)
 * @param x state variable
 * @return value  to charaterise the target set
 */
inline  double target (double * x)
{
        double res;
        if(sqrt(x[0]*x[0]+x[1]*x[1])<=0.2)
        {
                res=0.0;
        }
        else
        {
                res=PLUS_INF;
        }
        return res;
}
#include "zermelo_tmin_unused.h"

#endif /* TESTDATA_H_ */
