/*! \file  testZermelo.h
 *
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

#include "testZermelo_unused.h"


 void loadModelData();

/*! \var dim
 *  \brief State  dimension
 */
 const int dim=2;
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=2;
/*! \var dicret_type
 *  \brief the  type of discretization scheme used  for the dynamics
 *  EE or 1 = Euler Explicit  scheme
 *  EI or 2 = Euler Implicit
 *  RK2I or 3 = RK2 Implicit (RK2 for -F)
 */
const int discret_type=RK2I;

/*!
 * \var T maximum time horizon for the study
 */
double T=10.0;

/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MIN[dimc]={-1.1, -1.1};
/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX[dimc]={1.1,1.1};

/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */
double STATE_MIN[dim]={-6.0, -5.0};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={2.0, 5.0};

int sortieOK[dim]={0,0};


/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {11,11};
/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int nbPointsState[dim]    =  {101,101};
/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0};
unsigned long long int projection[dim]={0,0};

unsigned long long int trajProjection[dim]={0,0};
 double initPoint[dim]={-4.6, 0.01};

 int scaling[dim]={0,0};


 int saveCoupeBound=0;
 int saveCoupe=0;
 int sliceDirs[dim]={  1, 0};
 double sliceVals[dim]={  1.0, 0.0};

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

 string prefix="test2-";

/*!
 * \var globalDeltaT
 *  bool�en indique si le pas de temps  doit �re choisi globalement
 *  ou localement pour les algorithmes de viabilit�
 */
bool globalDeltaT=false;


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



  double a=0.25;
  double r_c=0.44;
  double r_u=1.0;

  void dynamics_continuous(double * x, double *u, double * image)
{

	image[0]=  (1-a*x[1]*x[1])+u[0];
	image[1]=  u[1];
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
	jacob[0][0]=0.0;
	jacob[0][1]=-2.0*a*x[1];
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
double LC=2*a*max(STATE_MIN[1], STATE_MAX[1]);
/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
inline void localDynBounds(double * x, double * res)
{
	res[0]=  abs(1-a*x[1]*x[1])+max(abs(CONTROL_MAX[0]),abs( CONTROL_MIN[0]));
	res[1]=  max(abs(CONTROL_MAX[1]),abs( CONTROL_MIN[1]));
}
/*!
 * \var M : global bound  for the dynamics if defined
 */
double M= max(max(abs(CONTROL_MAX[1]),abs( CONTROL_MIN[1])),abs(1.0)+max(abs(CONTROL_MAX[0]),abs( CONTROL_MIN[0])));
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
	double res;
	if(u[0]*u[0]+u[1]*u[1]<=r_u*r_u)
	{
		res=1.0;
	}
	else
	{
		res=PLUS_INF;
	}
	return res;
}


/*!
 * \brief Function  defining the state   constraints, corresonds  to k(x)
 *
 * This function defines the set K for admissible  states
 * @param x state variable
 * @return  value that caraterise the constraints set
 */

 double cx1=-1.5;
 double cy1=3.0;
 double l1=0.5;
 double h1=0.5;

 double cx2=-3.5;
 double cy2=0.0;
 double l2=0.5;
 double h2=1.0;

inline double constraintsX( double * x )
{
	double res1= max(abs(x[0]-cx1)-l1, abs(x[1]-cy1)-h1);
	double res2= max(abs(x[0]-cx2)-l2, abs(x[1]-cy2)-h2);

	double res=(min(res1,res2)<0)?PLUS_INF:1.0;

	return res;
}
/*!
 * Function that  characterise  the target set C, corresponds to c(x)
 * @param x state variable
 * @return value  to charaterise the target set
 */
inline  double target (double * x)
{
	double res;
	if(sqrt(x[0]*x[0]+x[1]*x[1])<r_c)
	{
		res=0.0;
	}
	else
	{
		res=PLUS_INF;
	}
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
