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
//control max
unsigned long long int Um=10;

//target

unsigned long long int N =20;

//Domain

unsigned long long int Xmax =2*N;
unsigned long long int Ymax= 2*N;

//time horizon

unsigned long long int Tmax = 20;




////////////////////////////////////


/*! \var dim
 *  \brief State  dimension
 */
const int dim=3;


// controles 
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=1;

/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
//double CONTROL_MIN[dimc]={-(r+alpha),-(r+alpha)};
double CONTROL_MIN[dimc]={0};

/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
//double CONTROL_MAX[dimc]={r+alpha,r+alpha};
double CONTROL_MAX[dimc]={(double) (Um)};

/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */

unsigned long long int nbPointsControl[dimc] =  {Um+1};




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
//const int dynType=DC;
const int dynType=DD;


/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */
double STATE_MIN[dim]={0.0, 0.0, 0.0};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={(double) Xmax, (double) Ymax, (double)Tmax};
// TODO mettre aire max des données

/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
//unsigned long long int nbPointsState[dim]    =  {401,401, 11};
unsigned long long int nbPointsState[dim]    =  {Xmax+1,Ymax+1, Tmax+1};


/*
 *
 *
 * Direction le long de laquelle l'ensemble sera représenté par des suites de bits
 */
unsigned long long int dirTramage =2;

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
int refine= 0; //2;

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0,0};

/*!
 * \var dbFileName
 * the  name for the  main data base file for the roject
 */

string prefix="CaptDiscret-" ;
//string prefix= std::to_string(i) + "_F_renouee_dispersion-" ;


/*
 * Definition of the dynamics  and associated functions and constants
 */

/*
 * Paramètre qui indique si l'ensemble doit être calculé
 * dans le cas contrainte il devra pouvoir être chargé depuis un fichier pour une étude de
 * trajectoires
 */
int computeSet=1;

//int saveBoundary=1;
int saveBoundary=1;

const int nbTrajs=0;
double initPoints[dim*nbTrajs]={};

int sortieOKinf[dim]={0,0};
int sortieOKsup[dim]={0,0};





inline double  constraintsX_fd( unsigned long long int * x )
{

	double res = 1.0;

	if(x[2]>=Tmax)
		res = (x[0]>=N && x[1]>=N)? 1.0: PLUS_INF;
	else
		res = (x[1]>=x[0])? 1.0: PLUS_INF;

	return res;
}

void dynamics_fd(unsigned long long int  * x, unsigned long long int *u, unsigned long long int * image)
{

if(x[2]<Tmax)
{
	image[0] = x[0]+min(x[1]-x[0],(unsigned long long int)floor(0.5*x[0])) ;

	// 1

	image[1] =x[1]+ u[0];

	image[2] = x[2]+1;
}
else
{
	image[0] = x[0];

		// 1

		image[1] =x[1];

		image[2] = x[2];
}
}


inline double  constraintsXU_fd( unsigned long long int * x, unsigned long long int * u )
{
	double res=1.0;
	return res;
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

////////////
void loadModelData();


void loadModelData(){

}



#include "test_capt_discret_unused.h"
#endif /* TESTDATA_H_ */
