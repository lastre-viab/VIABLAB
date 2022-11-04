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

double alpha = 1.0;
double r = 3.0;



////////////////////////////////////


/*! \var dim
 *  \brief State  dimension
 */
const int dim=2;


// controles 
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=2;

/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
//double CONTROL_MIN[dimc]={-(r+alpha),-(r+alpha)};
double CONTROL_MIN[dimc]={-alpha,-alpha};

/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
//double CONTROL_MAX[dimc]={r+alpha,r+alpha};
double CONTROL_MAX[dimc]={alpha,alpha};

/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {101,101};

 double pasC_0 = double( (CONTROL_MAX[0] - CONTROL_MIN[0]) / (nbPointsControl[0]-1))  ;
 double pasC_1 = double( (CONTROL_MAX[1] - CONTROL_MIN[1]) / (nbPointsControl[1]-1))  ;

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
  double norme_x = sqrt(x[0]*x[0] + x[1]*x[1]);
  double norme_u = sqrt(u[0]*u[0] + u[1]*u[1]);
  double res = 1.0;
  /*
  if ( norme_u > norme_x + alpha ) {res = PLUS_INF;
  } else {res = 1.0;};
   */
  if ( norme_u > alpha ) {res = PLUS_INF;} else {res = 1.0;};
  return res;
}




/*!
 * \var T maximum time horizon for the study
 */
double T= 100.0;
//int T= 100;



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
double STATE_MIN[dim]={-r, -r};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={r, r};
// TODO mettre aire max des données

/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
//unsigned long long int nbPointsState[dim]    =  {401,401, 11};
unsigned long long int nbPointsState[dim]    =  {100,100};


double pas_0 = double( (STATE_MAX[0] - STATE_MIN[0]) / (nbPointsState[0]-1))  ;
double pas_1 = double( (STATE_MAX[1] - STATE_MIN[1]) / (nbPointsState[1]-1))  ;


/*
 *
 *
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
int refine= 0; //2;

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0};

/*!
 * \var dbFileName
 * the  name for the  main data base file for the roject
 */

string prefix="F_test_tyche" ;
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







void dynamics(double * x, double *u, double * image);


void dynamics(double * x, double *u, double * image){

  // bornes tychastiques:
  double norme_x = sqrt(x[0]*x[0] + x[1]*x[1]);

  double temp_im0 = 0.0 ;
  double temp_im1 = 0.0 ;



  // indep de norme v pour commencer
  temp_im0 = 2*x[0] + u[0]  ;
  temp_im1 = 2*x[1] + u[1]  ;


  //image[0] = temp_im0 ;
  //image[1] = temp_im1 ;

  //////// projections

  // 0


  double pas_0 = double( (STATE_MAX[0] - STATE_MIN[0]) / nbPointsState[0] )  ;

  int k_0 = floor( (temp_im0 - STATE_MIN[0]) / pas_0 ) ;

  image[0] = STATE_MIN[0] + (k_0 * pas_0) ;



  // 1

  double pas_1 = double( (STATE_MAX[1] - STATE_MIN[1]) / nbPointsState[1] )  ;

  int k_1 = floor( (temp_im1 - STATE_MIN[1]) / pas_1) ;

  image[1] = STATE_MIN[1] + (k_1 * pas_1) ;

}

inline double  constraintsX_fd( unsigned long long int * x )
{
  double xReel0=STATE_MIN[0]+x[0]* pas_0;
  double xReel1=STATE_MIN[1]+x[1]* pas_1;

  double norme_x = sqrt(xReel0*xReel0 + xReel1*xReel1);

  double res = 1.0;
  if(norme_x <= r ) {res = 1.0;} else {res = PLUS_INF;};
  return res;
  //return 1.0;
  return  res;
}

void dynamics_fd(unsigned long long int  * x, unsigned long long int *u, unsigned long long int * image)
{
  double norme_x = sqrt(x[0]*x[0] + x[1]*x[1]);

  double temp_im0 = 0.0 ;
  double temp_im1 = 0.0 ;

  double xReel0=STATE_MIN[0]+x[0]* pas_0;
  double xReel1=STATE_MIN[1]+x[1]* pas_1;

  double uReel0=CONTROL_MIN[0]+u[0]* pasC_0;
  double uReel1=CONTROL_MIN[1]+u[1]* pasC_1;
  // indep de norme v pour commencer
  temp_im0 = 2*xReel0 + uReel0  ;
  temp_im1 = 2*xReel1 + uReel1  ;

  //image[0] = temp_im0 ;
  //image[1] = temp_im1 ;

  //////// projections

  // 0

  image[0] = (unsigned long long int)floor( (temp_im0 - STATE_MIN[0]) / pas_0 ) ;

  // 1

  image[1] = (unsigned long long int)floor( (temp_im1 - STATE_MIN[1]) / pas_1) ;
}


inline double  constraintsXU_fd( unsigned long long int * x, unsigned long long int * u )
{
  double res=1.0;

  double uReel0=CONTROL_MIN[0]+u[0]* pasC_0;
  double uReel1=CONTROL_MIN[1]+u[1]* pasC_1;
   double norme_u = sqrt(uReel0*uReel0 + uReel1*uReel1);

  if ( norme_u > alpha ) {res = PLUS_INF;} else {res = 1.0;};
   return res;
}

///////// tychastique
int dimc_ty=2;

double *CONTROL_MIN_ty;
/*! \var CONTROL_MAX_ty[dimc]
 *  \brief maximum values  for  tyches vector components
 *   \see controlParams
 */
double *CONTROL_MAX_ty;

unsigned long long *nbPointsControl_ty;

void dynamics_tych_fd(unsigned long long int  * x, unsigned long long int *u,unsigned long long int *v, unsigned long long int * image)
{
  double norme_x = sqrt(x[0]*x[0] + x[1]*x[1]);

   double temp_im0 = 0.0 ;
   double temp_im1 = 0.0 ;

   double xReel0=STATE_MIN[0]+x[0]* pas_0;
   double xReel1=STATE_MIN[1]+x[1]* pas_1;

   double uReel0=CONTROL_MIN[0]+u[0]* pasC_0;
   double uReel1=CONTROL_MIN[1]+u[1]* pasC_1;
   // indep de norme v pour commencer
   temp_im0 = 2*xReel0 + uReel0  ;
   temp_im1 = 2*xReel1 + uReel1  ;

   //image[0] = temp_im0 ;
   //image[1] = temp_im1 ;

   //////// projections

   // 0

   image[0] = (unsigned long long int)floor( (temp_im0 - STATE_MIN[0]) / pas_0 ) ;

   // 1

   image[1] = (unsigned long long int)floor( (temp_im1 - STATE_MIN[1]) / pas_1) ;

}

///////// valeurs
/*
CONTROL_MIN_ty = new double [dimc_ty];
CONTROL_MIN_ty={-r,-r};
 */
/// test perturbations


/*
void dynamics_ty(double * x, double *u,  double *v, double * image);
void dynamics_ty(double * x, double *u,  double *v, double * image)
{

  // bornes tychastiques:
  double norme_x = sqrt(x[0]*x[0] + x[1]*x[1]);
  double norme_v = sqrt(v[0]*v[0] + v[1]*v[1]);

  double temp_im0 = 0.0 ;
  double temp_im1 = 0.0 ;


  // indep de norme v pour commencer
  temp_im0 = 2*x[0] + u[0]  ;
  temp_im1 = 2*x[1] + u[1]  ;



  //////// projections

  // 0

  double pas_0 = double( (STATE_MAX[0] - STATE_MIN[0]) / nbPointsState[0] )  ;

  int k_0 = floor( (temp_im0 - STATE_MIN[0]) / pas_0 ) ;

  image[0] = STATE_MIN[0] + (k_0 * pas_0) ;



  // 1

  double pas_1 = double( (STATE_MAX[1] - STATE_MIN[1]) / nbPointsState[1] )  ;

  int k_1 = floor( (temp_im1 - STATE_MIN[1]) / pas_1) ;

  image[1] = STATE_MIN[1] + (k_1 * pas_1) ;

}

 */


/*
  if ( norme_v > norme_x) {
  // si perturbation hors de l'ensemble considéré
  temp_im0 = x[0] + u[0]  ;
  temp_im1 = x[1] + u[1]  ; } else {

  // si perturbation hors de l'ensemble considéré
  temp_im0 = x[0] + u[0] + v[0] ;
  temp_im1 = x[1] + u[1] + v[1] ;
  };
 */








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
  double norme_x = sqrt(x[0]*x[0] + x[1]*x[1]);

  double res = 1.0;
  if(norme_x <= r ) {res = 1.0;} else {res = PLUS_INF;};
  return res;
  //return 1.0;
}





////////////
void loadModelData();


void loadModelData(){


  int dimc_ty=2;

  CONTROL_MIN_ty = new double [dimc_ty];
  CONTROL_MIN_ty[0]= -r;
  CONTROL_MIN_ty[1]= -r;
  //CONTROL_MIN_ty={-r,-r};
  //CONTROL_MIN_ty[dimc_ty]={-r,-r};

  CONTROL_MAX_ty = new double [dimc_ty];
  CONTROL_MAX_ty[0]= r;
  CONTROL_MAX_ty[1]= r;
  //CONTROL_MAX_ty[dimc_ty]={r,r};

  nbPointsControl_ty = new unsigned long long [dimc_ty];
  nbPointsControl_ty[0] = 2;
  nbPointsControl_ty[1] = 2;
  //nbPointsControl_ty = {20,20};

}



#include "test_tycha_discret_unused.h"
#endif /* TESTDATA_H_ */
