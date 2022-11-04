/*! \file  ex3bis_Viabi2D_unused.h
 *
 */


#ifndef EX3BIS_VIABI2D_UNUSED_H_
#define EX3BIS_VIABI2D_UNUSED_H_

const int dimc_ty=1;
double CONTROL_MIN_ty[dimc_ty]={0.0};
/*! \var CONTROL_MAX_ty[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX_ty[dimc_ty]={0.0};
unsigned long long int nbPointsControl_ty[dimc_ty] =  {2};
int saveCoupeBound=0;
int saveCoupe=0;
int sliceDirs[dim]={  1, 0};
double sliceVals[dim]={  1.0, 0.0};
int saveProjection=0;
int sortieOK[dim]={0,0};
/*
 * sélection de type de reconstruction de trajectoire
 *  VD= 1, viable par défaut, on sélectionne le premier contrôle viable trouvé
 *  VL= 2, viable lourd: le contrôle reste constant tant qu'il viable;
 *  cette méthode nécessite une initialisation de contrôle
 */
int typeTraj=VD;

const int nbTrajs=0;
double initPoints[dim*nbTrajs]={};
double initControls[dimc*nbTrajs]={};
int scaling[dim]={0,0};
const int CaseCtrlDependsOnX=false; //true : u d�pend de x : lancement du programme KhiXU
void loadModelData();

/*!
 * Function  defining the  switch conditions  between the continueos and  discrete dynamics
 * @param x state variable
 * @param u control variable
 * @return  value that caraterise the restet  set
 */
inline double resetSet( double * x, double * u ){
  return 1.0;
}
/*!
 * \var dynType  defines the type of dynamics in the model
 *      1 or CC : continuous
 *      2 or DD : discrete
 *      3 or HD : hybrid
 */
const int dynType=CC;
int compteOnlyOptimalRetro=1;
 void loadModelData(){
}
/*!
 * Function that defines the  dynamical system
 */
void dynamics_hybrid(double * x, double *u, double * image){
}
/*!
 * Function that  characterise  the target set C, corresponds to c(x)
 * @param x state variable
 * @return value  to charaterise the target set
 */
inline  double target (double * x){
  return PLUS_INF;
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
inline double l(double * x, double * u ){
  return 1.0;
}
 /*!
 * Function  for optimisation criterion, corresponds to m(x, u)
 * @param x state variable
 * @param u control
 * @return value of m
 */
inline double m(double * x, double * u ){
  return 0.0;
}
 void postProcess();
void postProcess(){
}

int compute_tmin=0;
#endif /* EX3BIS_VIABI2D_H_ */
