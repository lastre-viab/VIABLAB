/*! \file  testZermelo_unused.h
 *

 *
 *
 */


#ifndef TESTZERMELO_UNUSED_H_
#define TESTZERMELO_UNUSED_H_


int compute_tmin=0;
const int dimc_ty=1;
double CONTROL_MIN_ty[dimc_ty]={0.0};
/*! \var CONTROL_MAX_ty[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX_ty[dimc_ty]={0.0};

unsigned long long int nbPointsControl_ty[dimc_ty] =  {2};

int ompThreads=1;

/*!
 * nature du probleme
 */
/*
 * sélection de type de reconstruction de trajectoire
 *  VD= 1, viable par défaut, on sélectionne le premier contrôle viable trouvé
 *  VL= 2, viable lourd: le contrôle reste constant tant qu'il viable;
 *  cette méthode nécessite une initialisation de contrôle
 */
int typeTraj=VD;

const int nbTrajs=0;
double initPoints[]={};
double initControls[]={};

                                             //    C : carr� [-1,1]x[-1,1]x[-1,1] de R3 ;
/*!
 * Function  defining the  switch conditions  between the continueos and  discrete dynamics
 * @param x state variable
 * @param u control variable
 * @return  value that caraterise the restet  set
 */
inline double resetSet( double * x, double * u )
{
        return 1.0;
}
/*!
 * \var dynType  defines the type of dynamics in the model
 *      1 or CC : continuous
 *      2 or DD : discrete
 *      3 or HD : hybrid
 */
 const int dynType=CC;

 /*!
  * Function that defines the  dynamical system
  */
   void dynamics_hybrid(double * x, double *u, double * image)
 {

 }


   unsigned int maxNbRetro=1048;
#endif /* TESTDATA_H_ */
