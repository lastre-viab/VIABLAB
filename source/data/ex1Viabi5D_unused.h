/*! \file  testZermelo_unused.h
 *

 *
 *
 */


#ifndef TESTZERMELO_UNUSED_H_
#define TESTZERMELO_UNUSED_H_



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
const int probNature=CP;

const int CaseCtrlDependsOnX=false; //true : u d�pend de x : lancement du programme KhiXU


const int FormeCorrespondance=CR3;      //{b : boule de R2  ; B : boule de R3 ; c : carr� [-1,1]x[-1,1] de R2 ;
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
 const int dynType=CD;

extern  void dynamics_continuous(double * x, double *u, double * image);
extern void dynamics_discrete(double * x, double *u, double * image);
 /*!
  * Function that defines the  dynamical system
  */
   void dynamics_hybrid(double * x, double *u, double * image)
 {
          if(resetSet(x,u))
                  dynamics_continuous(x,u,image);
          else
                  dynamics_discrete(x,u,image);
 }


   unsigned int maxNbRetro=1048;
#endif /* TESTDATA_H_ */
