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


int saveProjection=0;
int saveCoupeBound=1;
int saveCoupe=1;

int sliceDirs[dim]={  1, 1, 0, 0};
double sliceVals[dim]={ 0.5, 0.5, 0.0, 0.0};

int intermediate_savings = 0;
double l_Lip = 1.0;
double l_max=1.0;
int saveSubLevel=1;

double level=1000000.0;
int sortieOKinf[dim]={0};
int sortieOKsup[dim]={0};
int sortieOK[2*dim]={0,1, 0, 1, 0, 1, 0, 1};
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
   int compute_tmin=0;
#endif /* EX3BIS_VIABI2D_H_ */
