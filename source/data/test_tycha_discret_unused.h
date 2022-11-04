/*! \file  testZermelo_unused.h
 *

 *
 *
 */


#ifndef TESTZERMELO_UNUSED_H_
#define TESTZERMELO_UNUSED_H_
// double a=-0.202;
// double b =-0.787;

//double a=0.1;
//double b =0.613;

// double a=0.2725;
//double b=-0.575;

// double a=-0.361;   // julia 1
//  double b=0.63725;

int saveCoupeBound=0;
int saveCoupe=0;
int sliceDirs[dim]={  1, 0};
double sliceVals[dim]={  1.0, 0.0};
int saveProjection=0;
/*
 * Paramètre permettant d'indiquer que la borne sur du pavé de calcul par rapport
 * n'ets pas une contrainte "naturelle" du problème; Ainsi la sortie du domaine
 *  de calcul dans cette direction ne sera pas considéré comme non viable
 *  Permet de géré le fait que le domaine de certaines variables est infini
 */

/*
int sortieOKinf[dim]={0, 0};
int sortieOKsup[dim]={0,0};
int compute_tmin= 0;
*/

int saveSubLevel=0;
double level=T;
double l_Lip = 1.0;
double l_max=1.0;
int compute_tmin= 0;


int ompThreads=1;

/*
 * sélection de type de reconstruction de trajectoire
 *  VD= 1, viable par défaut, on sélectionne le premier contrôle viable trouvé
 *  VL= 2, viable lourd: le contrôle reste constant tant qu'il viable;
 *  cette méthode nécessite une initialisation de contrôle
 */
int typeTraj=VD;

double initControls[dimc*nbTrajs]={};

/*!
 * nature du probleme
 */
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
  * Function that defines the  dynamical system
  */

   void dynamics_hybrid(double * x, double *u, double * image)
 {

 }


   unsigned long long int projection[dim]={0};

   unsigned long long int trajProjection[dim]={0};
   double initPoint[dim]={0.5};

   int scaling[dim]={0};

   /*
    * target = 1
    * departure =0;
    * Ce parametre determine le sens des trajectoires
    */
   int target_or_departure_problem=1;

   /*!
    * \var globalDeltaT
    *  bool�en indique si le pas de temps  doit �re choisi globalement
    *  ou localement pour les algorithmes de viabilit�
    */
   bool globalDeltaT=false;


   int compteOnlyOptimalRetro=1;

   unsigned int maxNbRetro=1048;

/*
   void loadModelData()
   {

   }
*/
   /*!
    * jacobian matrix of the dynamics
    */
   inline void jacobian(double *x, double *u , double ** jacob)
   {
     jacob[0][0]=0;
     jacob[0][1]=0;
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
   double LC= 1.0;
   /*!
    * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
    * @param[in] x  the state variable
    * @param[out] res  the result
    */
   inline void localDynBounds(double * x, double * res)
   {
     res[0]=  1.0;
     res[1]=  1;
   }
   /*!
    * \var M : global bound  for the dynamics if defined
    */
   double M= 1.0;
   /*!
    * \var computeM : indicates which method  to use to compute the bound for dynamics
    * 0= analytic  global value used
    * 1= local calculation using the localDynBounds function
    * 2= local calculation using  explicit maximization on controls
    */
   const int computeM=1;
   
   /*!
    * Function that  characterise  the target set C, corresponds to c(x)
    * @param x state variable
    * @return value  to charaterise the target set
    */
   inline  double target (double * x)
   {
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
