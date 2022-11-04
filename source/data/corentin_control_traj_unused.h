/*! \file  Corentin4D_unused.h
 *
*  Ce fichier regroupe les paramètres qui doivent être définis mais qui ne
*  seront pas utilisés dans l'exemple considéré
 *
 *
 */


#ifndef CORENTIN4D_UNUSED_H_
#define CORENTIN4D_UNUSED_H_






/*!
 * \var rafine
 * \brief Paramètre qui définit le nmbre de raffinements successifs
 * à réaliser en partant de la grille initiale définie dasn le modèle
 *
 * Si ce paramètre est égal à zéro, un seul calcul sera réalisé
 * avec le maillage défini dans le modèle.
 *
 */
int refine=0;
/*!
 * \var saveCoupeBound
 * \brief Paramètre qui indique si on souhaite sauvegarder
 * la frontière de l'ensemble calculé après la fin d'un calcul
 *
 * Si ce paramètre est égal à 1 et si refine>0 alors la frontière sera sauvegardée après
 * chaque calcul  ( sinon, l'utilisateur doit modifier le programme principal standard pour
 * contrôler plus précisément les modalités de sauvegarde).
 */
int saveCoupeBound=1;

/*!
 * \var saveCoupe
 *
 * \brief Paramètre qui indique si on souhaite sauvegarder une coupe de l'ensemble
 * calculé.
 *
 * Si ce paramètre est égal à 1, l'utilisateur doit indiquer les
 * variables qu'il faut fixer et les valeurs à donner aux variables fixées
 */
int saveCoupe=0;
/*!
 * \var sliceDirs
 * \brief Tableau de dimension égale à la dimension d'état du système
 * qui indique quelles variables doivent être fixées dans l'enregistrement d'une coupe.
 *
 * Si sliceDirs[i]=1, la variable correspondante sera fixée : x[i]=slicevals[i].
 * La coupe enregistrée regroupera tous  les points du noyau de viabilité tels que
 * pour tout i tel que sliceDirs[i]=1, on a x[i]=slicevals[i]
 */
int sliceDirs[dim]={0, 0, 1, 0};
/*!
 * \var sliceVals[dim]
 * \brief Tableau qui indique les valeurs pour les variables qui doivent être fixées
 * lors d'un enregistrement de coupe.
 */
double sliceVals[dim]={0., 0, 1.0, 0.0};

unsigned long long int trajProjection[dim]={0,0,0,0};



int scaling[dim]={0,0,0,0};



double level=T;
double l_Lip = 1.0;
double l_max=1.0;
int compute_tmin= 0;

int sortieOKinf[dim]={0, 0, 0, 0};
int sortieOKsup[dim]={0,0, 0, 0};

int saveSubLevel=0;



int compteOnlyOptimalRetro=1;


int computeSet=1;


void loadModelData()
{

}


void dynamics_continuous(double * x, double *u, double * image);
void dynamics_discrete(double * x, double *u, double * image);


/*!
 * jacobian matrix of the dynamics
 */

inline void jacobian(double *x, double *u , double ** jacob)
{
  jacob[0][0]=0.0;
  jacob[0][1]=0.0;
  jacob[0][2]=0.0;
  jacob[0][3]=0.0;

  jacob[1][0]=0.0;
  jacob[1][1]=0.0;
  jacob[1][2]=0.0;
  jacob[1][3]=0.0;

  jacob[2][0]=0.0;
  jacob[2][1]=0.0;
  jacob[2][2]=0.0;
  jacob[2][3]=0.0;

  jacob[3][0]= sign(x[0]-sin(x[2]) )*( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) )+
             ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( sign(x[0]+sin(x[2]) )  );

  jacob[3][1]=( sign( x[1]-cos(x[2]) ) )* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) )+
              ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( sign( x[1]+cos(x[2]) ) );

  jacob[3][2]=( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* (sign(x[0]+sin(x[2]) )*cos(x[2])- sign( x[1]+cos(x[2]) )*sin(x[2]) )+
              ( -sign(x[0]-sin(x[2]) )*cos(x[2]) + sign( x[1]-cos(x[2]) ) * sin(x[2]))* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) );

  jacob[3][3]=0.0;

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
  res[0]=0.;
  res[1]=0.;
  res[2]=0.;
  res[3]= ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) );

}


void dynamics_tych_fd(unsigned long long int  * x, unsigned long long int *u,unsigned long long int *v, unsigned long long int * image)
      {

      }
inline double  constraintsX_fd( unsigned long long int * x )
{

	double res = 1.0;
	return res;
}

void dynamics_fd(unsigned long long int  * x, unsigned long long int *u, unsigned long long int * image)
{
}


inline double  constraintsXU_fd( unsigned long long int * x, unsigned long long int * u )
{
	return 1.0;
}



/*!
 * \brief Function  defining the state   constraints, corresponds  to k(x)
 *
 * This function defines the set K for admissible  states
 * @param x state variable
 * @return  value that characterize the constraints set
 */
inline double constraintsX( double * x )
{


	return 1.0;
}


/*!
 * Function that  characterise  the target set C, corresponds to c(x)
 * @param x state variable
 * @return value  to charaterise the target set
 */
inline  double target (double * x)
{
  double res =PLUS_INF;

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
   void loadModelData();

   unsigned int maxNbRetro=1048;
#endif /* TESTDATA_H_ */
