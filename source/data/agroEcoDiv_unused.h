/*! \file Zermelo_tmin_unused.h
 *
*  Ce fichier regroupe les paramètres qui doivent être définis mais qui ne
*  seront pas utilisés dans l'exemple considéré
 *
 *
 */


#ifndef AGROECODIV_UNUSED_H_
#define AGROECODIV_UNUSED_H_


/*!
 * \var sortieOK[dim]
 * \brief Tableau qui définit les modalités d'interprétation de
 * la sortie du domaine de calcul lors du calcul de noyau de viabilité. Permet d'indiquer
 * si le domaine réel est infini le long d'une direction , alors que le domaine de calcul est forcément borné.
 * Si la valeur correspondante est 1 alors  si le successeur d'un point x a  la composante correspondante au-delà
 * des bornes du domaine de calcul il sera considéré comme viable.
 */
int sortieOK[dim]={0,0};

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
int saveCoupeBound=0;

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
int sliceDirs[dim]={0, 0};
/*!
 * \var sliceVals[dim]
 * \brief Tableau qui indique les valeurs pour les variables qui doivent être fixées
 * lors d'un enregistrement de coupe.
 */
double sliceVals[dim]={0., 0};



int saveProjection=0;
double *initControls;

int sortieOKinf[dim]={0};
int sortieOKsup[dim]={0};
int intermediate_savings = 0;

int saveBoundary=0;
int saveSubLevel=1;
int gridMethod=MM;

double level=1000000.0;
/*
 * Sélection de l'ensemble à calculer
 */





double *initPoints;
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
  image[2] = 1.0;
}




void dynamics_continuous(double * x, double *u, double * image);
void dynamics_discrete(double * x, double *u, double * image);





/* *   *************************************************************************
 *      Definition of value function
 ***************************************************************************/




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


   /*!
    * jacobian matrix of the dynamics
    */
   inline void jacobian(double *x, double *u , double ** jacob)
   {
   }
   /*!
    * \var computeLC : indicates which method  to use to copute Lipschitz constant
    * 0= analytic  global value used
    * 1= local calculation using the jacobian matrix
    * 2= local calculation using the finite differences
    */
   int computeLC=1;

   /*!
    * \var LC Lipschitz constant if known analytically
    */
   double LC=  max(abs(CONTROL_MIN[0]), abs(CONTROL_MAX[0]));
   /*!
    * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
    * @param[in] x  the state variable
    * @param[../OUT] res  the result
    */
   inline void localDynBounds(double * x, double * res)
   {
   }
   /*!
    * \var M : global bound  for the dynamics if defined
    */
   double M= 0;
   /*!
    * \var computeM : indicates which method  to use to compute the bound for dynamics
    * 0= analytic  global value used
    * 1= local calculation using the localDynBounds function
    * 2= local calculation using  explicit maximization on controls
    */
   int computeM=1;
   /*!
    * Function  defining the  switch conditions  between the continueos and  discrete dynamics
    * @param x state variable
    * @param u control variable
    * @return  value that caraterise the restet  set
    */
   void dynamics_tych_fd(unsigned long long int  * x, unsigned long long int *u,unsigned long long int *v, unsigned long long int * image)
   {

   }
   /*
    * Definition of the dynamics  and associated functions and constants
    */
   void dynamics_tych_fd(unsigned long long int  * x, unsigned long long int *u,unsigned long long int *v, unsigned long long int * image);

   void dynamics_c(double * x, double *u, double * image);
   void dynamics_td(double * x, double *u, double * image);

   void dynamics_c(double * x, double *u, double * image)
   {
   	////cout<< " test dynamique x "<<x[0]<< " phys vect [0]"<<physVect[0]<< " image[0]"<<image[0]<<endl;
   }
   void dynamics_td(double * x, double *u, double * image)
   {
   }

   /*!
    * Function  for optimisation criterion, corresponds to l(x,u)
    * @param x state variable
    * @param u control
    * @return value of l
    */
   inline double l(double * x, double * u )
   {
   	return 0;
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



   /*!
    * Function that  characterise  the target set C, corresponds to c(x)
    * @param x state variable
    * @return value  to charaterise the target set
    */
   inline  double target (double * x)
   {
   	return  0.0;
   }


   inline double l_tych_fd(unsigned long long int  * x, unsigned long long int * u , unsigned long long int * v )
   {
   	return 0.0;
   }


   const int dimc_ty=0;
   double CONTROL_MIN_ty[dimc_ty]={};
   /*! \var CONTROL_MAX_ty[dimc]
    *  \brief maximum values  for  control vector components
    *   \see controlParams
    */
   double CONTROL_MAX_ty[dimc_ty]={};

   unsigned long long int nbPointsControl_ty[dimc_ty] =  {};

#endif /* TESTDATA_H_ */
