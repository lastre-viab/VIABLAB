/*! \file Zermelo_Lmin_unused.h
 *
*  Ce fichier regroupe les paramètres qui doivent être définis mais qui ne
*  seront pas utilisés dans l'exemple considéré
 *
 *
 */


#ifndef EQUILIBRES4D_UNUSED_H_
#define EQUILIBRES4D_UNUSED_H_


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

unsigned long long int trajProjection[dim]={0,0};


double initControls[dimc*nbTrajs]={};

/*
 * sélection de type de reconstruction de trajectoire
 *  VD= 1, viable par défaut, on sélectionne le premier contrôle viable trouvé
 *  VL= 2, viable lourd: le contrôle reste constant tant qu'il viable;
 *  cette méthode nécessite une initialisation de contrôle
 */
int typeTraj=VL;


int scaling[dim]={0,0};


/*
 * target = 1
 * departure =0;
 * Ce parametre determine le sens des trajectoires
 */
int target_or_departure_problem=1;

int compteOnlyOptimalRetro=1;


void loadModelData()
{

}


void dynamics_continuous(double * x, double *u, double * image);
void dynamics_discrete(double * x, double *u, double * image);

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

  return 1.0;
}



/* *   *************************************************************************
 *      Definition of value function
 ***************************************************************************/






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
