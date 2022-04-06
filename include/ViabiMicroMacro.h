/*
 * ViabiHJB.h
 *  *
 *    VIABLAB : a numerical library for Mathematical Viability Computations
 *    Copyright (C) <2020>  <Anna DESILLES, LASTRE>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *   
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Created on: 9 déc. 2013
 *      Author: ANYA
 */

/*!
 * \page introViabiHJB ViabiHJB:  Librairie d'algorithmes de viabilité épigraphiques
 *
 *
 * \section p1s1 Introduction
 * La classe ViabiHJB implémente des algorithmes  de viabilité
 * dédiés  aux calculs des ensembles  qui peuvent être représentés comme
 *  épigraphes de fonctions. Il peut s'agir de fonctions valeur de problèmes
 *  de contrôle optimal ou  de solution d'EDP de type HJB.
 *
 *  Ces algorithmes  fonctionnent en laissant une dimension, celle des valeurs de la fonction, libre
 *   et non discrétisée grâce au stockage   particulier des données, réalisé par la classe correspondante GridHJB. Selon la nature de la fonction
 *    dont on souhaite calculer l'épigraphe on utilisera un algorithme approprié de la classe ViabiHJB.
 *
 *    Chaque algorithme calcule une fonction valeur (ou son épigraphe) et la (ou les) fonction(s) de rétro-action associée(s).
 *    Chaque fonction valeur, quelle quoi soit la nature du problème de contrôle optimal associé,  sera stockée
 *    dans une table (base de données) dont les éléments sont  des paires  de la forme (clé, donnée).
 *    La clé, comme la donnée peuvent prendre des formes
 *    différentes, selon le problème. La clé sert à indexer les données :
 *     tout accès aux données se fait à partir des clés correspondantes.
 *
 *     Dans le cas de fonction valeur la clé sera de deux types:
 *     	- si la fonction valeur \f$ V(x)\f$ ne dépend pas de temps ( correspond aux problèmes
 *     	de contrôle optimal à horizon infini, ou problèmes de temps minimum, par exemple) la clé  sera
 *     	l'identifiant de \f$ x\f$ : numéro du point dans la grille. Ainsi dans ce cas les éléments de la table représentant
 *     	la fonction valeur  seront de la forme \f$ (x,V(x))\f$ .
 *
 *     	- si la fonction valeur  \f$ V(x,t)\f$ dépend de temps (problèmes à horizon fini) la clé  sera elle méme un couple
 *     	: \f$ (x,t)\f$ oé  \f$ x\f$ est représenté par son numéro dans la grille et $\f$ t \f$  est une valeur réelle.
 *
 *
 *   \section p1s2 Algorithmes de bassin de capture épigraphiques
 *  Le bassin de capture \f$  Capt_{F}(K,C)\f$   peut être  défini comme le
 *  domaine de la fonction temps minimal
 *  \f$  \omega(x)=\inf\{t\le 0,\ \ x(t)\in C, \ x(0)=x,\ \ \forall s\in[0,t],\ x(s)\in K\} \f$  .
 *  En utilisant cette définition, il suffit de calculer la fonction temps minimal pour déduire
 *  le bassin de capture d'une cible donnée \f$ C\f$  par un systéme dynamique \f$  F\f$  sous les contraintes
 *  \f$  x\in K \f$.
 *
 *   L'algorithme implémenté dans classe ViabiHJB pour le calcul des bassins de capture
 *   calcule l'épigraphe  de la fonction temps minimum. Deux versions de l'algorithme  sont proposées.
 *
 *  La fonction minTimeGlobalRho() est une version oé le pas de discrétisation
 *   temporelle \f$  \rho\f$   est  défini de façon globale par rapport l'état \f$  x \f$  .
 *   Il est déterminé  en fonction du schéma
 *    de discrétisation  en fonction des constantes  de régularité de la dynamique :
 *    \f$  L \f$  , la constante de Lipschitz  et \f$  M= sup \| F(x)\| \f$  .
 *
 *   La fonction minTimeLocalRho() est une version oé le pas de discrétisation temporelle \f$ \rho\f$   est calculé localement pour chaque point
 *    \f$  x\f$  en fonction du  schéma  de discrétisation  est des estimations locales  des constantes  de régularité de la dynamique :
 *    \f$  L\f$ , la constante de Lipschitz  et \f$  M= sup \| F(x)\| \f$.
 *
 *
 */


#ifndef VIABIHJB_H_
#define VIABIHJB_H_

#include "Viabi.h"
//#include "defs.h"
#include "GridMicroMacro.h"
#include "ParametersManager.h"


class ViabiMicroMacro: public Viabi {
public:
	ViabiMicroMacro(ParametersManager *pm);

	virtual ~ViabiMicroMacro();
	/*!
	 *  \brief Méthode de débuggage:  sert à afficher quelques informations  sur la classe
	 *  dans la console
	 */
	virtual void printViabiInfo();
	/*!
	 * \brief  Méthode permettant d'initialiser l'ensemble cible
	 *
	 * Cette classe pour l'instant  utilisera une méthode propriétaire
	 * \see initialiseTargetHJB()
	 */
	virtual void initialiseTarget();
	//virtual void initialiseConstraints() const;
	//virtual void initialiseTargetOrConstraints() const;

	virtual void ViabilityKernel( bool sortieOK,int nbArret);
	virtual void CaptureBasin( );
	virtual void GarantedViabilityKernel( bool sortieOK,int nbArret);

	virtual	void initialiseConstraints();
		void initialiseConstraints_CC();
		void initialiseConstraints_DD();
		void computeViableTrajectories();
		double computeViableTrajectory(double  *initPosition,  double initValue, string fileName, bool &succes);
		double computeViableTrajectory_DD(unsigned long long int *initPosition,  double initValue, string fileName, bool &succes);
		unsigned long long int  findViabControl_DD(double budget, unsigned long long int *currentPos,
				unsigned long long int *resPos,
				double & newBudget,
				bool &succes );
	virtual void computeTrajectories();
	virtual void loadViableSets();

	void viabKerValFunc();
	void viabKerValFunc_DD();

	double computeOptimalCaptTrajectory(double *initPosition, string fileName, bool &succes);


	/*!
	 *  \brief Cette foncion calcule l'épigraphe  de la fonction temps minimal
	 * dans l'hypothése que le pas de temps  est ajusté localement en fonction
	 *  d'une estimation  de la constante de Lipschitz et de M
	 *
	 * Le schéma de discrétisation est un schéma implicite, d'Euler ou RK ou autre
	 *   de la forme
	 *   \f[
	 *   x_{n+1}\in G_{imp}(x_n),\ \ G_{imp}(z)=\{y\in X,\ \ z\in \Phi_{imp}(y)\}
	 *   \f]
	 *
	 *   pour le schéma implicite d'Euler \f$\Phi_{imp}(y)=y-\rho F(y)\f$
	 *
	 *   pour le schéma RK2 : \f$\Phi_{imp}(y)=y-0.5 \rho (F(y)+F(y-\rho F(y))\f$
	 */


	void  captBasinEpi_omp();




	/*!
	 *   \brief  Fonction qui initialise la cible dans la base de données représentant
	 *   la fonction valeur à calculer
	 *
	 *   \todo Cette  fonction actuellement réalise l'initialisation d'un attribut de la
	 *   clsse : la liste des points de l'image en cours. Cette liste pour la toute première itération correspond tout simplement
	 *   à la liste des points de la cible.  Ce détail technique ne permet pas que cette fonction soit virtuelle,
	 *   car elle tocuhe à un  attribut spécifique à la classe ViabiHJB.  A modifier:  laisser la fonction virtuelle d'intialisation
	 *    de la cible telle que programmée, sans initialiser la liste  des points. Ensuit faire une
	 *    méthode privée  qui initialise cette liste ainsi que les bases de données des rétro actions, qui sont également
	 *      des attributs  spécifique. Cette méthode privée sera appelée par les méthodes
	 *      réalisant les différents algorithmes épigraphiques
	 *       au tout début avant la première itération
	 *
	 */
	void initialiseTargetHJB() ;
	void initialiseTargetHJB_DD() ;
	void initCnIndices(  );

	/*!
	 * Nombre de threads OMP disponible pour une parallélisation OPENMP
	 */
	int nbOMPthreads;

	double computeOptimalTrajectory(double * initPosition,string FileName);
	double computeOptimalTrajectory(unsigned long long int posX,string FileName);
	void computeEvaderCapturability(double ** trajectory, double * win, double * departureTimes,  unsigned long long int* capturePoints,int nbTrajPoints,string fileName);
	void computeEvaderCapturability_proj(double ** trajectory, double * win, double * departureTimes,  unsigned long long int* capturePoints,int nbTrajPoints,string fileName,unsigned long long int * proj);

	void computeAllPossibleCaptureTrajectories(  double * departureTimes, unsigned long long int* capturePoints, int  nbTrajPoints, string prefix );

private:

	GridMicroMacro * grid;
	/*!
	 *  \brief  Copie pour raisons de rapidté  de la valeur de dimension d'état
	 */
	int dim;
	/*!
	 *  \brief  Copie pour raisons de rapidté  de la valeur de dimension de contrôle
	 */
	int dimC;
	/*!
	 *  \brief Pointeur sur la base de données servant à enregister la rétroaction optimale
	 */

	double * vTab;
	double *vTab_tmp;

	void InitViabiMicroMacro(algoViabiParams avp);
	void saveValFunctions();
	void computeOptimalTrajectories();


	vector<unsigned long long int > indicesCn;
	vector<unsigned long long int > indicesCn_tmp;


	/* ******************************************************************************************************
	 * Attributs  servant de variables globales pour les différentes méthodes
	 *  de la classe. A utiliser avec prudence  car toutes les méthodes y ont accès
	 *  en écriture!  Attention en cas de parallélisation!
	 *********************************************************************************************************/

	unsigned long long int * intPointCoords, *intVect1;
	double * doublePointCoords, *doubleVect, *doubleVect1;
	unsigned long long int * intControlCoords;
	double * doubleControlCoords;
	unsigned long long int * imageCells;

	string filePrefix;

	void (ViabiMicroMacro::*computeCurrentImage)(int);

	/*!
	 *
	 *  \brief  Cette  fonction  réalise le calcul de
	 *  \f$  \Phi(C_{n+1}\setminus C_n)\f$ avec un calcul de rho local.
	 *
	 *  On parcourt  l'ensemble en construction
	 *   pour repérer les points qui appartiennent à \f$ C_n\f$  : pour cela on lit l'indicateur b
	 *    dans l'enregistrement de chaque point
	 *     si b=1  alors c'est un point ajouté à la dernière itération  et donc il appartient à \f$ C_n\f$
	 */
	void computeCurrIm_tmin(int iter);

	/*!
	 *
	 *  \brief  Cette  fonction  réalise le calcul de
	 *  \f$  \Phi(C_{n+1}\setminus C_n)\f$ avec un calcul de rho local.
	 *
	 *  On parcourt  l'ensemble en construction
	 *   pour repérer les points qui appartiennent à \f$ C_n\f$  : pour cela on lit l'indicateur b
	 *    dans l'enregistrement de chaque point
	 *     si b=1  alors c'est un point ajouté à la dernière itération  et donc il appartient à \f$ C_n\f$
	 */
	void computeCurrIm_Lmin(int iter);

	void computeCurrIm_DD( int iter);
	void computeDiscreteImageOfPoint_DD(unsigned long long int num);

	void addDataToCurrentImage(list<imageCell>::iterator *startIt, imageCell newCell,list<imageCell>::iterator *resIt );



	/*!
	 * \brief Cette fonction parcours la liste de points générée lors d'un calcul de l'image \f$  \Phi(C_{n+1}\setminus C_n)\f$
	 * et enregistre leurs données dans la base.
	 *
	 * Cette fonction intervient dans les méthodes
	 *  de calcul de la fonction temps minimal à la fin d'une étape de calcul.  Une fois  qu'une liste ordonnées et sans doublons de
	 *  tous le spoints de l'image \f$  \Phi(C_{n+1}\setminus C_n)\f$
	 *  est calculée, cette fonction va enregistrer  les données associées aux points dans différentes bases.
	 *  Elle enregistre les paires (x,V(x))  dans la base de la fonction valeur  (la base de la classe grid)
	 *  et elle enregistre les données de rétro-actions dans les bases correspondantes, associées à cette classe.
	 *
	 *   Tout  en parcourant la liste de points, cette fonction la modifie :
	 *
	 *    - si le point existe déjà dans la base mais que sa valeur est inférieure à celle trouvée dan l'image en cours,
	 *    le point est retiré de la liste qui sera ensuite utilisée pour le calcul de l'image suivante
	 *
	 *    -si le point n'existe pas ou si sa valeur dans la base est midifiée (parce que la valeur du point calculé est inférieure)
	 *    le point est maintenu dans la liste ets era utilisé pour le calcul de l'image suivante.
	 *
	 * @return nombre de points nouveaux : soit parce qu'ils n'existaient pas dans la base de données
	 * soit parce que leur fonction valeur a été modifiée par une plus petite.
	 */

	int addNewPoints();



	int (ViabiMicroMacro::*addNewPointsToSet)();


	/*!
	 * \brief Cette fonction transforme la liste de mailles représentant une image
	 * \f$  \Phi(C_{n+1}\setminus C_n)\f$  en une liste de points
	 *
	 * Le passage de la liste de mailles à la liste de points se fait en éliminant les doubons (car un point peut appartenir
	 * à plusieurs mailles de l'image)  et en calculant la valeur minimal
	 * dans le point, comme minimum  des valeurs de toutes les mailles auxquelles il appartient.
	 */

	void createPointsList();
	void createPointsList_DD();



	void (ViabiMicroMacro::*createCurrentPointsList)();



	/*!
	 * Méthode de debuggage: permet d'afficher sur la console
	 * les données des mailles  de la liste représentant  l'image en cours
	 * \see currentImageList
	 */
	void showCurrentImageList();

	/*!
	 *  \brief Méthode de debuggage: permet d'afficher sur la console
	 * les données des points  de la liste représentant  l'image en cours
	 * \see currentImagePointsList
	 */
	void showCurrentImagePointsList();

	/*!
	 *  \brief Cette fonction permet de mettre à jour les données d'une maille déjà présente dans
	 * l'image \f$ \Phi(C_{n}\setminus C_{n-1})\f$ en construction lorsque cette maille
	 * est à nouveau touchée par une évolution
	 *
	 *
	 * Les données mises à jour sont:
	 *  - la valeur minimale
	 *  - rétroaction optimale
	 *  - rétroaction viable
	 *
	 * @param itCell pointeur sur la maille à modifier dans la liste de mailles en construction
	 * @param newCell structure portant les données à ajouter à cette maille
	 */
	void addDataToCell(list<imageCell>::iterator itCell, imageCell newCell);



	void (ViabiMicroMacro::*addDataToCurrentCell)(list<imageCell>::iterator , imageCell );

	/*!
	 * \brief Cette fonction permet de mettre à jour les données d'un point déjà présent dans
	 * l'image \f$ \Phi(C_{n}\setminus C_{n-1})\f$ en construction lorsque ce point appartient
	 * à plusieurs mailles présetes dans l'image
	 *
	 *
	 * Les données mises à jour sont:
	 *  - la valeur minimale
	 *  - rétroaction optimale
	 *  - rétroaction viable
	 *
	 * @param itPoint pointeur sur le point à modifier dans la liste de points en construction
	 * @param newPoint structure portant les données à ajouter à cette maille
	 */
	void addDataToPoint(list<imagePoint>::iterator itPoint, imagePoint newPoint);

	void (ViabiMicroMacro::*addDataToCurrentPoint)(list<imagePoint>::iterator, imagePoint );
	/*!
	 *  \brief Cette fonction permet d'ajouter un nouveau point dans la liste ordonnée et sans doublons de points de l'image
	 *  \f$ \Phi(C_{n}\setminus C_{n-1})\f$ .
	 *
	 *  Si le point n'existe pas il est ajouté dans l'ordre croissant, sinon
	 *  ses données sont regroupées avec les données existantes.
	 * @param startIt pointeur sur le dernier point inserré dans la liste, permet de limiter la complexité de recherche
	 * @param newPoint le point à inserrer
	 * @param resIt pointeur sur le nouveau point inserré, sera utilisé pour
	 * la recherche suivante
	 */

	void  addTempDataToPointsList(list<imagePoint>::iterator *startIt, imageTempPoint newPoint,list<imagePoint>::iterator *resIt );
	void addDataToGivenPointsList(imagePointsList * tempImagePointsList, list<imagePoint>::iterator *startIt, imagePoint newPoint,list<imagePoint>::iterator *resIt );

	void addDataToPointsList(list<imagePoint>::iterator *startIt, imagePoint newPoint,list<imagePoint>::iterator *resIt );

	/*!
	 *  \brief  Une structure servant à stocker les données  de l'image discréte
	 *  d'un point au cours de parcours  de calcul  de \f$ \Phi(C_{n}\setminus C_{n-1}) \f$ .
	 *  \see computeDiscreteImageOfPoint
	 *  \see computeCurrentImageGlobalRho
	 *
	 *  Attention! Variable globale dans les méthodes de la classe!
	 */
	discretImageSet pointDI;

	/*!
	 *  \brief  Une structure servant à stocker les données  de l'image discréte
	 *  d'un point au cours de parcours  de calcul  de \f$ \Phi(C_{n}\setminus C_{n-1}) \f$ .
	 *  \see computeDiscreteImageOfPoint
	 *  \see computeCurrentImageGlobalRho
	 *
	 *  Attention! Variable globale dans les méthodes de la classe!
	 */
	discretImageSet_DD pointDI_DD;

	/*!
	 *  \brief La liste de structures représentant  des mailles d'une image en construction
	 * Attention! Les méthodes suivantes modifient  cette variable comme étant globale
	 * \see addDataToCurrentImage
	 * \see computeCurrentImageGlobalRho
	 */
	imageCellsList currentImageList;
	/*!
	 * \brief La liste de structures représentant  des points d'une image en construction
	 * Attention! Les méthodes suivantes modifient  cette variable comme étant globale
	 * \see createPointsListGlobalRho
	 * \see addDataToPointsList
	 */
	imagePointsList currentImagePointsList;

	std::list<imagePoint> * tempPointsList1;
	std::list<imagePoint> * tempPointsList2;
	int whichPointListToUse;

	/*!
	 * \brief Fonction qui calcule l'image discrére d'un point \f$ \Phi(x)\f$
	 * @param num numéro du point
	 */

	void  computeDiscreteImageOfPoint(unsigned long long int num);


	/*!
	 * \brief  Version parallélisée OpenMP  de la fonction qui calcule l'image discrére d'un point \f$ \Phi(x)\f$
	 * La parallélisation est faite ici en distribuant la boucle qui parcourt
	 * la liste de tous les controles.
	 *
	 * \todo Cette version est expérimentale, à tester. LE problème est que pour des
	 * dynamiques simples, la tâche confiée à chaque thread est trop courte par rapport à l'effort de
	 * distribution et de synchronisation; A voir plus tard d'autres façons de paralléliser  le calcul
	 *
	 * @param num numéro du point
	 */

	/*!
	 * \brief Fonction qui calcule la première itération d'un algorithme
	 *  de type bassin de capture "direct", issu d'une discrétisation implicite
	 *
	 *  Pour le calcul d'un bassin de capture on suppose que la dynaique contient zéro
	 *  sur la cible : autrement dit, on peut s'arréter, une fois arrivé sur la cible.
	 *  Pour l'application de  tous les théorémes de viabilité dans ce cas
	 *  nous devons avoir une dynamique à image convexe. Donc sur la cible, on définit la dynamique comme
	 *  \f$ \overline{Co}(\{0\} \cup F(x) )\f$. Comme dans es algorithmes issus de la discrétisation implicite
	 *  la première itération calcule l'image \f$ \Phi(C) \f$ on doit donc tenir compte que la
	 *  dynamique est "convexifiée" sur la cible. Cette fonction réalise  le calcul  de \f$ \Phi(C) \f$
	 *  en tenant compte de la convexification.
	 *
	 *  Le principe de calcul est le méme que pour la fonction qui calcule l'image  \f$ \Phi(C_{n}\setminus C_{n-1}) \f$ :
	 *  pour chaque point  de \f$ x\in C \f$ on calcule d'abord \f$ \Phi(x) \f$  sous forme de liste de mailles.
	 *  Ensuite pour chaque maille \f$ m_0\f$  d'une telle liste on ajoute dans l'image en construction la maille elle méme
	 *  avec sa valeur ainsi que toutes les mailles qui sont croisées par le segment reliant le point de départ \f$ x\f$
	 *  à \f$ m_0\f$ . La valeur associée à é chacune de ces ailles intermédiaires est une fraction du pas \f$ \rho(x)\f$
	 *  utilisé pour le calcul de l'image  \f$ \Phi(x) \f$ .
	 *
	 * @param iter :  numéro d'itération
	 *
	 * \see minTimeEpiLocalRho
	 */

	void  computeConvexifiedImage_tmin( int iter);
	void  computeConvexifiedImage_Lmin( int iter);
	void  computeConvexifiedImage_DD( int iter);

	void (ViabiMicroMacro::* computeFirstConvexifiedImage)( int iter);




	void  computeConvexifiedImage_tmin_omp( int iter);
	void  computeConvexifiedImage_Lmin_omp( int iter);


	void (ViabiMicroMacro::* computeFirstConvexifiedImage_omp)( int iter);





	/*!
	 * \brief Fonction qui permet d'ajouter dans l'image convexifiée en construction
	 * toutes les mailles croisées par un segment reliant un point de départ et une maille de
	 * son image discréte  \f$ \Phi(x) \f$ .
	 *
	 * @param posX numéro du point de départ
	 * @param numCell numéro d'une maille de son image
	 * @param tempImageCell pointeur sur une structure de type imageCell pour récupérer les données
	 * à inserrer dans l'image en construction
	 * @param rho pas de temps utilisé pour le calcul de  \f$ \Phi(x) \f$
	 */
	void addConvexCombinations(list<imagePoint>::iterator itPoint, unsigned long long int numCell, imageCell * tempImageCell,double rho ,list<imageCell>::iterator *itStart);

	bool testConstraintesForCell(unsigned long long int numCell);

	int targ_or_dep;
	int computeTmin;


	double computeOptimalTrajectory_tmin(double *initPosition, string fileName, bool &succes);

	int findOptiControl_tmin(double *currentPos,
			double &dt,
			int nbStepIter,
			double stepCoeff,
			double *resPos,
			int & nbViabVoisins,
			bool &succes );




	double computeOptimalTrajectory_Lmin(double *initPosition, string fileName, bool &succes);

	int findOptiControl_Lmin(double budget, double *currentPos,
			double &dt,
			int nbStepIter,
			double stepCoeff,
			double *resPos,
			double & newBudget,
			bool &succes );


};

#endif /* VIABIHJB_H_ */
