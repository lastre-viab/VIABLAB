/*
 * ViabiMicroMacroDiscrete.h
 *
 *  Created on: 1 ao�t 2024
 *      Author: adesi
 */

#ifndef VIABIMICROMACRODISCRETE_H_
#define VIABIMICROMACRODISCRETE_H_

#include "Viabi.h"
//#include "defs.h"
#include "GridMicroMacro.h"
#include "ParametersManager.h"
#include "ViabiMicroMacroTrajectoryHelper.h"

class ViabiMicroMacroDiscrete: public Viabi
    {
public:
    ViabiMicroMacroDiscrete(ParametersManager *pm);

    virtual ~ViabiMicroMacroDiscrete();
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

    virtual void ViabilityKernel(bool sortieOK, int nbArret);
    virtual void CaptureBasin();
    virtual void GarantedViabilityKernel(bool sortieOK, int nbArret);

    virtual void initialiseConstraints();
    virtual SysDyn* GetSysDynForViabProblem();
    GridMicroMacro* GetGridForViabProblem();
    void computeViableTrajectories();

    virtual void computeTrajectories();
    virtual void loadViableSets();
    virtual void saveViableSets();

    void viabKerValFunc(unsigned long long int nbArret);
    void viabKerGarantiValFunc(unsigned long long int nbArret);

    double computeOptimalCaptTrajectory(double *initPosition, string fileName,
	    bool &succes);


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
    void initialiseTargetHJB();
    void initCnIndices();

    /*!
     * Nombre de threads OMP disponible pour une parallélisation OPENMP
     */
    int nbOMPthreads;

    double computeOptimalTrajectory(double *initPosition, string FileName);
    double computeOptimalTrajectory(unsigned long long int posX,
	    string FileName);
    void computeEvaderCapturability(double **trajectory, double *win,
	    double *departureTimes, unsigned long long int *capturePoints,
	    int nbTrajPoints, string fileName);
    void computeEvaderCapturability_proj(double **trajectory, double *win,
	    double *departureTimes, unsigned long long int *capturePoints,
	    int nbTrajPoints, string fileName, unsigned long long int *proj);

    void computeAllPossibleCaptureTrajectories(double *departureTimes,
	    unsigned long long int *capturePoints, int nbTrajPoints,
	    string prefix);

private:

    ViabiMicroMacroTrajectoryHelper *trajectoryHelper;

    GridMicroMacro *grid;
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

    double *vTab;
    double *vTab_tmp;

    void InitViabiMicroMacroDiscrete(algoViabiParams avp);
    void saveValFunctions();
    void computeOptimalTrajectories();

    vector<unsigned long long int> indicesCn;
    vector<unsigned long long int> indicesCn_tmp;

    /* ******************************************************************************************************
     * Attributs  servant de variables globales pour les différentes méthodes
     *  de la classe. A utiliser avec prudence  car toutes les méthodes y ont accès
     *  en écriture!  Attention en cas de parallélisation!
     *********************************************************************************************************/

    unsigned long long int *intPointCoords, *intVect1;
    double *doublePointCoords, *doubleVect, *doubleVect1;
    unsigned long long int *intControlCoords;
    double *doubleControlCoords;
    unsigned long long int *imageCells;

    string filePrefix;

    void (ViabiMicroMacroDiscrete::*computeCurrentImage)(int);

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

    void computeCurrIm(int iter);
    void computeDiscreteImageOfPoint(unsigned long long int num);

    void addDataToCurrentImage(list<imageCell>::iterator *startIt,
	    imageCell newCell, list<imageCell>::iterator *resIt);

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

    int (ViabiMicroMacroDiscrete::*addNewPointsToSet)();

    /*!
     * \brief Cette fonction transforme la liste de mailles représentant une image
     * \f$  \Phi(C_{n+1}\setminus C_n)\f$  en une liste de points
     *
     * Le passage de la liste de mailles à la liste de points se fait en éliminant les doubons (car un point peut appartenir
     * à plusieurs mailles de l'image)  et en calculant la valeur minimal
     * dans le point, comme minimum  des valeurs de toutes les mailles auxquelles il appartient.
     */

    void createPointsList();

    void (ViabiMicroMacroDiscrete::*createCurrentPointsList)();

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

    void (ViabiMicroMacroDiscrete::*addDataToCurrentCell)(list<imageCell>::iterator,
	    imageCell);

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
    void addDataToPoint(list<imagePoint>::iterator itPoint,
	    imagePoint newPoint);

    void (ViabiMicroMacroDiscrete::*addDataToCurrentPoint)(list<imagePoint>::iterator,
	    imagePoint);
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

    void addTempDataToPointsList(list<imagePoint>::iterator *startIt,
	    imageTempPoint newPoint, list<imagePoint>::iterator *resIt);
    void addDataToGivenPointsList(imagePointsList *tempImagePointsList,
	    list<imagePoint>::iterator *startIt, imagePoint newPoint,
	    list<imagePoint>::iterator *resIt);

    void addDataToPointsList(list<imagePoint>::iterator *startIt,
	    imagePoint newPoint, list<imagePoint>::iterator *resIt);


    /*!
     *  \brief  Une structure servant à stocker les données  de l'image discréte
     *  d'un point au cours de parcours  de calcul  de \f$ \Phi(C_{n}\setminus C_{n-1}) \f$ .
     *  \see computeDiscreteImageOfPoint
     *  \see computeCurrentImageGlobalRho
     *
     *  Attention! Variable globale dans les méthodes de la classe!
     */
    discretImageSet_DD pointDI;

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

    std::list<imagePoint> *tempPointsList1;
    std::list<imagePoint> *tempPointsList2;
    int whichPointListToUse;


    bool testConstraintesForCell(unsigned long long int numCell);

    int targ_or_dep;
    int computeTmin;

    };


#endif /* VIABIMICROMACRODISCRETE_H_ */
