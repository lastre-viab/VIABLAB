/*
 * Grid.h
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
 *  Created on: 10 sept. 2013
 *      Author: Anna DESILLES
 */

#ifndef GRID_H_
#define GRID_H_
#include "defs.h"
/*! \class Grid
 *
 * \brief  Classe abstraite pour représenter les grilles
 *
 * La classe abstraite qui définit une grille sans spécifier
 * la méthode de représentation  et de stockage
 * Cette classe regroupe les attributs communs  à toutes les implémetations et
 * un certain nombre de méthodes qui sont communes ou doivent être implémentées
 * par chaque classe qui hérite selon la spécification donnée ici;
 * Une grille représente  toujours un ensemble, par defaut vide
 *
 */
class Grid
{

protected:
    Grid();
public:

    int dim; /*!<\brief Dimension  	 */

    double *limInf; /*!< \brief Valeurs minimales pour X, tableau  */

    double *limSup; /*!< \brief Valeurs maximales pour X, tableau 	 */

    double *step; /*!< \brief  Pas de discretisation pour X , tableau	 */

    double maxStep; /*!< \brief Valeur maximale de pas , scalaire	 */

    bool *periodic; /*!< \brief Indicateur booléen de périodicité, tableau*/
    /*!<
     *  1= variable périodique
     *  0=variable non périodique
     */
    bool arePointsGridCenters;
    int pow3;
    int nbPointsCube;
    unsigned long long int *nbPoints; /*!< \brief Nombre de points de grille par axe ,tableau	 */

    unsigned long long int nbTotalPoints;

    unsigned long long int *vectUnsigIntTemp; /*!< \brief Une mémoire tampon pour les calculs ,tableau	 */
    unsigned long long int *vectInt; /*!< \brief Une mémoire tempon pour les calculs ,tableau	 */
    unsigned long long int getDim() const;

    const long long int* getIndicesDecalCell() const; /*!<\brief accès  aux indices de déclage pour les sommets d'une maille*/

    const long long int *getIndicesDecal() const;
        
    /*!
     * \fn virtual void Grid::printGrid(void) const=0
     *
     * fonction abstraite permettant une impression console des
     * informations de la grille
     * Doit être impéematée par chaque classe qui hérite
     */
    virtual void printGrid(void) const;
    /*!
     *  \brief Destructeur
     *
     *  Chaque classe doit prévoir son propre destructeur en fonction
     *  des choix de gestions des représentations en mémoire
     */
    virtual ~Grid();

    virtual bool isInSet(const unsigned long long int *coords) const;

    /*!
     * \fn virtual void savePointsList(string fileName) const=0
     *
     * fonction abstraite permettant l'enregistrement
     * dans un fichier de la liste des points appartenant
     * à l'ensemble représenté. Les points de la grille
     * qui ne sont pas dans l'ensemble ne sont pas enregistrés
     *
     * \param[in] fileName : chaine de caracteres, le nom du fichier
     * pour l'enregistrement
     * Doit être impéematée par chaque classe qui hérite
     */
    virtual void savePointsList(const string &fileName) const;

    /*!
     * \fn virtual void saveValOnGrid(string fileName) const=0
     *
     * fonction abstraite permettant  d'enregistrer dans un fichier
     * la valeur définissant l'ensemble sur la grille
     * Dans le cas de l'ensemble épigraphique la valeur réelle est enregistrée
     * Dans le cas d'ensemble booléen, la fonction caractéristique est enregistrée
     *
     * \param[in] fileName : chaine de caracteres, le nom du fichier
     * pour l'enregistrement
     * Doit être impéematée par chaque classe qui hérite
     */
    virtual void saveValOnGrid(const string &fileName) const;
    /*!
     * \fn virtual void saveValOnGridLight(string fileName) const=0
     *
     * fonction abstraite permettant  d'enregistrer dans un fichier
     * la valeur définissant l'ensemble sur la grille
     * La version light permet d'enregistrer seulement les points
     * appartenant au noyau de viabilit�
     * Dans le cal d'ensemble minces par exemple
     * cela permet de diminuer la taille des fichiers
     * Dans le cas de l'ensemble épigraphique la valeur réelle est enregistrée
     * Dans le cas d'ensemble booléen, la fonction caractéristique est enregistrée
     *
     * \param[in] fileName : chaine de caracteres, le nom du fichier
     * pour l'enregistrement
     * Doit être impéematée par chaque classe qui hérite
     */
    virtual void saveValOnGridLight(const string &fileName) const;
    /*!
     * \brief Fonction générique qui permet de calculer les coordonnées entiéres d'un point de la grille é
     * partir de son numéro
     * @param[in] num numéro du point
     * @param[out] res pointeur sur l'espace mémoire dans lequel les coordonnées devront être stockées
     *
     * Remarques:
     *  -  Supposons que les nombres de points de grille par axe sont représentés par le vecteur
     *  \f$(n_0,n_1,\dots, n_{d-1})\f$ . La numérotation est supposée être  dans l'ordre alpha-numérique des
     *  coordonnées entiéres des points
     *   \f$ (i_0,i_1,...i_{d-1})\f$ avec pour tout \f$j=0,\dots, d-1\f$ , \f$ i_j = 0,\dots, n_j-1\f$.
     *   Cette  fonction calcule l'inverse de la transformation suivante
     *   \f[
     *   (i_0,i_1,...i_{d-1}) \mapsto i_0+i_1*n_0+i_2*n_1*n_0+\dots +i_{d_1}\prod_{j=0}^{d-2}n_j
     *   \f]
     *
     *  - attention! é réserver la mémoire correctement  en accord avec la dimension du probléme avant de passer l'adresse
     *  en argument é cette fonction
     */
    void numToIntCoords(unsigned long long int num,
                        unsigned long long int *res) const;

    /*!
     * \brief Fonction générique qui permet de calculer les coordonnées entiéres et réelles
     *  d'un point de la grille é
     * partir de son numéro
     * @param[in] num numéro du point
     * @param[out] resI pointeur sur l'espace mémoire dans lequel les coordonnées entiéres devront être stockées
     * @param[out] resD pointeur sur l'espace mémoire dans lequel les coordonnées  réelles devront être stockées
     *
     * Remarques:
     *  -  Supposons que les nombres de points de grille par axe sont représentés par le vecteur
     *  \f$(n_0,n_1,\dots, n_{d-1})\f$ . La numérotation est supposée être  dans l'ordre alpha-numérique des
     *  coordonnées entiéres des points
     *   \f$ (i_0,i_1,...i_{d-1})\f$ avec pour tout \f$j=0,\dots, d-1\f$ , \f$ i_j = 0,\dots, n_j-1\f$.
     *   Cette  fonction calcule l'inverse de la transformation suivante
     *   \f[
     *   (i_0,i_1,...i_{d-1}) \mapsto i_0+i_1*n_0+i_2*n_1*n_0+\dots +i_{d_1}\prod_{j=0}^{d-2}n_j
     *   \f]
     *
     *   Puis, sachant  \f$ (li_0,\dots, li_{d-1}\f$ les coordonnées du coin inférieur  de la grille et du pas de discrétisation
     *     \f$ (h_0,\dots, h_{d-1})\f$ on reconstruit
     *   les coordonnées réelles :
     *   \f[
     *   x_j=li_j+i_j*h_j,\ \ j=0,\dots, d-1
     *   \f]
     *
     *  - attention! é réserver la mémoire correctement  en accord avec la dimension du probléme avant de passer l'adresse
     *  en argument é cette fonction
     */
    void numToIntAndDoubleCoords(unsigned long long int num,
                                 unsigned long long int *resI, double *resD) const;

    /*!
     * \brief Fonction qui calcule le numéro d'un point é partir de ses coordonnées entiéres
     *
     * @param coords pointeur sur l'adresse mémoire oé sont stockées les  coordonnées entiéres d'un point
     * @param res numéro obtenu
     *
     * Remarques:
     *  -  Supposons que les nombres de points de grille par axe sont représentés par le vecteur
     *  \f$(n_0,n_1,\dots, n_{d-1})\f$ . La numérotation est supposée être  dans l'ordre alpha-numérique des
     *  coordonnées entiéres des points
     *   \f$ (i_0,i_1,...i_{d-1})\f$ avec pour tout \f$j=0,\dots, d-1\f$ , \f$ i_j = 0,\dots, n_j-1\f$.
     *   Cette  fonction calcule  la transformation suivante
     *   \f[
     *   (i_0,i_1,...i_{d-1}) \mapsto i_0+i_1*n_0+i_2*n_1*n_0+\dots +i_{d_1}\prod_{j=0}^{d-2}n_j
     *   \f]
     */
    unsigned long long int intCoordsToNum(const unsigned long long int *coords) const;


    void intCoordsToDoubleCoords(const unsigned long long int *coords, double *coordsDouble) const;
    
    /*!
     * \brief Fonction  qui identifie la maille  de la grille qui contient le point
     * de coordonnées réelles données
     *
     * @param coords  pointeur sur l'adresse mémoire oé sont les coordonnées  réelles du point é localiser
     * @return le numéro du coin inférieur de la maille qui contient le point
     */
    unsigned long long int localizePoint(const double *coords) const;

    /*!
     * Fonction  technqiue qui précalcule certaines données utiles aux reprérages classiques dans
     * la grille :
     *
     *  -indicesDecalCell : décalages de numéro entre un coin inférieur d'une maille et tous les autres vertex
     *  d'une méme maille  (\f$ 2^d\f$ voisins)
     *  -indicesDecalAxes  : décalages de numéros pour les voisins par axe d'un point donné  (\f$ 2d\f$ voisins)
     *  -indicesDecal      : décalages vers tous les voisins dans la grille d'un point donné (\f$ 3^d\f$ voisins)
     * 	-lesDecalagesCell  : décalcages de coordonnées entiéres entre le coins inférieur d'une maille
     * 	 et  tous les autres vertex
     * 	-lesDecalagesAxes  : décalcage de coordonnées entiéres entre un point et tous ses voisins le long des axes
     *
     */

    virtual unsigned long long int getNearestPointInSet(const double *coords) const;

    void computeGridShifts();
    bool ArePointsInTheSameCell(const double *coords1, const double *coords2) const;

    /*!
     * \brief Fonction qui permet de corriger les dépassements éventuels de coordonnées qui
     * sont des variables périodiques dans le modéle
     *
     * @param vect coordonnées réelles
     *
     * La fonction vérifie quelles sont les variables périodiques , en utilisant l'attribut periodic
     */
    void periodizePoint(double *vect) const;

    /*!
     * Fonction qui vérifie si un vecteur de coordonénes données appartient au pavé délimitant la grille
     *
     * @param coords coordonnées réelles du vecteur
     * @return booléen selon l'appartenance ou non de vecteur au pavé de la grille
     */
    bool isPointInGrid(const double *coords) const;
    bool isPointInGrid_fd(const unsigned long long int *coords) const;



    /*!
     * Méthode d'accès
     * @return Le nombre total de mailles de la grille
     */
    unsigned long long int getNbTotalCells() const;

    /*!
     * Méthode d'accès
     * @return Le nombre total de points de la grille
     */
    unsigned long long int getNbTotalPoints() const;

    /*!
     * Méthode d'accès
     * @return Le pas de discrétisation maximal de la grille
     */
    double getMaxStep() const;

    /*!
     * Nombre total de mailles de la grille
     */
    unsigned long long int nbTotalCells;
    /*!
     * Nombre de mailles par axe
     */
    unsigned long long int *nbCells;

    /***********************************************************
     * Accesseurs  communs
     ***************************************************************/

    long long int *indicesDecal;

    string filePrefix;
    const unsigned long long int* getNbPoints() const;

    vector<int> dirs;
    vector<double> values;
    vector<unsigned long long int> values_fd;

    long long int *indicesDecalCell;

    unsigned long long int *indicesDecalAxes;

    unsigned long long int **lesDecalagesCell;
    unsigned long long int **lesDecalagesAxes;

    bool *sortieOKinf;
    bool *sortieOKsup;
    bool isPointInGridWithConstr(const double *coords) const;

    bool unboundedDomain;

};

#endif /* GRID_H_ */
