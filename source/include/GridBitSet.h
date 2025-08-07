/*
 * GridBitSet.h
 *
 * *    VIABLAB : a numerical library for Mathematical Viability Computations
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
 *      Author: ANYA
 */

#ifndef GRIDBITSET_H_
#define GRIDBITSET_H_

#include "Grid.h"
#include "Params.h"

class Grid_BitSet final : public Grid
{
public:
    Grid_BitSet(const gridParams &gp);
    virtual ~Grid_BitSet();
    virtual void printGrid(void) const;

    void addPointToSet(const unsigned long long int *coords, double value);
    void addPointToSet(unsigned long long int pos, double value);
    int getNumberOfViableNeighbors(unsigned long long int pos) const;
    virtual bool isInSet(const unsigned long long int *coords) const;
    virtual unsigned long long int getNearestPointInSet(const double *coords) const;

    void findNearestViabPointInCell(const double *startPoint, const double *currentPoint,
	double *newPoint, double (*dynConstraints)(const double*, double*)) const;
    virtual void savePointsList(const string &fileName) const;
    virtual void saveValOnGrid(const string &fileName) const;
    virtual void saveValOnGridLight(const string &fileName) const;

    boost::dynamic_bitset<> analyseTrameMasque(unsigned long long int posX) const;
    unsigned long long int getDirTram() const;
    unsigned long long int getLongTrame() const;
    boost::dynamic_bitset<>**getGridTab();
    boost::dynamic_bitset<>**getGridTabNew();
    const boost::dynamic_bitset<>* const *getGridTab() const;
    const boost::dynamic_bitset<>* const *getGridTabNew() const;
    void setPoint(const unsigned long long int *coords, bool val);
    void setPoint(unsigned long long int posX, unsigned long long int posTrame,
	    bool val);

    void saveCoupe(const string &nomFichier) const;
    void saveBoundary(const string &nomFichier) const;
    void saveCoupeBoundary(const string &nomFichier) const;
    void saveProjetion(const string &fileName, const unsigned long long int *projection) const;
    void loadSet(const string &fileName);
    void refineVoidGridTab();

    const unsigned long long int* getNbPointsSubGrid() const;
    unsigned long long int getNbPointsTotalSubGrid() const;
    boost::dynamic_bitset<> analyseTrameMasqueBis(unsigned long long int posX,
	    bool print) const;
    void refineTrameMasque(unsigned long long int posX, int parite);

    void copyGrid(const boost::dynamic_bitset<> * const *grIn,
	    boost::dynamic_bitset<> **grOut) const;

    void refine();

private:
    int cp1; //pas de calcul des d�riv�es partielles: x+hx/pas1 , y+hy/p1
    int pp1;
    void computeSubGridShifts();

    int pas1;       //pas de calcul des d�riv�es partielles: x+hhx/pas1

    int nbPointsCubeSub;
    int pow3Sub;
    long long int *indicesDecalSub;
    long long int *indicesDecalCellSub;

    unsigned long long int *indicesDecalAxesSub;

    unsigned long long int **lesDecalagesCellSub;
    unsigned long long int **lesDecalagesAxesSub;
    /*!
     * Nombre total de mailles de la grille
     */
    unsigned long long int nbTotalCellsSub;
    /*!
     * Nombre de mailles par axe
     */
    unsigned long long int *nbCellsSub;

    int dirTramage;
    unsigned long long int *twoPowers;
    unsigned long long int nbP_dm1;
    unsigned long long int longTrame;
    double li_trame;
    double ls_trame;
    boost::dynamic_bitset<> **gridTab;
    boost::dynamic_bitset<> **gridTabNew;
    int nbOMPThreads;

    unsigned long long int *nbPointsSubGrid;
    unsigned long long int *px;

    unsigned long long int nbPointsTotalSubGrid;

    unsigned long long int *tempIntCoords;

    /*!
     * \brief Fonction g�n�rique qui permet de calculer les coordonn�es enti�res d'un point de la grille �
     * partir de son num�ro
     * @param[in] num num�ro du point
     * @param[out] res pointeur sur l'espace m�moire dans lequel les coordonn�es devront �tre stock�es
     *
     * Remarques:
     *  -  Supposons que les nombres de points de grille par axe sont repr�sent�s par le vecteur
     *  \f$(n_0,n_1,\dots, n_{d-1})\f$ . La num�rotation est suppos�e �tre  dans l'ordre alpha-num�rique des
     *  coordonn�es enti�res des points
     *   \f$ (i_0,i_1,...i_{d-1})\f$ avec pour tout \f$j=0,\dots, d-1\f$ , \f$ i_j = 0,\dots, n_j-1\f$.
     *   Cette  fonction calcule l'inverse de la transformation suivante
     *   \f[
     *   (i_0,i_1,...i_{d-1}) \mapsto i_0+i_1*n_0+i_2*n_1*n_0+\dots +i_{d_1}\prod_{j=0}^{d-2}n_j
     *   \f]
     *
     *  - attention! � r�server la m�moire correctement  en accord avec la dimension du probl�me avant de passer l'adresse
     *  en argument � cette fonction
     */
    inline void numToIntCoords_dm1(unsigned long long int num,
	    unsigned long long int *res) const;
    inline void numToSubGridCoords(unsigned long long int num,
	    unsigned long long int *res) const;

    /*!
     * \brief Fonction g�n�rique qui permet de calculer les coordonn�es enti�res et r�elles
     *  d'un point de la grille �
     * partir de son num�ro
     * @param[in] num num�ro du point
     * @param[out] resI pointeur sur l'espace m�moire dans lequel les coordonn�es enti�res devront �tre stock�es
     * @param[out] resD pointeur sur l'espace m�moire dans lequel les coordonn�es  r�elles devront �tre stock�es
     *
     * Remarques:
     *  -  Supposons que les nombres de points de grille par axe sont repr�sent�s par le vecteur
     *  \f$(n_0,n_1,\dots, n_{d-1})\f$ . La num�rotation est suppos�e �tre  dans l'ordre alpha-num�rique des
     *  coordonn�es enti�res des points
     *   \f$ (i_0,i_1,...i_{d-1})\f$ avec pour tout \f$j=0,\dots, d-1\f$ , \f$ i_j = 0,\dots, n_j-1\f$.
     *   Cette  fonction calcule l'inverse de la transformation suivante
     *   \f[
     *   (i_0,i_1,...i_{d-1}) \mapsto i_0+i_1*n_0+i_2*n_1*n_0+\dots +i_{d_1}\prod_{j=0}^{d-2}n_j
     *   \f]
     *
     *   Puis, sachant  \f$ (li_0,\dots, li_{d-1}\f$ les coordonn�es du coin inf�rieur  de la grille et du pas de discr�tisation
     *     \f$ (h_0,\dots, h_{d-1})\f$ on reconstruit
     *   les coordonn�es r�elles :
     *   \f[
     *   x_j=li_j+i_j*h_j,\ \ j=0,\dots, d-1
     *   \f]
     *
     *  - attention! � r�server la m�moire correctement  en accord avec la dimension du probl�me avant de passer l'adresse
     *  en argument � cette fonction
     */
    void numToIntAndDoubleCoords_dm1(unsigned long long int num,
	    unsigned long long int *resI, double *resD) const;

    /*!
     * \brief Fonction qui calcule le num�ro d'un point � partir de ses coordonn�es enti�res
     *
     * @param coords pointeur sur l'adresse m�moire o� sont stock�es les  coordonn�es enti�res d'un point
     * @param res num�ro obtenu
     *
     * Remarques:
     *  -  Supposons que les nombres de points de grille par axe sont repr�sent�s par le vecteur
     *  \f$(n_0,n_1,\dots, n_{d-1})\f$ . La num�rotation est suppos�e �tre  dans l'ordre alpha-num�rique des
     *  coordonn�es enti�res des points
     *   \f$ (i_0,i_1,...i_{d-1})\f$ avec pour tout \f$j=0,\dots, d-1\f$ , \f$ i_j = 0,\dots, n_j-1\f$.
     *   Cette  fonction calcule  la transformation suivante
     *   \f[
     *   (i_0,i_1,...i_{d-1}) \mapsto i_0+i_1*n_0+i_2*n_1*n_0+\dots +i_{d_1}\prod_{j=0}^{d-2}n_j
     *   \f]
     */
    unsigned long long int intCoordsToNum_dm1(const unsigned long long int *coords) const;
    /*!
     * \brief Fonction  qui identifie la maille  de la grille qui contient le point
     * de coordonn�es r�elles donn�es
     *
     * @param coords  pointeur sur l'adresse m�moire o� sont les coordonn�es  r�elles du point � localiser
     * @return le num�ro du coin inf�rieur de la maille qui contient le point
     */

    };

#endif /* GRIDBITSET_H_ */
