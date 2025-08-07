/*
 * ViabiBitSet.cpp
 *
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
 *  Created on: reated on: 27 jul. 2013
 *      Author: Anna DESILLES
 */

#include "../include/ViabiBitSet.h"
#include "../include/WeakDeclarations.h"
#include "../include/ControlPickerBitSet.h"
#include "../include/ViabiBitSetTrajectoryHelper.h"

ViabiBitSet::~ViabiBitSet()
    {
    // TODO Auto-generated destructor stub
    }

void ViabiBitSet::initialiseTarget()
    {

    /*!
     *  cette fonction initialise la base de donn�es pour la dynamique. Elle
     *  ajoute dnas la base les premiers points pour lesquels la fonction valeur
     *  est r�elle. Au d�but de l'algorithme de bassin de capture
     *    seuls les points de la  cible  ont une fonction valeur r�elle
     *
     */
    initialiseTargetPointList();
    }

void ViabiBitSet::setK0_fd()
    {
    int dirTramage = grid->getDirTram();

    unsigned long long int tailleTrame = grid->getNbPointsTotalSubGrid();

    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();
    unsigned long long int longTrame = grid->getLongTrame();
    boost::dynamic_bitset<> masque;
    unsigned long long int indice[dim - 1];
    unsigned long long int coordDiscretes[dim];
    bool testK;

    for (unsigned long long int posX = 0; posX < tailleTrame; posX++)
	{

	numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

	for (int j = 0; j < dirTramage; j++)
	    {
	    coordDiscretes[j] = indice[j];

	    }

	for (int j = dirTramage + 1; j < dim; j++)
	    {
	    coordDiscretes[j] = indice[j - 1];
	    }

	for (unsigned long long int k = 0; k < longTrame; k++)
	    {
	    coordDiscretes[dirTramage] = k;
	    testK = (dynsys->constraintsX_fd(coordDiscretes) < PLUS_INF);
	    grid->setPoint(posX, k, testK);
	    }

	}	// fin de for de parcours de la trame
    }

void ViabiBitSet::setK0()
    {

    cout << " Initialisation de contraintes" << endl;
    int dirTramage = grid->getDirTram();

    unsigned long long int tailleTrame = grid->getNbPointsTotalSubGrid();

    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    //calcul de la taille e prevoir pour les coordonnees des indices de debut de parcours
    //que la methode GPU va renvoyer

    double *limInf = grid->limInf;
    double *gridStep = grid->step;

    double xCoordsDouble[dim];

    //    cout<<"masque points enleves cree"<<masquePointsEnleves;
    unsigned long long int longTrame = grid->getLongTrame();
    boost::dynamic_bitset<> masque;
    unsigned long long int indice[dim - 1];
    bool testK;

    unsigned long long int posX = 0;

    for (posX = 0; posX < tailleTrame; posX++)
	{

	numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

	for (int j = 0; j < dirTramage; j++)
	    {
	    xCoordsDouble[j] = limInf[j] + indice[j] * gridStep[j];
	    }

	for (int j = dirTramage + 1; j < dim; j++)
	    {
	    xCoordsDouble[j] = limInf[j] + indice[j - 1] * gridStep[j];
	    }

	for (unsigned long long int k = 0; k < longTrame; k++)
	    {
	    xCoordsDouble[dirTramage] = limInf[dirTramage] + k * gridStep[dirTramage];

	    testK = (dynsys->constraintsX(xCoordsDouble) < PLUS_INF);

	    grid->setPoint(posX, k, testK);
	    }
	}	// fin de for de parcours de la trame
    }

void ViabiBitSet::initialiseConstraints()
    {
    if ((this->dynsys->getDynType() == CC) || (this->dynsys->getDynType() == DC))
	{
	setK0();
	}
    else
	{
	setK0_fd();
	}
    }

void ViabiBitSet::testK0()
    {
    int dirTramage = grid->getDirTram();

    unsigned long long int tailleTrame = grid->getNbPointsTotalSubGrid();

    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    double *limInf = grid->limInf;
    double *gridStep = grid->step;

    unsigned long long int longTrame = grid->getLongTrame();

    unsigned long long int posX = 0;

    boost::dynamic_bitset<> **gridTab = grid->getGridTab();
#pragma omp parallel for num_threads(nbOMPThreads)  private(posX) shared( gridTab, limInf, gridStep, nbPointsSub,dirTramage , longTrame, tailleTrame) default(none)

    for (posX = 0; posX < tailleTrame; posX++)
	{
	unsigned long long int indice[dim - 1];
	double xCoordsDouble[dim];
	bool testK;

	numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

	for (int j = 0; j < dirTramage; j++)
	    {
	    xCoordsDouble[j] = limInf[j] + indice[j] * gridStep[j];
	    }

	for (int j = dirTramage + 1; j < dim; j++)
	    {
	    xCoordsDouble[j] = limInf[j] + indice[j - 1] * gridStep[j];
	    }
	for (unsigned long long int k = 0; k < longTrame; k++)
	    {
	    xCoordsDouble[dirTramage] = limInf[dirTramage] + k * gridStep[dirTramage];
	    testK = (dynsys->constraintsX(xCoordsDouble) < PLUS_INF);
	    if (!testK)
		grid->setPoint(posX, k, testK);
	    }
	}	// fin de for de parcours de la trame
    }

ViabiBitSet::ViabiBitSet(ParametersManager *pm)
    {

    /*
     * instanciation du système dynamique
     */
    modelParams = pm;

    const systemParams *sp = pm->getSystemParameters();
    const algoViabiParams *avp = pm->getAlgoParameters();
    const controlParams *cp = pm->getControlParameters();
    const gridParams *gp = pm->getGridParameters();

    dim = gp->DIM;
    grid = new Grid_BitSet(*gp);
    dynsys = new SysDyn(*sp, dim, *cp, grid);
    
    InitViabiBitSet(*avp);
    }

void ViabiBitSet::InitViabiBitSet(const algoViabiParams &avbp)
    {

    /*!
     * On initialise la structure  qui servira au stockage temporaire de l'image discrete d'un point de grille
     * Les tableaux sont initialis�s  avec leur taille maximale possible :
     * �gale au nb  de controles.
     * C'est la valeur effective  de nbImageCells , calcul� � chaque evaluation de l'image discrete
     * d'un point qui  servira  � lire et remplir  correctement
     * la bonne partie de ces tableaux
     */

    dim = grid->dim;
    dimC = dynsys->getDimC();
    filePrefix = avbp.FILE_PREFIX;
    nbOMPThreads = avbp.NB_OMP_THREADS;

    /*!
     * On initialise la structure  qui servira au stockage temporaire de l'image discrete d'un point de grille
     * Les tableaux sont initialis�s  avec leur taille maximale possible :
     * �gale au nb  de controles.
     * C'est la valeur effective  de nbImageCells , calcul� � chaque evaluation de l'image discrete
     * d'un point qui  servira  � lire et remplir  correctement
     * la bonne partie de ces tableaux
     */

    pointDI.tabImageCells = new unsigned long long int[dynsys->getTotalNbPointsC()];

    doubleVect = new double[dim];
    doubleVect1 = new double[dim];
    imageCells = new unsigned long long int[dynsys->getTotalNbPointsC()];    
    //cout<< " viab bit set contructeur ok \n";
    }

void ViabiBitSet::printViabiInfo()
    {
    grid->printGrid();
    }

void ViabiBitSet::noyauViabi_FD(int nbArret)
    {
    int dirTramage = grid->getDirTram();
    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    unsigned long long int tailleTrame = grid->getNbPointsTotalSubGrid();

    bool testF;
    unsigned long long int longTrame = grid->getLongTrame();

    unsigned long long int nbC = dynsys->getTotalNbPointsC();
    boost::dynamic_bitset<> masque;
    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(longTrame, 0);
    boost::dynamic_bitset<> **gridTab = grid->getGridTab();
    unsigned long long int indice[dim];
    unsigned long long int **controlIntCoords = dynsys->getControlIntCoords();

    unsigned long long int coordDiscretes[dim], imageXU[dim];

    unsigned long long int compteComm;
    bool testNonVide;
    int nbIter = 0;
    unsigned long long int comptEtats = 0, comptEnleves = nbArret + 1;
    while (comptEnleves > (unsigned long long int) nbArret)
	{

	comptEnleves = 0;
	comptEtats = 0;
	int dirParcours = -1;
	//for( unsigned long long  int posX=0+((1-dirParcours)/2)*(tailleTrame-1);dirParcours*(posX-((1+dirParcours)/2)*(tailleTrame-1))<=0;posX+=dirParcours)
	for (unsigned long long int posX = 0; posX < tailleTrame; posX++)
	    {
	    if (!(*gridTab[posX]).none())
		{
		comptEtats++;
		testF = false;
		numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

		masque = grid->analyseTrameMasque(posX);

		masquePointsEnleves->set();

		if (masque.none() | testF)
		    {
		    //cout<<" rien e analyser posx= "<<posX<<" \n";
		    }
		else
		    {
		    for (int j = 0; j < dirTramage; j++)
			{
			coordDiscretes[j] = indice[j];

			}

		    for (int j = dirTramage + 1; j < dim; j++)
			{
			coordDiscretes[j] = indice[j - 1];

			}

		    for (unsigned long long int k = 0; k < longTrame; k++)
			{

			if (masque[k])
			    {

			    coordDiscretes[dirTramage] = k;
			    testNonVide = false;
			    compteComm = 0;
			    while (!testNonVide && (compteComm < nbC))
				{
				if (dynsys->constraintsXU_fd(coordDiscretes, controlIntCoords[compteComm]) < PLUS_INF)
				    {
				    //  cout<< "  controle "<<controlIntCoords[compteComm][0]<<" "<<controlIntCoords[compteComm][1]<< " est autorise \n";
				    dynsys->dynamics_fd(coordDiscretes, controlIntCoords[compteComm], imageXU);
				    // cout<< "  image "<<imageXU[0]<<" "<<imageXU[1]<< " est autorise \n";

				    if (grid->isPointInGrid_fd(imageXU))
					{
					testNonVide = grid->isInSet(imageXU);
					}
				    else
					{
					testNonVide = false;
					}
				    }
				compteComm++;
				}
			    if (!testNonVide)
				{
				masquePointsEnleves->set(k, false);
				comptEnleves++;
				}
			    }	// fin de if masque[k]
			}	// fin de for  de parcours de masque

		    if (masquePointsEnleves->count() < (unsigned long long int) longTrame)
			{
			(*gridTab[posX]) &= (*masquePointsEnleves);
			}
		    }

		}	//fin de if la trame n'est pas vide
	    }	// fin de for de parcours de la trame
	nbIter++;
	cout << " nbPoint enlevee " << comptEnleves << endl;
	}
    cout << "fini nbIter=" << nbIter;
    }

bool ViabiBitSet::findViabImagePoint_noControl(double *xCoordsDouble, bool print)
    {

    //cout<< " findImage no control debut"<<endl;
    print = false;
    int nbPointsCube = (int) pow(2.0, dim);		//pow(2.0, dim);
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */
    unsigned long long int testI[dim];
    double testV[dim];
    /*
     * tableau de coordonnées  définissant
     * la grille de contrôles
     */
    double **controlCoords = dynsys->getControlCoords();
    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    //cout<< " nb controls total "<<nbCTotal<<endl;
    /*
     * indices de déplacement pour parcourir les sommets d'une malle à partir
     * du sommet inf
     */
    const long long int *indicesDecalCell = grid->getIndicesDecalCell();
    /*
     * coordonnées du points corent de la trajectoire
     */

    /*
     * numéros de mailles
     */
    int cellNum, posTemp;

    double rho = dynsys->calculRho_local(xCoordsDouble);

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;
    double doubleVect1[dim];
    double hMax = grid->maxStep;

    unsigned long long int cu = 0;
    /*
     * On recherche la plus grande puissance de 2 plus grande  que le nombre totale de
     * contrôles à parcourir
     */

    int stepCu = 1;

    //cout<< " parcours controles debut = "<<cu<< " step= "<<stepCu<< " nb controls total "<<nbCTotal<<endl;
    while ((cu < nbCTotal) && !testNonVide)
	{
	/*
	 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	 * au point en cours
	 */
	if (dynsys->constraintsXU(xCoordsDouble, controlCoords[cu]) < PLUS_INF)
	    {
	    //cout<< " constr XU OK rho = "<<rho<<endl;
	    //cout<< " xCoords "; printVector(xCoordsDouble, dim);
	    /*
	     * calcul du successeur  du point en cours
	     */
	    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[cu], doubleVect1, rho);
	    //cout<< "image "; printVector(doubleVect1, dim);

	    if (grid->isPointInGrid(doubleVect1))
		{
		/*
		 *  le sucesseur est dans la grille de calcul
		 */
		if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
		    {

		    /*
		     * le successeur vérifie les contraintes
		     * On identifie la maille où il se trouve
		     */
		    cellNum = grid->localizePoint(doubleVect1);

		    /*
		     * On parcours les sommets de la maille
		     * autour du sucesseur pour voir s'il y a des
		     * points viables
		     */
		    int ii = 0;
		    while (ii < nbPointsCube && !testNonVide)
			{
			posTemp = cellNum + indicesDecalCell[ii];
			grid->numToIntAndDoubleCoords(posTemp, testI, testV);
			if (dynsys->dynConstraintsForTraj(xCoordsDouble, testV) < PLUS_INF)
			    {
			    double dist = 0.0;
			    for (int k = 0; k < dim; k++)
				{
				dist = max(dist, abs(testV[k] - doubleVect1[k]));
				}
			    testNonVide = grid->isInSet(testI) && (dist <= hMax / 2.0);

			    }
			ii++;
			}
		    }
		}
	    else
		{
		//cout<< " sortie du domaine : autorisation : "<<grid->unboundedDomain<< "point "<< doubleVect1[0]<<" "<<doubleVect1[1]<<endl;
		testNonVide = grid->unboundedDomain && grid->isPointInGridWithConstr(doubleVect1) && (dynsys->constraintsX(doubleVect1) < PLUS_INF);
		}
	    }
	cu += stepCu;
	}			//fin de parcours de tous les contrôles

    return testNonVide;

    }

bool ViabiBitSet::findViabImagePoint(double *xCoordsDouble, bool print)
{
    print = false;
    int nbPointsCube = (int) pow(2.0, dim);			//pow(2.0, dim);
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */
    unsigned long long int *testI = new unsigned long long int[dim];
    double *testV = new double[dim];
    /*
     * tableau de coordonnées  définissant
     * la grille de contrôles
     */
    double **controlCoords = dynsys->getControlCoords();
    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();

    /*
     * indices de déplacement pour parcourir les sommets d'une malle à partir
     * du sommet inf
     */
    const long long int *indicesDecalCell = grid->getIndicesDecalCell();
    /*
     * coordonnées du points corent de la trajectoire
     */

    /*
     * numéros de mailles
     */
    int cellNum;

    int posTemp;
    double rho = dynsys->calculRho_local(xCoordsDouble);

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;

    double *doubleVect1 = new double[dim];

    double hMax = grid->maxStep;

    unsigned long long int cu = 0;

    /*
     * On recherche la plus grande puissance de 2 plus grande  que le nombre totale de
     * contrôles à parcourir
     */

    int p = 0;
    unsigned long long int pow2 = 1;
    while (pow2 <= nbCTotal)
	{
        p++;
        pow2 *= 2;
	}
    p--;
    pow2 = pow2 / 2;

    int first = 0;
    int stepCu = pow2;
    cu = first;

    while ((cu < nbCTotal) && !testNonVide)
	{
        /*
         * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
         * au point en cours
         */
        if (dynsys->constraintsXU(xCoordsDouble, controlCoords[cu]) < PLUS_INF)
	    {
            /*
             * calcul du successeur  du point en cours
             */
            (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[cu], doubleVect1, rho);
            if (grid->isPointInGrid(doubleVect1))
            {
                /*
                 *  le sucesseur est dans la grille de calcul
                 */
                if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
                {
                    /*
                     * le successeur vérifie les contraintes
                     * On identifie la maille où il se trouve
                     */
                    cellNum = grid->localizePoint(doubleVect1);
                    // cout<< " num cellule "<<cellNum<<endl;

                    /*
                     * On parcours les sommets de la maille
                     * autour du sucesseur pour voir s'il y a des
                     * points viables
                     */
                    int ii = 0;
                    while (ii < nbPointsCube && !testNonVide)
                    {
                        posTemp = cellNum + indicesDecalCell[ii];
                        grid->numToIntAndDoubleCoords(posTemp, testI, testV);
                        if (dynsys->dynConstraintsForTraj(xCoordsDouble, testV) < PLUS_INF)
                        {
                            double dist = 0.0;
                            for (int k = 0; k < dim; k++)
                            {
                                dist = max(dist, abs(testV[k] - doubleVect1[k]));
                            }
                            testNonVide = grid->isInSet(testI) && (dist <= hMax / 2.0);
                            ii++;
                        }
                    }
                }
            }
            else
            {
                testNonVide = grid->unboundedDomain && grid->isPointInGridWithConstr(doubleVect1) && (dynsys->constraintsX(doubleVect1) < PLUS_INF);
            }

	    }
        cu += stepCu;
	}			//fin de parcours de tous les contrôles

    if (!testNonVide)
	{
	first = stepCu / 2;
	while (stepCu > 1 && !testNonVide)
	    {
	    //    cout<< " parcours controles debut = "<<first<< " step= "<<stepCu<<endl;
	    cu = first;
	    while (cu < nbCTotal && !testNonVide)
		{
		/*
		 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
		 * au point en cours
		 */
		if (dynsys->constraintsXU(xCoordsDouble, controlCoords[cu]) < PLUS_INF)
		    {
		    /*
		     * calcul du successeur  du point en cours
		     */
		    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[cu], doubleVect1, rho);
		    if (grid->isPointInGrid(doubleVect1))
			{
			/*
			 *  le sucesseur est dans la grille de calcul
			 */
			if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
			    {
			    /* cout<< " image du point ";
			     for(int k=0;k<dim;k++)
			     {
			     cout<< " "<<doubleVect1[k];
			     }*/
			    /*
			     * le successeur vérifie les contraintes
			     * On identifie la maille où il se trouve
			     */
			    cellNum = grid->localizePoint(doubleVect1);
			    // cout<< " num cellule "<<cellNum<<endl;

			    /*
			     * On parcours les sommets de la maille
			     * autour du sucesseur pour voir s'il y a des
			     * points viables
			     */
			    int ii = 0;
			    while (ii < nbPointsCube && !testNonVide)
				{
				posTemp = cellNum + indicesDecalCell[ii];
				grid->numToIntAndDoubleCoords(posTemp, testI, testV);
				if (dynsys->dynConstraintsForTraj(xCoordsDouble, testV) < PLUS_INF)
				    {
				    double dist = 0.0;
				    for (int k = 0; k < dim; k++)
					{
					dist = max(dist, abs(testV[k] - doubleVect1[k]));
					}
				    //testNonVide= grid->isInSet(testI) && (dist<=hMax/2.0);

				    testNonVide = grid->isInSet(testI);

				    ii++;
				    }
				}

			    }
			}
		    }
		cu += stepCu;
		}			//fin de parcours de tous les contrôles

	    first = first / 2;
	    stepCu = stepCu / 2;
	    }
	}

    delete[] doubleVect1;
    delete[] testV;
    delete[] testI;

    return testNonVide;

    }

bool ViabiBitSet::findGuarantedViabImagePoint(double *xCoordsDouble, bool print)
    {
    print = false;
    int nbPointsCube = (int) pow(2.0, dim);			//pow(2.0, dim);
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */
    unsigned long long int *testI = new unsigned long long int[dim];
    double *testV = new double[dim];
    /*
     * tableau de coordonnées  définissant
     * la grille de contrôles
     */
    double **controlCoords = dynsys->getControlCoords();
    double **tychCoords = dynsys->getTychCoords();
    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int nbTyTotal = dynsys->getTotalNbPointsTy();

    /*
     * indices de déplacement pour parcourir les sommets d'une malle à partir
     * du sommet inf
     */
    const long long int *indicesDecalCell = grid->getIndicesDecalCell();
    /*
     * coordonnées du points corent de la trajectoire
     */

    /*
     * numéros de mailles
     */
    int cellNum;

    int posTemp;
    double rho = dynsys->calculRho_local(xCoordsDouble);

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;

    double *doubleVect1 = new double[dim];

    double hMax = grid->maxStep;

    unsigned long long int cu = 0;

    /*
     * On recherche la plus grande puissance de 2 plus grande  que le nombre totale de
     * contrôles à parcourir
     */

    int p = 0;
    unsigned long long int pow2 = 1;
    while (pow2 <= nbCTotal)
	{
	p++;
	pow2 *= 2;
	}
    p--;
    pow2 = pow2 / 2;

    int first = 0;
    int stepCu = pow2;
    cu = first;
    bool testAllTych;
    while ((cu < nbCTotal) && !testNonVide)
	{
        /*
         * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
         * au point en cours
         */
        if (dynsys->constraintsXU(xCoordsDouble, controlCoords[cu]) < PLUS_INF)
	    {
            testAllTych = true;
            long long unsigned int ty = 0;
            while (ty < nbTyTotal && testAllTych)
            {
                if(dynsys->constraintsXV_tych(xCoordsDouble, tychCoords[ty]) < PLUS_INF)
                {
                    /*
                     * calcul du successeur  du point en cours
                     */
                    (dynsys->*(dynsys->discretDynamics_tych))(xCoordsDouble, controlCoords[cu], tychCoords[ty], doubleVect1, rho);
                    if (grid->isPointInGrid(doubleVect1))
                    {
                        /*
                         *  le sucesseur est dans la grille de calcul
                         */
                        if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
                        {
                            /*
                             * le successeur vérifie les contraintes
                             * On identifie la maille où il se trouve
                             */
                            cellNum = grid->localizePoint(doubleVect1);
                            // cout<< " num cellule "<<cellNum<<endl;

                            /*
                             * On parcours les sommets de la maille
                             * autour du sucesseur pour voir s'il y a des
                             * points viables
                             */
                            int ii = 0;
                            bool imageInSet = false;
                            while (ii < nbPointsCube && !imageInSet)
                            {
                                posTemp = cellNum + indicesDecalCell[ii];
                                grid->numToIntAndDoubleCoords(posTemp, testI, testV);
                                double dist = 0.0;
                                for (int k = 0; k < dim; k++)
                                {
                                    dist = max(dist, abs(testV[k] - doubleVect1[k]));
                                }
                                imageInSet = grid->isInSet(testI) && (dist <= hMax / 2.0);
                                ii++;
                            }
                            testAllTych &= imageInSet;
                        }
                    }
                    else
                    {
                        testAllTych = grid->unboundedDomain && grid->isPointInGridWithConstr(doubleVect1)
                            && (dynsys->constraintsX(doubleVect1) < PLUS_INF);
                    }
                }

                ty++;
            }
            testNonVide = testAllTych;
	    }

        cu += stepCu;
	}			//fin de parcours de tous les contrôles

    if (!testNonVide)
	{
        first = stepCu / 2;
        while (stepCu > 1 && !testNonVide)
	    {
            //    cout<< " parcours controles debut = "<<first<< " step= "<<stepCu<<endl;
            cu = first;
            while (cu < nbCTotal && !testNonVide)
            {
                if (dynsys->constraintsXU(xCoordsDouble, controlCoords[cu]) < PLUS_INF)
                {
                    /*
                     * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
                     * au point en cours
                     */
                    testAllTych = true;
                    unsigned long long int ty = 0;
                    while (ty < nbTyTotal && testAllTych)
                    {
                        if(dynsys->constraintsXV_tych(xCoordsDouble, tychCoords[ty]) < PLUS_INF)
                        {
                            /*
                             * calcul du successeur  du point en cours
                             */
                            (dynsys->*(dynsys->discretDynamics_tych))(xCoordsDouble, controlCoords[cu], tychCoords[ty], doubleVect1, rho);
                            if (grid->isPointInGrid(doubleVect1))
                            {
                                /*
                                 *  le sucesseur est dans la grille de calcul
                                 */
                                if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
                                {
                                    /* cout<< " image du point ";
                                       for(int k=0;k<dim;k++)
                                       {
                                       cout<< " "<<doubleVect1[k];
                                       }*/
                                    /*
                                     * le successeur vérifie les contraintes
                                     * On identifie la maille où il se trouve
                                     */
                                    cellNum = grid->localizePoint(doubleVect1);
                                    // cout<< " num cellule "<<cellNum<<endl;

                                    /*
                                     * On parcours les sommets de la maille
                                     * autour du sucesseur pour voir s'il y a des
                                     * points viables
                                     */
                                    int ii = 0;
                                    bool imageInSet = false;
                                    while (ii < nbPointsCube && !imageInSet)
                                    {
                                        posTemp = cellNum + indicesDecalCell[ii];
                                        grid->numToIntAndDoubleCoords(posTemp, testI, testV);
                                        double dist = 0.0;
                                        for (int k = 0; k < dim; k++)
                                        {
                                            dist = max(dist, abs(testV[k] - doubleVect1[k]));
                                        }
                                        imageInSet = grid->isInSet(testI) && (dist <= hMax / 2.0);
                                        ii++;
                                    }
                                    testAllTych &= imageInSet;
                                }
                            }
                            else
                            {
                                testAllTych = grid->unboundedDomain && grid->isPointInGridWithConstr(doubleVect1)
                                    && (dynsys->constraintsX(doubleVect1) < PLUS_INF);
                            }

                        }
                        ty++;
                    }
                }
                testNonVide = testAllTych;
                cu += stepCu;
            }			//fin de parcours de tous les contrôles

            first = first / 2;
            stepCu = stepCu / 2;
	    }
	}

    delete[] doubleVect1;
    delete[] testV;
    delete[] testI;

    return testNonVide;

    }

void ViabiBitSet::initialiseTargetPointList()
    {
    /*
     *  cette fonction initialise la base de donn�es pour . Elle
     *  ajoute dnas la base les premiers points pour lesquels la fonction valeur
     *
     *   est r�elle. Au d�but de l'algorithme de bassin de capture
     *    seuls le spoints de la  cible  ont une fonction valeur r�elle
     *
     */

    //ouverture de la base
    unsigned long long int *x = new unsigned long long int[dim];

    double *xReel = new double[dim];

    int totalPointsC = 0;  // nombre de points  de l'espace des commandes

    int totalPointsX = grid->getNbTotalPoints();

    imagePoint currentPoint;
        currentPointsImage.clear();
        currentCellsImage.clear();

    unsigned long long int pos;


    /*
     *  on parcourt  tous les points de l'espace discret  fini
     *   et  on choisit les points o� la fonction cible  renvoie une faleur finie
     */
    for (pos = 0; pos < (unsigned long long int) totalPointsX; pos++)
	{

	// cout<< " pos= "<<pos<<endl;
	/*!
	 * le compteur pos  ets l'unique num�ro entier du point en cours
	 * dans la num�rotation alphab�tique : on parcourt axe par axe
	 */

	/*!
	 * on restitue les  coordonn�es neti�res  du point � partir de son num�ro
	 * ainsi que ses coordonn�es r�elles
	 */
	grid->numToIntAndDoubleCoords(pos, x, xReel);

	if ((*(dynsys->constraintsX))(xReel) < PLUS_INF)
	    {
	    double c = (*(dynsys->target))(xReel);
	    if ((*(dynsys->target))(xReel) < PLUS_INF)
		{
		//				cout<<" x dans la cible posX="<<pos;
		//				for(int i=0;i<dim;i++)
		//				{
		//					cout<<" "<<xReel[i]<<" ";
		//				}
		//				cout<<"\n";

		totalPointsC++;

		/*!
		 * Si la fonction  target() renvoie un r�el <INF on ajoute le point dans la base de donn�es
		 * en appelant la fonction  addPointToSet()
		 */
		//cout<< " coucou\n";
		grid->addPointToSet(x, 1.0);
		//cout<< " coucou1\n";
		/*!
		 * On ajoute �galement ce point � la liste des points qui servira pour le calcul de la premi�re image
		 * \f$ \Phi(C)\f$, voir les m�thodes minTimeEpiLocalRho() et minTimeEpiGlobalRho().
		 */

		//        cout<<" fini ajout point retour creation de base \n";//////system("pause");
		currentPointsImage[pos] = c;

		}
	    }
	else
	    {
	    //cout<<"le point est dans l'obstacle\n";
	    }

	}
    //cout<<"remplissage fini\n";
    //cout<<"total points Cible est "<<totalPointsC<<"\n";

    }


void ViabiBitSet::computeConvexifiedImage(int iter)
    {
    spdlog::info("Computing convexified image of the target set. number of points {}", currentPointsImage.size());
    int posX;

    list<imageCell>::iterator itStart, itNew;
    currentCellsImage.clear();

    double rho;
    map<unsigned long long int, double>::iterator itPoint = currentPointsImage.begin(), itLastPoint = currentPointsImage.end(), itTemp;
    map<unsigned long long int, list<unsigned long long int> >::iterator itCell, itLastCell, itCellsTemp;
    unsigned long long int intPointCoords[dim];
    double doublePointCoords[dim];

    DiscretPointImage *tempPointImage = new DiscretPointImage((unsigned long long int) 0, dim, dynsys, grid);
    map<unsigned long long int, list<unsigned long long int> > *cellsWithControls = tempPointImage->GetImageCellsWithControls();

    while (itPoint != itLastPoint)
	{
	posX = itPoint->first;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	rho = dynsys->calculRho_local(doublePointCoords);
	//	spdlog::info("[Min Time problem] : Timestep =  {}", rho);

	/*!
	 * On calcule l'image discrète  du point
	 */
	double tempVal = itPoint->second + rho;

	tempPointImage->setPointNumAndRebuild(posX);

	itCell = cellsWithControls->begin();
	itLastCell = cellsWithControls->end();
	while (itCell != itLastCell)
	    {
	    unsigned long long int cellNum = itCell->first;
	    // cout<< "convexified image  cellnum = "<< cellNum<< " nb total de pojnts " << grid->getNbTotalPoints()<<" ";
	    if(cellNum < grid->getNbTotalPoints())
		{
		if (auto result = currentCellsImage.find(cellNum); result != currentCellsImage.end())
		    {
		    result->second = min((double) result->second, tempVal);
		    }
		else
		    {
		    currentCellsImage[cellNum] = tempVal;
		    }
		addConvexCombinations(posX, itPoint->second, tempVal, itCell->first, rho);
		}
	    itCell++;
	    }
	itPoint++;
	}
    spdlog::info("[Min Time problem] : Computing convexified image finished. number of points {}", currentPointsImage.size());
    }

void ViabiBitSet::addConvexCombinations(unsigned long long int posX, double pointVal, double newCellVal, unsigned long long int numCell,
	double rho)
    {

    /*!
     * On appelle d'abord deux fois la fonction numToIntAndDoubleCoords()
     * pour calculer les coordonnées réelles du point de départ x  et  du coin inférieur de la maille
     *  appartenant é l'image \f$ \Phi(x)\f$ y.
     */

    double LVal = (newCellVal - pointVal) / rho;

    double doublePointCoords[dim], doubleCellCoords[dim], doubleVect[dim], doubleVect1[dim];
    unsigned long long int tempCoords[dim];

    grid->numToIntAndDoubleCoords(posX, tempCoords, doublePointCoords);
    grid->numToIntAndDoubleCoords(numCell, tempCoords, doubleCellCoords);
    double *gridStep = grid->step;

    list<imageCell>::iterator itNew;
    double dist = 0.;
    /*!
     * Ensuite on calcule le vecteur différence \f$ z=y-x\f$ et sa norme \f$ \|z\|_2 \f$ .
     */
    for (int i = 0; i < dim; i++)
	{
	doubleCellCoords[i] = doubleCellCoords[i] + 0.5 * gridStep[i];
	doubleVect[i] = doubleCellCoords[i] - doublePointCoords[i];
	dist += doubleVect[i] * doubleVect[i];
	}
    /*!
     * Le segment reliant \f$ x\f$ et \f$ y\simeq x-\rho F(x,u)\f$ pour un certain \f$ u\in U(x)\f$  a pour équation
     * \f[
     * x+tz,\ t\in[0,1]
     * \f]
     * On détermine alors un pas de progression \f[
     * \Delta t=\min(0.1, \frac{h_{max}}{\|z\|_2})
     * \f]
     * de faéon é pouvoir passer d'une maille é l'autre en avanéant avec  ce pas le long du segment
     * On parcourt ensuite le segment avec le pas calculé et on ajoute é l'image les mailles croisées
     * avec pour valeur une partie du pas de temps \f$\rho\f$ .
     */
    double deltat = min(0.1, grid->getMaxStep() / (2.0 * sqrt(dist)));
    double t = deltat;
    unsigned long long int newCellNum, lastVisitCellNum = grid->getNbTotalCells() + 1;
    while (t < 1.0)
	{

	for (int i = 0; i < dim; i++)
	    {
	    doubleVect1[i] = doublePointCoords[i] + t * doubleVect[i];
	    }

	if (grid->isPointInGrid(doubleVect1))
	    {
	    /*!
	     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
	     */

	    if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
		{

		/*!
		 * Si l'image  est dans l'ensemble de contraintes sur l'état \f$ K \f$
		 * on calcule le numéro de maille qui contient cette image
		 */
		newCellNum = grid->localizePoint(doubleVect1);
		double tempVal = pointVal + t * rho * LVal;
		if (auto result = currentCellsImage.find(newCellNum); result != currentCellsImage.end())
		    {
		    result->second = min((double) result->second, tempVal);
		    }
		else
		    {
		    currentCellsImage[newCellNum] = tempVal;
		    }
		}
	    }
	t = t + deltat;
	}
    }


void ViabiBitSet::createPointsList()
    {

    spdlog::debug("Create Current Point list. Number of cells  {}", currentCellsImage.size());
    double t1, t2, elapsed_time;

    timeval tim;
    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    unsigned long long int posX, cellNum;
    double testV[dim];
    unsigned long long int testI[dim];
    currentPointsImage.clear();

    int nbPointsCube = (int) pow(2.0, dim);	//pow(2.0, dim);
    map<unsigned long long int, double>::iterator itCell = currentCellsImage.begin(), itLast = currentCellsImage.end();

    const long long int *indicesDecalCell = grid->getIndicesDecalCell();

    while (itCell != itLast)
	{
	cellNum = itCell->first;
	// cout<< "create points liste  cellnum = "<< cellNum<< " nb total de pojnts " << grid->getNbTotalPoints()<<" ";
	double tempVal = itCell->second;
	for (int i = 0; i < nbPointsCube; i++)
	    {
	    posX = cellNum + indicesDecalCell[i];
	    grid->numToIntAndDoubleCoords(posX, testI, testV);
	    if (dynsys->constraintsX(testV) < PLUS_INF)
		{
		if (auto result = currentPointsImage.find(posX); result != currentPointsImage.end())
		    {
		    result->second = min((double) result->second, tempVal);
		    }
		else
		    {
		    currentPointsImage[posX] = tempVal;
		    }
		}
	    }
	itCell++;
	}
    //cout<<endl;
    currentCellsImage.clear();

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));
    spdlog::debug("Create Current Point listfinished done {} sec; Number of points is {} ",
	    elapsed_time, currentPointsImage.size());

    }


void ViabiBitSet::computeTrajectories()
{
    dynsys->setDynamicsForward();
    switch (modelParams->getAlgoParameters()->SET_TYPE) {
    case VIAB:
        computeViableTrajectories();
        break;
    case CAPT:
        spdlog::error("Trajectories for capture bassin problems not implemented yet");
        break;
    case VIABG:
        computeTychasticViableTrajectories();
        break;
    }
}

void ViabiBitSet::computeViableTrajectories()
{
    ostringstream os;
    string fileName;

    const int nbTrajs = modelParams->getNbTrajectories();
    
    for (int tr = 0; tr < nbTrajs; tr++)
    {
        TrajectoryParametersManager tpm{modelParams, tr};
        const trajectoryParams *tp = tpm.getTrajectoryParameters();
        TypeTraj typeTraj = tp->TRAJECTORY_TYPE;
        double T = tp->maxTime;
        ViabiBitSetTrajectoryHelper trajHelper(grid, dynsys, &tpm);
            
        if (isViableDefault(typeTraj))
        {
            os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1 << ".dat";
            fileName = os.str();
            os.str("");

            trajHelper.computeViableTrajectory(tp->INIT_POINT, T, fileName);
        }
        else if (typeTraj == VL)
        {
            os << "../OUTPUT/" << filePrefix << "-traj-H-" << tr + 1 << ".dat";
            fileName = os.str();
            os.str("");

            trajHelper.computeViableTrajectoryHeavy(tp->INIT_POINT, tp->INIT_CONTROL, T, fileName);
        }
        else if (typeTraj == CAUTIOUS || typeTraj == WEIGHTED_CONTROLS_CAUTIOUS) {
            os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1 << ".dat";
            fileName = os.str();
            os.str("");

            trajHelper.computeCautiousTrajectory(tp->INIT_POINT, T, fileName);
        }
        else if (typeTraj == CAUTIOUS_HEAVY || typeTraj == WEIGHTED_CONTROLS_CAUTIOUS_HEAVY) {
            os << "../OUTPUT/" << filePrefix << "-traj-H-" << tr + 1 << ".dat";
            fileName = os.str();
            os.str("");

            trajHelper.computeCautiousTrajectoryHeavy(tp->INIT_POINT, tp->INIT_CONTROL, T, fileName);
        }
        else if (typeTraj == STRATEGY_LIST) {
            os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1 << ".dat";
            fileName = os.str();
            os.str("");
                
            trajHelper.computeStrategyTrajectory(tp->INIT_POINT, T, fileName);
        }
        else {
            spdlog::warn("Unimplemeted trajectory type. No trajectory calculation.");
        }
	    }	
}

void ViabiBitSet::computeTychasticViableTrajectories() {
    ostringstream os;
    string fileName;

    const int nbTrajs = modelParams->getNbTrajectories();
    
    for (int tr = 0; tr < nbTrajs; tr++)
    {
        TrajectoryParametersManager tpm{modelParams, tr};
        const trajectoryParams *tp = tpm.getTrajectoryParameters();
        TypeTraj typeTraj = tp->TRAJECTORY_TYPE;
        double T = tp->maxTime;
        ViabiBitSetTrajectoryHelper trajHelper(grid, dynsys, &tpm);
        TychePicker tychePicker(&tpm, dynsys);

        if (typeTraj == STRATEGY_LIST) {
            os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1 << ".dat";
            fileName = os.str();
            os.str("");
                
            trajHelper.computeTychasticStrategyTrajectory(tp->INIT_POINT, T, fileName, tychePicker);
        }
        else {
            spdlog::error("Unimplemented trajectory type. No trajectory calculation");
        }        
    }
}

void ViabiBitSet::CaptureBasin()
    {
    dynsys->setDynamicsBackward();
    const algoViabiParams *avp = modelParams->getAlgoParameters();
    int refine = avp->GRID_REFINMENTS_NUMBER;

    ostringstream os;
    string fileName;
    if ((dynsys->getDynType() == CC) || (dynsys->getDynType() == DC))
	{

	CaptureBasin_ContinuousDynamics();
	}
    else
	{
	if ((dynsys->getDynType() == DD))
	    {
	    CaptureBasin_DiscreteDynamics();
	    }
	}

    this->saveViableSets();

    }

void ViabiBitSet::CaptureBasin_DiscreteDynamics()
    {

    }
void ViabiBitSet::CaptureBasin_ContinuousDynamics()
    {

    cout << " capture basin local rho \n";
    /*!
     * \var nbNewPoints : nb de nouveaux points ajout�s � l'�tape n
     */
    int nbNewPoints = 1;

    int iter = 0;
    /*!
     * On calcule la premi�re it�ration, en tanant compte de la convexification de la dynamique sur la cible
     * On appelle pour cela la fonction computeConvexifiedImage().
     */
    computeConvexifiedImage(iter);
    //     this->showCurrentImageList();
    ////////system("pause");
    /*!
     * -A partir de la liste des mailles on cr�e une liste ordonn�e de points repr�sentant l'image � l'aide de la fonction
     * createPointsList().
     */
    this->createPointsList();
    /*!
     * - A partir de la liste des points repr�sentant l'image  trois bases  de donn�es sont aliment�es: La base
     * principale, associ�e � la grille  et repr�sentant la fonction valeur, la base de r�troaction optimale et la base de
     *  r�tro-action viable. On appelle pour cela la fonction addNewPoints().
     */
    nbNewPoints = addNewPoints();
    cout << "  Premiers point points ajoutes =  " << nbNewPoints << endl;
    //system("pause");
    iter++;

    /*!
     * Tant qu'il y a de nouveaux points ajout�s on r�p�te les op�rations suivantes.
     */
    while ((nbNewPoints > 0))
	{
	cout << "  nbNewPoints=  " << nbNewPoints << endl;

	/*!
	 * -Calcul de l'image de \f^C_{n}\setminus C_n\f$  par la fonction computeCurrentImageLocalRho(): l'image est enregistr�e en �moire vive sous forme de liste
	 * de r�f�rences de mailles dans lesquelles arrive au moins une �volution. Chaque r�f�rence de maille contient
	 * des informations sur tous les ant�c�dants de cette maille ainsi  que la valeur minimale de
	 * temps. Dans cette version o� \f$ \rho\f$ est global, la fonction valeur prend la m�me valeur
	 * � chaque �tape : \f$ \rho \cdot n \f$.
	 */
	this->computeCurrentImage(iter);
	//this->showCurrentImageList();

	/*!
	 * -A partir de la liste des mailles on cr�e une liste ordonn�e de points repr�sentant l'image en appelant la fonction createPointsList().
	 *  Comme pour le mailles, chaque r�f�rence de point regroupe les informatons (regroup�es � partir de diff�rentes
	 * mailles dont est vertex)  sur les ant�c�dents  de ce point. Le but de la cr�ation de cette liste est d'�liminer
	 * les doublons afin de minimiser les acc�s � la base de donn�es
	 */
	this->createPointsList();
	//this->showCurrentImagePointsList();
	/*!
	 * - A partir de la liste des points repr�sentant l'image  trois bases  de donn�es sont aliment�es: La base
	 * principale, associ�e � la grille  et repr�sentant la fonction valeur, la base de r�troaction optimale et la base de
	 *  r�tro-action viable. On appelle ici la fonction addNewPoints().
	 */
	nbNewPoints = addNewPoints();
	// cout<< " points ajoutes =  "<<nbNewPoints<<endl;
	iter++;
	}

    }

void ViabiBitSet::computeCurrentImage(int iter)
    {
    spdlog::info("[Min Time problem] : Computing current image. number of points {}", currentPointsImage.size());
    int posX;

    list<imageCell>::iterator itStart, itNew;

    currentCellsImage.clear();
    double rho;
    map<unsigned long long int, double>::iterator itPoint = currentPointsImage.begin(), itLastPoint = currentPointsImage.end(), itTemp;
    map<unsigned long long int, list<unsigned long long int> >::iterator itCell, itLastCell, itCellsTemp;
    unsigned long long int intPointCoords[dim];
    double doublePointCoords[dim];
    DiscretPointImage *tempPointImage = new DiscretPointImage((unsigned long long int) 1, dim, dynsys, grid);
    map<unsigned long long int, list<unsigned long long int> > *cellsWithControls = tempPointImage->GetImageCellsWithControls();

    while (itPoint != itLastPoint)
	{
	posX = itPoint->first;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	rho = dynsys->calculRho_local(doublePointCoords);
	//spdlog::info("[Min Time problem] : Timestep =  {}", rho);

	/*!
	 * On calcule l'image discrète  du point
	 */
	double tempVal = itPoint->second + rho;

	tempPointImage->setPointNumAndRebuild(posX);

	itCell = cellsWithControls->begin();
	itLastCell = cellsWithControls->end();
	while (itCell != itLastCell)
	    {
	    unsigned long long int cellNum = itCell->first;
	    //cout<< "current im  cellnum = "<< cellNum<< " nb total de pojnts " << grid->getNbTotalPoints()<<" ";
	    if(cellNum < grid->getNbTotalPoints())
		{
		if (auto result = currentCellsImage.find(cellNum); result != currentCellsImage.end())
		    {
		    result->second = min((double) result->second, tempVal);
		    }
		else
		    {
		    currentCellsImage[cellNum] = tempVal;
		    }
		}
	    itCell++;
	    }
	itPoint++;
	}

    spdlog::info("[Min Time problem] : Computing current image finished. number of cells {}", currentCellsImage.size());
    currentPointsImage.clear();
    }
int ViabiBitSet::addNewPoints()
    {
    /*!
     * Cette fonction  réalise l'ajout  dans la base de données  de la nouvelle couche
     * \f$ \Phi(C_{n}\setminus C_{n-1})\f$  dans la base de données.
     * Elle doit également aouter aux bases de données de rétroaction optimale et viable
     * les rétroactions des points calculés
     *
     * Les points sont stockés sous forme de liste de structures imagePoint ordonnée par numéro de point x.
     *  Chaque structure contient le numéro du  point, la fonction valeur et
     *  les listes de triplets \f$ (x,u,\rho)\f$ représentant les rétrocations du point, optimale et viable.
     *
     *  Pour chaque point  de la liste on vérifie  d'abord s'il est déjé  dans la base de données. S'il n'y est pas,
     *  on l'ajoute é l'ensemble en construction et en méme temps on ajoute ses rétroactions aux deux bases de rétroaction
     *  correspondantes.
     *
     *  Si le point existe déjé dans l'ensemble, il posséde déjé une rétro-action égaement. Dans ce cas, on
     *  vérifie les fonctions valeurs. Si celle du nuveau point est inférieure, alors on modifie le point existant
     *  et la rétroaction optimale. Sinon,  on modifie  selement les rétroactons viables.
     *
     *  Important!  Il est é noter ici que c'est cette méme liste de points qui doit étre ajoutée  é la base
     *  par cette fonction qui servira ensuite é la construction de la couche suivante.
     *  Ainsi, si un point de cette liste exste déjé dans la base et que sa fonction valeur de la base n'est pas modifiée,
     *  il n'est pas considéré comme nouveau, il n'appartient pas é \f$ C_{n+1}\setminus C_n\f$ . On doit donc le supprimer
     *   de la liste pour ne pas  calculer son image plutard.
     *
     */

    double t1, t2, elapsed_time;

    timeval tim;
    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

    int nbNewPoints = 0;
    unsigned long long int intCoords[dim];
    map<unsigned long long int, double>::iterator itPoint = currentPointsImage.begin(), itLastPoint = currentPointsImage.end(), itTemp;

    imagePoint tempPoint;
    while (itPoint != itLastPoint)
	{
	unsigned long long int pointNum = itPoint->first;
	double value = itPoint->second;

	grid->numToIntCoords(pointNum, intCoords);
	if (!grid->isInSet(intCoords))
	    {

	    // cout << " le point est nouveau on l'ajout dans la base  de fonc val \n";
	    grid->addPointToSet(pointNum, 1.0);
	    //     cout<< " on l'ajoute maintenant  dans la base  de retro \n *******************************\n";

	    /*
	     * Ajout de r�troactions  dans les bases de donn�es correspondantes
	     */
	    itPoint++;
	    nbNewPoints++;
	    }
	else
	    {
	    itPoint = currentPointsImage.erase(itPoint);
	    }
	}

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));

    cout << "Elapsed time : " << elapsed_time << " sec." << endl << std::flush;
    return nbNewPoints;
    }

void ViabiBitSet::noyauViabi_sansControle(int nbArret)
    {

    int dirTramage = grid->getDirTram();

    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    bool testF;

    //    cout<<"masque points enleves cree"<<masquePointsEnleves;
    unsigned long long int longTrame = grid->getLongTrame();
    boost::dynamic_bitset<> masque;
    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(longTrame, 0);
    boost::dynamic_bitset<> **gridTab = grid->getGridTab();

    //pointeur sur la fonction d'analyse d'un point
    //permet de changer de fonction selon la dimension de l'espace de commandes

    double *limInf = grid->limInf;
    double *gridStep = grid->step;

    unsigned long long int indice[dim];

    double xCoordsDouble[dim];
    bool testNonVide;

    int nbIter = 0;
    unsigned long long int comptEtats = 0, comptEnleves = nbArret + 1;
    unsigned long long int subGridSize = grid->getNbPointsTotalSubGrid();
    testK0();

    while (comptEnleves > (unsigned long long int) nbArret)
	{
	cout << "nouvelle boucle while\n";

	comptEnleves = 0;
	comptEtats = 0;

	for (unsigned long long int posX = 0; posX < subGridSize; posX++)
	    {
	    //cout<< " posX="<<posX<<endl;
	    if (!gridTab[posX]->none())
		{
		comptEtats++;
		testF = false;
		numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

		//cout<<  (*gridTab[posX])<<endl;

		masque = grid->analyseTrameMasqueBis(posX, 0);

		//	cout<<"masque d'analyse  "<<masque<<endl;

		masquePointsEnleves->set();

		if (masque.none() | testF)
		    {
		    //		cout<<" rien e analyser posx= "<<posX<<" \n";
		    }
		else
		    {
		    for (int j = 0; j < dirTramage; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j] * gridStep[j];
			}

		    for (int j = dirTramage + 1; j < dim; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j - 1] * gridStep[j];
			}
		    ///cout<< " Analyse en cours "<<endl;
		    for (unsigned long long int k = 0; k < longTrame; k++)
			{
			if (masque[k])
			    {
			    xCoordsDouble[dirTramage] = limInf[dirTramage] + k * gridStep[dirTramage];
			    testNonVide = this->findViabImagePoint_noControl(xCoordsDouble, false);
			    //cout<<" "<<testNonVide;
			    if (!testNonVide)
				{
				masquePointsEnleves->set(k, false);
				comptEnleves++;
				}
			    }	// fin de if masque[k]
			}	// fin de for  de parcours de masque
		    //cout<< " fini"<<endl;
		    //cout<< " masque points enleves "<<*masquePointsEnleves<<endl;
		    if (masquePointsEnleves->count() < (unsigned long long int) longTrame)
			{

			*gridTab[posX] &= (*masquePointsEnleves);
			}

		    }

		}			//fin de if la trame n'est pas vide

	    }				// fin de for de parcours de la trame

    cout << "End of iteration " << nbIter << " Number of points removed: " << comptEnleves << "\n";

	nbIter++;

	}
    cout << "fini nbIter=" << nbIter;
    //foncCarNoyau->printTrame();

    }

void ViabiBitSet::noyauViabi_sansControle_omp(int nbArret)
    {

    boost::dynamic_bitset<> **gridTab = grid->getGridTab();
    boost::dynamic_bitset<> **gridTabNew = grid->getGridTabNew();
    unsigned long long int subGridSize = grid->getNbPointsTotalSubGrid();
    unsigned long long int comptEnleves = nbArret + 1;

    testK0();
    grid->copyGrid(gridTab, gridTabNew);

    int nbIter = 0;

    while (comptEnleves > (unsigned long long int) nbArret)
	{
	comptEnleves = 0;

	//cout<<"dim etat>2\n";

	// cout<<"dim etat  " <<dim<< " taille de trame est " <<tailleTrame<<"\n";
	unsigned long long int posX = 0;
	double *xCoordsDouble;
	unsigned long long int *indice;
#pragma omp parallel  num_threads(nbOMPThreads)  reduction(+:comptEnleves) private(posX, xCoordsDouble, indice)  shared( gridTab, gridTabNew,subGridSize) default(none)
	{
	xCoordsDouble = new double[dim];

	indice = new unsigned long long int[dim];
#pragma omp for
	for (posX = 0; posX < subGridSize; posX++)
	    {
	    int tid = omp_get_thread_num();
	    int dirTramage = grid->getDirTram();

	    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

	    unsigned long long int longTrame = grid->getLongTrame();

	    double *limInf = grid->limInf;
	    double *gridStep = grid->step;

	    bool testNonVide;
	    //printf("thread numero %d nouvelle boucle while\n", tid);
	    boost::dynamic_bitset<> masque;
	    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(longTrame, 0);

	    //  cout<< " posX="<<posX<< " size "<<gridTab->size()<<endl;
	    if (!gridTab[posX]->none())
		{
		numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);
		masque = grid->analyseTrameMasqueBis(posX, 1 - tid);

		masquePointsEnleves->set();
		if (!masque.none())
		    {
		    for (int j = 0; j < dirTramage; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j] * gridStep[j];
			}

		    for (int j = dirTramage + 1; j < dim; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j - 1] * gridStep[j];
			}
		    for (unsigned long long int k = 0; k < longTrame; k++)
			{
			if (masque[k])
			    {
			    xCoordsDouble[dirTramage] = limInf[dirTramage] + k * gridStep[dirTramage];
			    testNonVide = this->findViabImagePoint_noControl(xCoordsDouble, false);

			    if (!testNonVide)
				{

				masquePointsEnleves->set(k, false);
				comptEnleves++;
				}
			    }			// fin de if masque[k]
			}		// fin de for  de parcours de masque

		    if (masquePointsEnleves->count() < (unsigned long long int) longTrame)
			{
			(*gridTabNew[posX]) &= (*masquePointsEnleves);
			}
		    }
		}			//fin de if la trame n'est pas vide

	    }					//fin de for OMP
	delete[] xCoordsDouble;
	delete[] indice;
	}
	grid->copyGrid(gridTabNew, gridTab);

	cout << "Itération " << nbIter << " terminée. Nombre de points  points enlevés: " << comptEnleves << "\n";
	nbIter++;

	}
    cout << "fini nbIter=" << nbIter;
    }

void ViabiBitSet::GarantedViabilityKernel(int nbArret)
    {
    dynsys->setDynamicsForward();
    const algoViabiParams *avp = modelParams->getAlgoParameters();
    int refine = avp->GRID_REFINMENTS_NUMBER;

    ostringstream os;
    string fileName;
    double t1, t2, elapsed_time, t1_glob, t2_glob, elapsed_time_glob;
    timeval tim, tim_glob;
    //  raffinements successifs ( y compris le premier calcul si refine=0)
    int refIter = -1;
    int seuilArret = nbArret;

    cout << " debut de la boucle de rafinements nb refine = " << refine << endl;
    while (refIter < refine)
	{
	gettimeofday(&tim, NULL);              //mesure le temps d'execution
	t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	//  Calcul du noyau de viabilité
	cout << "=============================Debut viab kernel ==============================" << endl;
	GuarantedViabilityKernelSimple(seuilArret);
	cout << "==============================================================================" << endl;
	gettimeofday(&tim, NULL);
	t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	elapsed_time = (double) ((t2 - t1));
	cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;
	//   cout<< "=============================================================================="<<endl;
	/*  * sauvegardes
	 */
	if (avp->INTERMEDIATE_SAVINGS)
	    {
	    this->saveIntermediateViableSets(refIter + 1);
	    }

	refIter++;
	seuilArret *= 2;
	/*
	 * raffinement du maillage
	 */

	if (refIter < refine)
	    {
	    grid->refine();
	    }
	}

    this->saveViableSets();
    }

SysDyn* ViabiBitSet::GetSysDynForViabProblem()
    {
    return this->dynsys;
    }

string ViabiBitSet::getSetName(SetType type)
{
    string setName = toString(type);
    transform(setName.begin(), setName.end(), setName.begin(), [](char c) {
        return tolower(c);
    });
    // On veut que VIABG renvoie "viabG" (et non pas "viabg")
    if (type == VIABG) setName[4] = 'G';
    return setName;
}

void ViabiBitSet::loadViableSets()
    {

        const algoViabiParams *avp = modelParams->getAlgoParameters();
    int refine = avp->GRID_REFINMENTS_NUMBER;

    ostringstream os;
    string fileName;
    string setName = getSetName(avp->SET_TYPE);

    // on charge dans la mémoire l'ensemble calculé et enregistré
    // correspondant au dernier raffinement
    os << "../OUTPUT/" << filePrefix << "-" << setName << ".dat";
    fileName = os.str();
    os.str("");
    grid->loadSet(fileName);

    }

void ViabiBitSet::saveViableSets(const string &baseFilenameSuffix) {
    ostringstream os;
    string fileName;
    const algoViabiParams *avp = modelParams->getAlgoParameters();
    string setName = getSetName(avp->SET_TYPE);

    os << "../OUTPUT/" << filePrefix << "-" << setName << baseFilenameSuffix << ".dat";
    fileName = os.str();
    os.str("");
    if (avp->SAVE_VIAB_LIGHT)
	{
        grid->saveValOnGridLight(fileName);
	}
    else
	{
        grid->saveValOnGrid(fileName);
	}
    if (avp->SAVE_SLICE)
	{
        os << "../OUTPUT/" << filePrefix << "-" << setName << baseFilenameSuffix << "-Slice.dat";
        fileName = os.str();
        os.str("");
        grid->saveCoupe(fileName);
	}
    if (avp->SAVE_SLICE_BOUND)
	{
        os << "../OUTPUT/" << filePrefix << "-" << setName << baseFilenameSuffix << "-SliceBound.dat";
        fileName = os.str();
        os.str("");
        grid->saveCoupeBoundary(fileName);
	}
    if (avp->SAVE_BOUNDARY)
	{
        os << "../OUTPUT/" << filePrefix << "-" << setName << baseFilenameSuffix << "-bound.dat";
        fileName = os.str();
        os.str("");
        grid->saveBoundary(fileName);
	}

    if (avp->SAVE_PROJECTION)
	{
        os << "../OUTPUT/" << filePrefix << "-" << setName << baseFilenameSuffix << "-proj.dat";
        fileName = os.str();
        os.str("");
        /*
         *  calcul et sauvegarde  de la projection du  noyau
         */
        grid->saveProjetion(fileName, avp->PROJECTION);
	}
}

void ViabiBitSet::saveViableSets()
{
    saveViableSets("");
}

void ViabiBitSet::saveIntermediateViableSets(int refine)
{
    std::ostringstream os;
    os << "-" << refine;
    saveViableSets(os.str());
}

void ViabiBitSet::ViabilityKernel(int nbArret)
    {
    dynsys->setDynamicsForward();
    const algoViabiParams *avp = modelParams->getAlgoParameters();
    int refine = avp->GRID_REFINMENTS_NUMBER;

    ostringstream os;
    string fileName;
    double t1, t2, elapsed_time, t1_glob, t2_glob, elapsed_time_glob;
    timeval tim, tim_glob;
    //  raffinements successifs ( y compris le premier calcul si refine=0)
    int refIter = -1;
    int seuilArret = nbArret;

    cout << " debut de la boucle de rafinements nb refine = " << refine << endl;
    while (refIter < refine)
	{
	gettimeofday(&tim, NULL);              //mesure le temps d'execution
	t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	//  Calcul du noyau de viabilité
	cout << "=============================Debut viab kernel ==============================" << endl;
	ViabilityKernelSimple(seuilArret);
	cout << "==============================================================================" << endl;
	gettimeofday(&tim, NULL);
	t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	elapsed_time = (double) ((t2 - t1));
	cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;
	//   cout<< "=============================================================================="<<endl;
	/*  * sauvegardes
	 */
	if (avp->INTERMEDIATE_SAVINGS)
	    {
	    this->saveIntermediateViableSets(refIter + 1);
	    }

	refIter++;
	seuilArret *= 2;
	/*
	 * raffinement du maillage
	 */

	if (refIter < refine)
	    {
	    grid->refine();
	    }
	}

    this->saveViableSets();
    }

void ViabiBitSet::ViabilityKernelSimple(int nbArret)
{
    if (dynsys->getDimC() == 0)
	{
        if (nbOMPThreads > 1)
	    {
            noyauViabi_sansControle_omp(nbArret);
	    }
        else
	    {
            cout << "viab kernel sans controle" << endl;
            noyauViabi_sansControle(nbArret);
	    }
	}
    else
	{
        if ((dynsys->getDynType() == CC) || (dynsys->getDynType() == DC))
	    {
            if (nbOMPThreads > 1)
            {
                noyauViabi_omp(nbArret);
            }
            else
            {
                noyauViabi(nbArret);
            }
	    }
        else
	    {
            if ((dynsys->getDynType() == DD))
            {
                noyauViabi_FD(nbArret);
            }
	    }
	}
}

void ViabiBitSet::GuarantedViabilityKernelSimple(int nbArret)
    {

    if ((dynsys->getDynType() == CC) || (dynsys->getDynType() == DC))
	{

	noyauViabiGuaranti(nbArret);
	}
    else
	{
	if ((dynsys->getDynType() == DD))
	    {
	    noyauViabiGaranti_FD(nbArret);
	    }
	}
    }

void ViabiBitSet::noyauViabi(int nbArret)
    {

    int dirTramage = grid->getDirTram();

    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    bool testF;

    unsigned long long int longTrame = grid->getLongTrame();
    boost::dynamic_bitset<> masque;
    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(longTrame, 0);
    boost::dynamic_bitset<> **gridTab = grid->getGridTab();

    //pointeur sur la fonction d'analyse d'un point
    //permet de changer de fonction selon la dimension de l'espace de commandes

    double *limInf = grid->limInf;
    double *gridStep = grid->step;

    unsigned long long int indice[dim];

    double xCoordsDouble[dim];
    bool testNonVide;
    int nbIter = 0;
    unsigned long long int comptEtats = 0, comptEnleves = nbArret + 1;
    unsigned long long int subGridSize = grid->getNbPointsTotalSubGrid();

    testK0();

    while (comptEnleves > (unsigned long long int) nbArret)
	{
	cout << "nouvelle boucle while\n";

	comptEnleves = 0;

	//cout<<"dim etat>2\n";
	comptEtats = 0;
	for (unsigned long long int posX = 0; posX < subGridSize; posX++)
	    {
	    //        cout<< " posX="<<posX<< " size "<<gridTab->size()<<endl;
	    if (!gridTab[posX]->none())
		{
		comptEtats++;
		testF = false;
		numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

		// cout<<  (*gridTab[posX])<<endl;

		masque = grid->analyseTrameMasque(posX);

		// cout<<"masque d'analyse  "<<masque<<endl;

		masquePointsEnleves->set();

		if (masque.none() | testF)
		    {
		    cout << " rien e analyser posx= " << posX << " \n";
		    }
		else
		    {
		    for (int j = 0; j < dirTramage; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j] * gridStep[j];
			}

		    for (int j = dirTramage + 1; j < dim; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j - 1] * gridStep[j];
			}

		    for (unsigned long long int k = 0; k < longTrame; k++)
			{
			if (masque[k])
			    {
			    xCoordsDouble[dirTramage] = limInf[dirTramage] + k * gridStep[dirTramage];
			    testNonVide = this->findViabImagePoint(xCoordsDouble, false);

			    if (!testNonVide)
				{
				masquePointsEnleves->set(k, false);
				comptEnleves++;
				}
			    }	// fin de if masque[k]
			}	// fin de for  de parcours de masque

		    if (masquePointsEnleves->count() < (unsigned long long int) longTrame)
			{
			//cout<< " on enleve sur posX = "<<posX<<endl;
			//cout<< " grid tab ici    "<< *gridTab[posX] << endl;
			//cout<< " masque          "	<< 	*masquePointsEnleves <<endl;
			*gridTab[posX] &= (*masquePointsEnleves);
			//cout<< " grid tab  APRES "<< *gridTab[posX] << endl;
			//cout<< " =========================================================\n";
			}

		    }

		}			//fin de if la trame n'est pas vide

	    }				// fin de for de parcours de la trame

	cout << "Itération " << nbIter << " terminée. Nombre de points  points enlevés: " << comptEnleves << "\n";

	nbIter++;

	}
    cout << "fini nbIter=" << nbIter;
    //foncCarNoyau->printTrame();

    }

void ViabiBitSet::noyauViabiGuaranti(int nbArret)
    {

    int dirTramage = grid->getDirTram();

    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    bool testF;

    unsigned long long int longTrame = grid->getLongTrame();
    boost::dynamic_bitset<> masque;
    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(longTrame, 0);
    boost::dynamic_bitset<> **gridTab = grid->getGridTab();

    //pointeur sur la fonction d'analyse d'un point
    //permet de changer de fonction selon la dimension de l'espace de commandes

    double *limInf = grid->limInf;
    double *gridStep = grid->step;

    unsigned long long int indice[dim];

    double xCoordsDouble[dim];
    bool testNonVide;
    int nbIter = 0;
    unsigned long long int comptEtats = 0, comptEnleves = nbArret + 1;
    unsigned long long int subGridSize = grid->getNbPointsTotalSubGrid();

    testK0();

    while (comptEnleves > (unsigned long long int) nbArret)
	{
	cout << "nouvelle boucle while\n";

	comptEnleves = 0;

	//cout<<"dim etat>2\n";
	comptEtats = 0;
	for (unsigned long long int posX = 0; posX < subGridSize; posX++)
	    {
	    //        cout<< " posX="<<posX<< " size "<<gridTab->size()<<endl;
	    if (!gridTab[posX]->none())
		{
		comptEtats++;
		testF = false;
		numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

		// cout<<  (*gridTab[posX])<<endl;

		masque = grid->analyseTrameMasque(posX);

		// cout<<"masque d'analyse  "<<masque<<endl;

		masquePointsEnleves->set();

		if (masque.none() | testF)
		    {
		    cout << " rien e analyser posx= " << posX << " \n";
		    }
		else
		    {
		    for (int j = 0; j < dirTramage; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j] * gridStep[j];
			}

		    for (int j = dirTramage + 1; j < dim; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j - 1] * gridStep[j];
			}

		    for (unsigned long long int k = 0; k < longTrame; k++)
			{
			if (masque[k])
			    {
			    xCoordsDouble[dirTramage] = limInf[dirTramage] + k * gridStep[dirTramage];
			    testNonVide = this->findGuarantedViabImagePoint(xCoordsDouble, false);

			    if (!testNonVide)
				{
				masquePointsEnleves->set(k, false);
				comptEnleves++;
				}
			    }	// fin de if masque[k]
			}	// fin de for  de parcours de masque

		    if (masquePointsEnleves->count() < (unsigned long long int) longTrame)
			{
			//cout<< " on enleve sur posX = "<<posX<<endl;
			//cout<< " grid tab ici    "<< *gridTab[posX] << endl;
			//cout<< " masque          "	<< 	*masquePointsEnleves <<endl;
			*gridTab[posX] &= (*masquePointsEnleves);
			//cout<< " grid tab  APRES "<< *gridTab[posX] << endl;
			//cout<< " =========================================================\n";
			}

		    }

		}			//fin de if la trame n'est pas vide

	    }				// fin de for de parcours de la trame

	cout << "Itération " << nbIter << " terminée. Nombre de points  points enlevés: " << comptEnleves << "\n";

	nbIter++;

	}
    cout << "fini nbIter=" << nbIter;
    //foncCarNoyau->printTrame();

    }

void ViabiBitSet::noyauViabi_omp(int nbArret)
    {
    boost::dynamic_bitset<> **gridTab = grid->getGridTab();
    boost::dynamic_bitset<> **gridTabNew = grid->getGridTabNew();
    unsigned long long int subGridSize = grid->getNbPointsTotalSubGrid();
    unsigned long long int comptEnleves = nbArret + 1;

    testK0();
    grid->copyGrid(gridTab, gridTabNew);

    int nbIter = 0;

    while (comptEnleves > (unsigned long long int) nbArret)
	{
	comptEnleves = 0;

	//cout<<"dim etat>2\n";

	// cout<<"dim etat  " <<dim<< " taille de trame est " <<tailleTrame<<"\n";

	unsigned long long int posX = 0;
	double *xCoordsDouble;
	unsigned long long int *indice;
#pragma omp parallel  num_threads(nbOMPThreads)  reduction(+:comptEnleves) private(posX, xCoordsDouble, indice)  shared( gridTab, gridTabNew,subGridSize) default(none)
	{
	xCoordsDouble = new double[dim];
	indice = new unsigned long long int[dim];
#pragma omp for
	for (posX = 0; posX < subGridSize; posX++)
	    {
	    //int  tid = omp_get_thread_num();

	    int dirTramage = grid->getDirTram();

	    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

	    unsigned long long int longTrame = grid->getLongTrame();

	    double *limInf = grid->limInf;
	    double *gridStep = grid->step;
	    // printf("thread numero %d nouvelle boucle while\n", nbTh);

	    if (!gridTab[posX]->none())
		{
		// cout<<  (*gridTab[posX])<<endl;

		boost::dynamic_bitset<> masque = grid->analyseTrameMasqueBis(posX, false);

		//   cout<<"masque d'analyse  "<<masque<<endl;

		if (!masque.none())
		    {
		    bool testNonVide;

		    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(longTrame, 0);
		    numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

		    masquePointsEnleves->set();

		    for (int j = 0; j < dirTramage; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j] * gridStep[j];
			}

		    for (int j = dirTramage + 1; j < dim; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j - 1] * gridStep[j];
			}

		    for (unsigned long long int k = 0; k < longTrame; k++)
			{

			if (masque[k])
			    {
			    xCoordsDouble[dirTramage] = limInf[dirTramage] + k * gridStep[dirTramage];
			    testNonVide = this->findViabImagePoint(xCoordsDouble, false);

			    //cout<< " fini\n";
			    if (!testNonVide)
				{
				masquePointsEnleves->set(k, false);
				comptEnleves++;
				}
			    }			// fin de if masque[k]
			}		// fin de for  de parcours de masque

		    if (masquePointsEnleves->count() < (unsigned long long int) longTrame)
			{
			(*gridTabNew[posX]) &= (*masquePointsEnleves);
			}
		    }
		}			//fin de if la trame n'est pas vide

	    }						//fin de for OMP
	delete[] xCoordsDouble;
	delete[] indice;
	}
	//cout<<"Itération "<<nbIter<< " AVANT COPIE "<<comptEnleves<<"\n";
	grid->copyGrid(gridTabNew, gridTab);

	cout << "Itération " << nbIter << " terminée. Nombre de points  points enlevés: " << comptEnleves << "\n";

	nbIter++;

	}
    cout << "Calcul de noyau fini. Nb total d'itérations: " << nbIter;

    }

void ViabiBitSet::noyauViabiGaranti_FD(int nbArret)
    {
    int dirTramage = grid->getDirTram();

    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    // vrai si  pendant l'iteration en cours on a enleve au moins un point du noyau en construction
    //parcours de trame monodirectionnel

    int tailleTrame = grid->getNbPointsTotalSubGrid();
    //cout<<"dim etat  " <<dim<< " taille de trame est " <<tailleTrame<<"\n";

    //calcul de la taille e prevoir pour les coordonnees des indices de debut de parcours
    //que la methode GPU va renvoyer

    bool testF;

    //	cout<<"masque points enleves cree"<<masquePointsEnleves;
    unsigned long long int longTrame = grid->getLongTrame();

    unsigned long long int nbC = dynsys->getTotalNbPointsC();

    unsigned long long int nbTy = dynsys->getTotalNbPointsTy();

    // cout<< "  noayu viab garanti FD bitset nbTy= "<<nbTy<<endl;

    boost::dynamic_bitset<> masque;
    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(longTrame, 0);
    boost::dynamic_bitset<> **gridTab = grid->getGridTab();
    //pointeur sur la fonction d'analyse d'un point
    //permet de changer de fonction selon la dimension de l'espace de commandes

    unsigned long long int indice[dim];

    unsigned long long int **controlIntCoords = dynsys->getControlIntCoords();
    unsigned long long int **tychIntCoords = dynsys->getTychIntCoords();

    unsigned long long int coordDiscretes[dim], imageXU[dim];

    unsigned long long int compteComm, compteTych;
    bool testNonVide, testAllTych;
    int nbIter = 0;
    int comptEtats = 0, comptEnleves = nbArret + 1;

    while (comptEnleves > nbArret)
	{
	cout << "nouvelle boucle while\n";

	comptEnleves = 0;

	//cout<<"dim etat>2\n";
	comptEtats = 0;
	//cout<<"dim etat  " <<dim<< " taille de trame est " <<tailleTrame<<"\n";
	for (int posX = 0; posX < tailleTrame; posX++)
	    {
	    // cout<< " posX="<<posX<<endl;
	    if (!(*gridTab[posX]).none())
		{
		comptEtats++;
		testF = false;
		numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);
		//cout<<  (*gridTab[posX])<<endl;

		masque = grid->analyseTrameMasque(posX);

		//  	cout<<"masque d'analyse  "<<masque<<endl;

		masquePointsEnleves->set();

		if (masque.none() | testF)
		    {
		    //cout<<" rien e analyser posx= "<<posX<<" \n";
		    }
		else
		    {
		    for (int j = 0; j < dirTramage; j++)
			{
			coordDiscretes[j] = indice[j];

			}

		    for (int j = dirTramage + 1; j < dim; j++)
			{
			coordDiscretes[j] = indice[j - 1];

			}

		    if (!((masque.count() == 2) && (masque[0] && masque[longTrame - 1])))
			{

			for (unsigned long long int k = 0; k < longTrame; k++)
			    {

			    if (masque[k])
				{

				coordDiscretes[dirTramage] = k;
				//cout<< "  coords discretes du point ";
				// printVector(coordDiscretes,dim);
				testNonVide = false;
				compteComm = 0;

				// cout<<"\n analyse de point nbC= "<<nbC<<endl;
				//	cout<<" test non vide =  ";
				while (!testNonVide && compteComm < nbC)
				    {
				    // calcul  de vecteur de controle u

				    // calcul de l'image F(x,u)
				    /*cout<< "  control coords";
				     for(int jj=0;jj<3;jj++)
				     cout<< " "<<controlIntCoords[compteComm][jj];
				     cout<<endl;*/
				    if (dynsys->constraintsXU_fd(coordDiscretes, controlIntCoords[compteComm]) < PLUS_INF)
					{
					// cout<< "  controle "<<controlIntCoords[compteComm][0]<<" "<<controlIntCoords[compteComm][1]<< " est autorise \n";

					testAllTych = true;
					for (compteTych = 0; compteTych < nbTy; compteTych++)
					    {
					    //   cout<< "  tychet "<<tychIntCoords[compteTych][0]<<endl;
					    dynsys->dynamics_tych_fd(coordDiscretes, controlIntCoords[compteComm], tychIntCoords[compteTych],
						    imageXU);
					    // cout<<  " image ";
					    //                                             printVector(imageXU,dim);

					    if (grid->isPointInGrid_fd(imageXU))
						{
						testAllTych &= grid->isInSet(imageXU);

						//   cout<< "image est dans la grille de calcul ok et test no vide est "<<testAllTych<<endl;
						}
					    else
						{
						testAllTych = false;
						// cout<< "  on sort de la grille \n";
						}
					    }
					testNonVide = testAllTych;
					}

				    //	cout<<"   "<<testNonVide;
				    compteComm++;
				    }
				//cout<< " fini\n";
				if (!testNonVide)
				    {

				    masquePointsEnleves->set(k, false);
				    comptEnleves++;
				    }
				}			// fin de if masque[k]
			    }		// fin de for  de parcours de masque

			if (masquePointsEnleves->count() < (unsigned long long int) longTrame)
			    {

			    (*gridTab[posX]) &= (*masquePointsEnleves);
			    }
			}
		    else
			{
			int k = 0;
			coordDiscretes[dirTramage] = k;

			testNonVide = false;
			compteComm = 0;

			//	cout<<"\n analyse de point nbC= "<<nbC<<endl;
			// 	cout<<" test non vide =  ";
			while (!testNonVide && compteComm < nbC)
			    {
			    if (dynsys->constraintsXU_fd(coordDiscretes, controlIntCoords[compteComm]) < PLUS_INF)
				{
				testAllTych = true;
				for (compteTych = 0; compteTych < nbTy; compteTych++)
				    {
				    dynsys->dynamics_tych_fd(coordDiscretes, controlIntCoords[compteComm], tychIntCoords[compteTych], imageXU);
				    if (grid->isPointInGrid_fd(imageXU))
					{
					testAllTych &= grid->isInSet(imageXU);
					//cout<< " point est dans la grille de calcul ok et test no vide est "<<testNonVide<<endl;
					}
				    else
					{
					testAllTych = false;
					//cout<< "  on sort de la grille \n";
					}

				    }
				testNonVide = testAllTych;
				}
			    //	cout<<"   "<<testNonVide;
			    compteComm++;
			    }
			if (!testNonVide)
			    {
			    masquePointsEnleves->set(k, false);
			    comptEnleves++;
			    }
			k = longTrame - 1;
			coordDiscretes[dirTramage] = k;

			testNonVide = false;
			compteComm = 0;

			while (!testNonVide && compteComm < nbC)
			    {
			    if (dynsys->constraintsXU_fd(coordDiscretes, controlIntCoords[compteComm]) < PLUS_INF)
				{
				testAllTych = true;
				for (compteTych = 0; compteTych < nbTy; compteTych++)
				    {
				    dynsys->dynamics_tych_fd(coordDiscretes, controlIntCoords[compteComm], tychIntCoords[compteTych], imageXU);
				    if (grid->isPointInGrid_fd(imageXU))
					{
					testAllTych &= grid->isInSet(imageXU);
					//cout<< " point est dans la grille de calcul ok et test no vide est "<<testNonVide<<endl;
					}
				    else
					{
					testAllTych = false;
					//cout<< "  on sort de la grille \n";
					}
				    }
				testNonVide = testAllTych;
				}
			    //	cout<<"   "<<testNonVide;
			    compteComm++;
			    }
			if (!testNonVide)
			    {
			    masquePointsEnleves->set(k, false);
			    comptEnleves++;
			    }

			if (masquePointsEnleves->count() < (unsigned long long int) longTrame)
			    {
			    (*gridTab[posX]) &= (*masquePointsEnleves);
			    }
			}
		    }

		}			//fin de if la trame n'est pas vide
	    }				// fin de for de parcours de la trame

	nbIter++;
	//cout<<"  points enleves: "<<comptEnleves<<"\n";
	//cout<<" nbIter="<<nbIter<<"\n";

	}
    cout << "Calcul fini. Nombre total d'itérations : " << nbIter;

    }

void ViabiBitSet::saveViabRetro(string fileName)
    {
    //cout<<"ecriture  de l'ensemble dans un fichier \n";
    // instructions

    ofstream fichierB(fileName.c_str());

    if (fichierB)  // si l'ouverture a r�ussi
	{

	int dirTramage = grid->getDirTram();

	const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

	unsigned long long int tailleTrame = grid->getNbPointsTotalSubGrid();
	//	cout<<"masque points enleves cree"<<masquePointsEnleves;
	unsigned long long int longTrame = grid->getLongTrame();

	unsigned long long int nbC = dynsys->getTotalNbPointsC();
	boost::dynamic_bitset<> masque;
	boost::dynamic_bitset<> **gridTab = grid->getGridTab();

	unsigned long long int indice[dim - 1];

	unsigned long long int **controlIntCoords = dynsys->getControlIntCoords();

	unsigned long long int coordDiscretes[dim], imageXU[dim];

	unsigned long long int compteComm;
	bool testNonVide;
	int nbIter = 0;
	vector<unsigned long long int> viabControls;
	vector<unsigned long long int> images;
	unsigned long long int posIm;

	for (unsigned long long int posX = 0; posX < tailleTrame; posX++)
	    {
	    if (!(*gridTab[posX]).none())
		{
		numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

		for (int j = 0; j < dirTramage; j++)
		    {
		    coordDiscretes[j] = indice[j];
		    }

		for (int j = dirTramage + 1; j < dim; j++)
		    {
		    coordDiscretes[j] = indice[j - 1];
		    }

		for (unsigned long long int k = 0; k < longTrame; k++)
		    {

		    if ((*gridTab[posX])[k])
			{
			viabControls.clear();
			images.clear();
			coordDiscretes[dirTramage] = k;

			//  cout<< " point x "; printVector(coordDiscretes,dim);
			// system("pause");
			posIm = grid->intCoordsToNum(coordDiscretes);
			// cout<< "  la position du point  x "<< posIm<<endl;
			fichierB << posIm << " ";
			testNonVide = false;
			compteComm = 0;
			//	cout<<"\n analyse de point nbC= "<<nbC<<endl;
			//	cout<<" test non vide =  ";
			while (compteComm < nbC)
			    {
			    // calcul  de vecteur de controle u

			    // calcul de l'image F(x,u)

			    /* cout<< "  control coords";
			     for(int jj=0;jj<3;jj++)
			     cout<< " "<<controlIntCoords[compteComm][jj];
			     cout<<endl;*/
			    if (dynsys->constraintsXU_fd(coordDiscretes, controlIntCoords[compteComm]) < PLUS_INF)
				{
				//cout<< "  controle "<<controlIntCoords[compteComm][0]<<" "<<controlIntCoords[compteComm][1]<< " est autorise \n";

				dynsys->dynamics_fd(coordDiscretes, controlIntCoords[compteComm], imageXU);
				// cout<< " image "; printVector(imageXU,dim);
				//cout<<  " image ";
				//printVector(imageXU,dim);
				if (grid->isPointInGrid_fd(imageXU))
				    {
				    testNonVide = grid->isInSet(imageXU);
				    //cout<< " point est dans la grille de calcul ok et test no vide est "<<testNonVide<<endl;
				    }
				else
				    {
				    testNonVide = false;
				    //cout<< "  on sort de la grille \n";
				    }

				}
			    else
				{
				testNonVide = false;
				}
			    //	 cout<<"   "<<testNonVide;
			    if (testNonVide)
				{
                    posIm = grid->intCoordsToNum(imageXU);
				//	 	cout<< " numero image "<<posIm<<endl;
				viabControls.push_back(compteComm);
				images.push_back(posIm);
				}
			    compteComm++;
			    }
			// system("pause");
			fichierB << images.size() << " ";
			for (unsigned long long int jj = 0; jj < images.size(); jj++)
			    {
			    fichierB << viabControls.at(jj) << " " << images.at(jj) << " ";
			    }
			fichierB << endl;
			//cout<< " fini\n";

			}				// fin de if masque[k]
		    }			// fin de for  de parcours de masque
		//system("pause");

		}			//fin de if la trame n'est pas vide
	    }				// fin de for de parcours de la trame

	nbIter++;

	fichierB.close();
	// je referme le fichier

	}
    else
	// sinon
	cerr << "Erreur � l'ouverture !" << endl;

    //cout<<"fichier fini\n";

    }

void ViabiBitSet::saveViabGarantiRetro(string fileName)
    {
    //cout<<"ecriture  de l'ensemble dans un fichier \n";
    // instructions

    ofstream fichierB(fileName.c_str());

    if (fichierB)  // si l'ouverture a r�ussi
	{
	int dirTramage = grid->getDirTram();

	const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

	int tailleTrame = grid->getNbPointsTotalSubGrid();
	//	cout<<"masque points enleves cree"<<masquePointsEnleves;
	unsigned long long int longTrame = grid->getLongTrame();

	unsigned long long int nbC = dynsys->getTotalNbPointsC();
	unsigned long long int nbTy = dynsys->getTotalNbPointsTy();

	boost::dynamic_bitset<> masque;
	boost::dynamic_bitset<> **gridTab = grid->getGridTab();
	unsigned long long int indice[dim - 1];

	unsigned long long int **controlIntCoords = dynsys->getControlIntCoords();
	unsigned long long int **tychIntCoords = dynsys->getTychIntCoords();

	unsigned long long int coordDiscretes[dim], imageXU[dim];

	unsigned long long int compteComm, compteTych;
	bool testNonVide, testAllTych;

	int nbIter = 0;

	vector<unsigned long long int> viabControls;
	vector<unsigned long long int> images;
	unsigned long long int posIm;

	for (int posX = 0; posX < tailleTrame; posX++)
	    {
	    //	cout<< " posX="<<posX<<endl;
	    if (!(*gridTab[posX]).none())
		{

		numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);
		for (int j = 0; j < dirTramage; j++)
		    {
		    coordDiscretes[j] = indice[j];

		    }

		for (int j = dirTramage + 1; j < dim; j++)
		    {
		    coordDiscretes[j] = indice[j - 1];

		    }

		for (unsigned long long int k = 0; k < longTrame; k++)
		    {

		    if ((*gridTab[posX])[k])
			{
			viabControls.clear();
			images.clear();
			coordDiscretes[dirTramage] = k;

			//  cout<< " point x "; printVector(coordDiscretes,dim);
			// system("pause");
			posIm = grid->intCoordsToNum(coordDiscretes);
			// cout<< "  la position du point  x "<< posIm<<endl;
			fichierB << posIm << " ";
			testNonVide = false;
			compteComm = 0;
			//	cout<<"\n analyse de point nbC= "<<nbC<<endl;
			//	cout<<" test non vide =  ";
			while (compteComm < nbC)
			    {
			    // calcul  de vecteur de controle u

			    // calcul de l'image F(x,u)

			    /* cout<< "  control coords";
			     for(int jj=0;jj<3;jj++)
			     cout<< " "<<controlIntCoords[compteComm][jj];
			     cout<<endl;*/
			    if (dynsys->constraintsXU_fd(coordDiscretes, controlIntCoords[compteComm]) < PLUS_INF)
				{
				//cout<< "  controle "<<controlIntCoords[compteComm][0]<<" "<<controlIntCoords[compteComm][1]<< " est autorise \n";

				testAllTych = true;
				for (compteTych = 0; compteTych < nbTy; compteTych++)
				    {
				    dynsys->dynamics_tych_fd(coordDiscretes, controlIntCoords[compteComm], tychIntCoords[compteTych], imageXU);
				    if (grid->isPointInGrid_fd(imageXU))
					{
					testAllTych &= grid->isInSet(imageXU);
					//cout<< " point est dans la grille de calcul ok et test no vide est "<<testNonVide<<endl;
					}
				    else
					{
					testAllTych = false;
					//cout<< "  on sort de la grille \n";
					}
				    }
				testNonVide = testAllTych;

				}
			    else
				{
				testNonVide = false;
				}
			    //	 cout<<"   "<<testNonVide;
			    if (testNonVide)
				{

				//	 	cout<< " numero image "<<posIm<<endl;
				viabControls.push_back(compteComm);
				for (compteTych = 0; compteTych < nbTy; compteTych++)
				    {
				    dynsys->dynamics_tych_fd(coordDiscretes, controlIntCoords[compteComm], tychIntCoords[compteTych], imageXU);
				    posIm = grid->intCoordsToNum(imageXU);
				    images.push_back(posIm);
				    }

				}
			    compteComm++;
			    }
			// system("pause");
			fichierB << viabControls.size() << " ";
			int cptIm = 0;
			for (unsigned long long int jj = 0; jj < viabControls.size(); jj++)
			    {
			    fichierB << viabControls.at(jj) << " ";
			    for (compteTych = 0; compteTych < nbTy; compteTych++)
				{
				fichierB << images.at(cptIm) << " ";
				cptIm++;
				}
			    }
			fichierB << endl;
			//cout<< " fini\n";

			}				// fin de if masque[k]
		    }			// fin de for  de parcours de masque
		//system("pause");

		}			//fin de if la trame n'est pas vide
	    }				// fin de for de parcours de la trame

	nbIter++;

	fichierB.close();
	// je referme le fichier

	}
    else
	// sinon
	cerr << "Erreur � l'ouverture !" << endl;

    //cout<<"fichier fini\n";

    }
