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
    //	int totalPointsC=0;
    //
    //	unsigned long long int dim=grid->dim;
    //
    //	unsigned long long int  * x=new unsigned long long int [dim];
    //
    //	double * xReel =new double[dim];
    //
    //	int totalPointsX=grid->getNbTotalPoints();
    //	unsigned long long int pos ;
    /*
     *  on parcourt  tous les points de l'espace discret  fini
     *   et  on choisit les points o� la fonction cible  renvoie une faleur finie
     */
    //for( pos=0;pos<(unsigned long long int)totalPointsX;pos++)
    //{
    /*
     * le compteur pos  ets l'unique num�ro entier du point en cours
     * dans la num�rotation alphab�tique : on parcourt axe par axe
     */

    /*
     * on restitue les  coordonn�es neti�res  du point � partir de son num�ro
     * ainis que ses coordonn�es r�elles
     */
    //		grid->numToIntAndDoubleCoords(pos,x,xReel);
    //
    //		if((*(dynsys->target))(xReel)<PLUS_INF)
    //		{
    //			totalPointsC++;
    //			grid->addPointToSet(x,1.0);
    //		}
    //	}
    initialiseTargetPointList();
    //cout<<"Total points cible est "<<totalPointsC<<"\n";
    }

void ViabiBitSet::setK0_fd()
    {
    int dirTramage = grid->getDirTram();

    unsigned long long int tailleTrame = grid->getNbPointsTotalSubGrid();

    unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();
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

    unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

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
	    xCoordsDouble[dirTramage] = limInf[dirTramage]
		    + k * gridStep[dirTramage];

	    testK = (dynsys->constraintsX(xCoordsDouble) < PLUS_INF);

	    grid->setPoint(posX, k, testK);
	    }
	}	// fin de for de parcours de la trame
    }

void ViabiBitSet::initialiseConstraints()
    {
    if ((this->dynsys->getDynType() == CC)
	    || (this->dynsys->getDynType() == DC))
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

    unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

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
	    xCoordsDouble[dirTramage] = limInf[dirTramage]
		    + k * gridStep[dirTramage];
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

    systemParams sp = *(pm->getSystemParameters());
    sp.L_LIP = 1.0;
    sp.L_MF = 1.0;
    algoViabiParams avp = *(pm->getAlgoParameters());
    controlParams cp = *(pm->getControlParameters());
    gridParams gp = *(pm->getGridParameters());

    dim = gp.DIM;
    grid = new Grid_BitSet(gp);
    dynsys = new SysDyn(sp, dim, cp, grid);

    InitViabiBitSet(avp);
    }

void ViabiBitSet::InitViabiBitSet(algoViabiParams avbp)
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

    pointDI.tabImageCells =
	    new unsigned long long int[dynsys->getTotalNbPointsC()];

    doubleVect = new double[dim];
    doubleVect1 = new double[dim];
    imageCells = new unsigned long long int[dynsys->getTotalNbPointsC()];
    //cout<< " viab bit set contructeur ok \n";
    }

void ViabiBitSet::printViabiInfo()
    {
    grid->printGrid();
    }

void ViabiBitSet::noyauViabi_FD(bool sortieOK, int nbArret)
    {
    int dirTramage = grid->getDirTram();
    unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    unsigned long long int tailleTrame = grid->getNbPointsTotalSubGrid();

    bool testF;
    unsigned long long int longTrame = grid->getLongTrame();

    unsigned long long int nbC = dynsys->getTotalNbPointsC();
    boost::dynamic_bitset<> masque;
    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(
	    longTrame, 0);
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
				if (dynsys->constraintsXU_fd(coordDiscretes,
					controlIntCoords[compteComm]) < PLUS_INF)
				    {
				    //  cout<< "  controle "<<controlIntCoords[compteComm][0]<<" "<<controlIntCoords[compteComm][1]<< " est autorise \n";
				    dynsys->dynamics_fd(coordDiscretes,
					    controlIntCoords[compteComm],
					    imageXU);
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

		    if (masquePointsEnleves->count()
			    < (unsigned long long int) longTrame)
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

int ViabiBitSet::computeViableTrajectory(double *initPosition, double finalTime,
	string fileName)
    {

    double succes = 0;
    int nbPointsCube = (int) pow(2.0, dim);	//pow(2.0, dim);
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
     * indices de déplacement pour parcourir les sommets d'une malle à partir
     * du sommet inf
     */
    long long int *indicesDecalCell = grid->getIndicesDecalCell();
    /*
     * coordonnées du points corent de la trajectoire
     */
    double xCoordsDouble[dim], imageVect[dim];
    /*
     * le pas de temps déterminé localement pour chaque point de la trajectoire
     */
    double rho;
    /*
     * numéros de mailles
     */
    int cellNum;
    /*
     * listes  pour contenir la trajectoire ( temps-position) et les contrôles
     */
    list<valarray<double> > traj, trajC;

    /*
     * structures accumulables dansune liste
     */
    valarray<double> newTrajPoint(dim + 1);
    valarray<double> trajControlCoords(dimC);

    int posTemp;

    cout << " calcul de traj a partir de coords \n";
    cout << " Postion initiale = ";

    for (int l1 = 0; l1 < dim; l1++)
	{
	cout << " " << initPosition[l1];
	}
    cout << " \n";
    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;
    int cptOK = 0;
    if (grid->isPointInGrid(initPosition))
	{
	if (dynsys->constraintsX(initPosition) < PLUS_INF)
	    {
	    cellNum = grid->localizePoint(initPosition);

	    for (int ii = 0; ii < nbPointsCube; ii++)
		{
		posTemp = cellNum + indicesDecalCell[ii];
		grid->numToIntAndDoubleCoords(posTemp, testI, testV);
		if (grid->isInSet(testI))
		    {
		    cptOK++;
		    }
		}

	    testNonVide = (cptOK > 0);

	    if (!testNonVide)
		{
		cout
			<< " La position initiale sélectionnée n'appartiant pas au noyau de viabilité\n";
		succes = 0;
		}
	    else
		{
		/*
		 * la position initiale se trouve dans le noyau de viabilité
		 * on initialise le temps à 0  et recopie la pos initiale
		 * dans le coordonnées temporaires du point en cours de la trajectoire
		 */
		double time = 0.0;
		for (int i = 0; i < dim; i++)
		    {
		    xCoordsDouble[i] = initPosition[i];
		    }
		int nbIter = 0;

		bool testviabInt = false;
		int maxnbViabPoints;
		int bestCu;
		while (time < finalTime && nbIter <= NB_MAX_TRAJ_ITER)
		    {
		    for (int i = 0; i < dim; i++)
			{
			newTrajPoint[i] = xCoordsDouble[i];
			}
		    newTrajPoint[dim] = time;
		    traj.push_back(newTrajPoint);

		    rho = dynsys->calculRho_local(xCoordsDouble);

		    rho = min(rho, finalTime - time);

		    bestCu = this->findViabControl(xCoordsDouble, rho, 1, 1.0,
			    imageVect, maxnbViabPoints, testNonVide);
		    time += rho;
		    testviabInt = (maxnbViabPoints == nbPointsCube);

		    // la boucle s'arête ici u premier contrôle
		    // qui donne un successeur viable

		    if (testNonVide)
			{
			// contrôle viable trouvé
			// on recopie ce contrôle dans la liste et
			// le successeur devient le point  courent
			if (testviabInt)
			    {

			    for (int dc = 0; dc < dimC; dc++)
				{
				trajControlCoords[dc] =
					controlCoords[bestCu][dc];
				}
			    trajC.push_back(trajControlCoords);
			    for (int i = 0; i < dim; i++)
				{
				xCoordsDouble[i] = imageVect[i];
				}
			    }
			else
			    {

			    rho = dynsys->calculRho_local(xCoordsDouble);

			    rho = min(rho, finalTime - time);

			    bestCu = this->findViabControl(xCoordsDouble, rho,
				    10, 0.95, imageVect, maxnbViabPoints,
				    testNonVide);

			    for (int dc = 0; dc < dimC; dc++)
				{
				trajControlCoords[dc] =
					controlCoords[bestCu][dc];
				}
			    trajC.push_back(trajControlCoords);
			    for (int i = 0; i < dim; i++)
				{
				xCoordsDouble[i] = imageVect[i];
				}
			    }
			}
		    else
			{
			cout << "   Echec! Sortie de l'ensemble viable \n";
			break;

			}
		    nbIter++;
		    }
		if (time >= finalTime)
		    {
		    succes = 1.0;
		    }
		}		//fin de else (reconstruction de trajectoire)

	    }
	else
	    {
	    printf(" Point initial hors de l'ensemble de contraintes. Arret\n");
	    succes = 0;
	    }
	}
    else
	{
	printf(" Point initial hors de grille de calcul. Arret\n");
	succes = 0;
	}

    // if(succes)
	{

	printf(" trajectoire trouvée. Enregistrement\n");

	FILE *fi;
	fi = fopen(fileName.c_str(), "w");
	if (fi == NULL)
	    {
	    printf("** error: impossible to open the file %s.\n",
		    fileName.c_str());

	    }
	else
	    {
	    list<valarray<double> >::iterator it = traj.begin();
	    list<valarray<double> >::iterator itc = trajC.end();
	    itc--;
	    trajC.push_back(*itc);
	    itc = trajC.begin();

	    while (it != traj.end())
		{

		for (int l1 = 0; l1 < dim; l1++)
		    {
		    fprintf(fi, "%15.8f ", (*it)[l1]);
		    cout << " " << (*it)[l1];
		    }
		fprintf(fi, "%15.8f ", (*it)[dim]);
		cout << " " << (*it)[dim];
		for (int dc = 0; dc < dimC; dc++)
		    {
		    fprintf(fi, "%15.8f ", (*itc)[dc]);
		    cout << " " << (*itc)[dc] << endl;
		    }
		fprintf(fi, "\n");
		it++;
		itc++;
		//   traj.pop_front();
		//   trajC.pop_front();
		}
	    fclose(fi);
	    }
	}
    return succes;

    }

int ViabiBitSet::computeViableTrajectoryHeavy(double *initPosition,
	double *initControl, double finalTime, string fileName)
    {

    double succes = 0;
    int nbPointsCube = (int) pow(2.0, dim);			//pow(2.0, dim);
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
     * indices de déplacement pour parcourir les sommets d'une malle à partir
     * du sommet inf
     */
    long long int *indicesDecalCell = grid->getIndicesDecalCell();
    /*
     * coordonnées du points corent de la trajectoire
     */
    double xCoordsDouble[dim], imageVect[dim], currentControl[dimC];
    /*
     * le pas de temps déterminé localement pour chaque point de la trajectoire
     */
    double rho;
    /*
     * numéros de mailles
     */
    int cellNum;
    /*
     * listes  pour contenir la trajectoire ( temps-position) et les contrôles
     */
    list<valarray<double> > traj, trajC;

    /*
     * structures accumulables dansune liste
     */
    valarray<double> newTrajPoint(dim + 1);
    valarray<double> trajControlCoords(dimC);

    int posTemp;

    cout << " calcul de traj a partir de coords \n";
    cout << " Postion initiale = ";

    for (int l1 = 0; l1 < dim; l1++)
	{
	cout << " " << initPosition[l1];
	}
    cout << " \n";
    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;
    int cptOK = 0;
    if (grid->isPointInGrid(initPosition))
	{
	if (dynsys->constraintsX(initPosition) < PLUS_INF)
	    {
	    cellNum = grid->localizePoint(initPosition);

	    for (int ii = 0; ii < nbPointsCube; ii++)
		{
		posTemp = cellNum + indicesDecalCell[ii];
		grid->numToIntAndDoubleCoords(posTemp, testI, testV);
		if (grid->isInSet(testI))
		    {
		    cptOK++;
		    }
		}

	    testNonVide = (cptOK > 0);

	    if (!testNonVide)
		{
		cout
			<< " La position initiale sélectionnée n'appartiant pas au noyau de viabilité\n";
		succes = 0;
		}
	    else
		{
		/*
		 * la position initiale se trouve dans le noyau de viabilité
		 * on initialise le temps à 0  et recopie la pos initiale
		 * dans le coordonnées temporaires du point en cours de la trajectoire
		 */
		double time = 0.0;
		for (int i = 0; i < dim; i++)
		    {
		    xCoordsDouble[i] = initPosition[i];
		    }

		for (int i = 0; i < dimC; i++)
		    {
		    currentControl[i] = initControl[i];
		    }

		int nbIter = 0;
		/*
		 * On itère tant que le temps n'a pas dépassé l'horizon donné
		 */

		bool testviabInt = false;
		int maxnbViabPoints;
		int bestCu;
		while (time < finalTime && nbIter <= NB_MAX_TRAJ_ITER)
		    {
		    // cout<< " point en cours ";
		    for (int i = 0; i < dim; i++)
			{
			newTrajPoint[i] = xCoordsDouble[i];
			// cout<< " "<<newTrajPoint[i];
			}
		    newTrajPoint[dim] = time;
		    //  cout<< " temps= "<<newTrajPoint[dim]<<endl;
		    traj.push_back(newTrajPoint);

		    rho = dynsys->calculRho_local(xCoordsDouble);

		    rho = min(rho, finalTime - time);
		    dynsys->setRho(rho);
		    // cout<< " rho= "<<rho<<endl;
		    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble,
			    currentControl, doubleVect1, rho);
		    testNonVide = false;
		    if (grid->isPointInGrid(doubleVect1))
			{
			/*
			 *  le sucesseur est dans la grille de calcul
			 */
			if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
			    {
			    cellNum = grid->localizePoint(doubleVect1);
			    // cout<< " num cellule "<<cellNum<<endl;

			    cptOK = 0;
			    for (int ii = 0; ii < nbPointsCube; ii++)
				{
				posTemp = cellNum + indicesDecalCell[ii];
				grid->numToIntAndDoubleCoords(posTemp, testI,
					testV);
				if (grid->isInSet(testI))
				    {
				    cptOK++;
				    }
				}
			    testviabInt = (cptOK == nbPointsCube);
			    testNonVide = (cptOK > 0);
			    }
			}
		    if (testNonVide)
			{
			// contrôle viable trouvé
			// on recopie ce contrôle dans la liste et
			// le successeur devient le point  courent
			if (testviabInt)
			    {

			    //   cout<<  " image interieure tourvee \n";

			    for (int dc = 0; dc < dimC; dc++)
				{
				trajControlCoords[dc] = currentControl[dc];
				}
			    trajC.push_back(trajControlCoords);
			    for (int i = 0; i < dim; i++)
				{
				xCoordsDouble[i] = doubleVect1[i];
				}
			    }
			else
			    {
			    cout
				    << " ======================= Recalage pendant la phase controle constant  =======================\n";
			    grid->findNearestViabPointInCell(xCoordsDouble,
				    doubleVect1, xCoordsDouble,
				    dynsys->dynConstraintsForTraj);
			    for (int dc = 0; dc < dimC; dc++)
				{
				trajControlCoords[dc] = currentControl[dc];
				}
			    trajC.push_back(trajControlCoords);

			    }
			}
		    else
			{
			cout
				<< " ======================= CHANGEMENT DE CONTROLE =======================\n";

			bestCu = this->findViabControl(xCoordsDouble, rho, 1,
				1.0, imageVect, maxnbViabPoints, testNonVide);

			testviabInt = (maxnbViabPoints == nbPointsCube);

			// la boucle s'arête ici u premier contrôle
			// qui donne un successeur viable

			//  cout<<   " Premiere recherche de controle viable  fini parcours de controles on a test interieur = "<<testviabInt<<
			//      " test non vide "<<testNonVide<< " maxnbViabPoints =  "<<maxnbViabPoints<< " bes c u= "<<bestCu<<endl;
			if (testNonVide)
			    {
			    // contrôle viable trouvé
			    // on recopie ce contrôle dans la liste et
			    // le successeur devient le point  courent
			    if (testviabInt)
				{

				//   cout<<  " image interieure tourvee \n";

				for (int dc = 0; dc < dimC; dc++)
				    {
				    trajControlCoords[dc] =
					    controlCoords[bestCu][dc];
				    currentControl[dc] =
					    controlCoords[bestCu][dc];
				    }
				trajC.push_back(trajControlCoords);
				for (int i = 0; i < dim; i++)
				    {
				    xCoordsDouble[i] = imageVect[i];
				    }
				}
			    else
				{
				cout
					<< " ======================= Recalage =======================\n";
				grid->findNearestViabPointInCell(xCoordsDouble,
					imageVect, xCoordsDouble,
					dynsys->dynConstraintsForTraj);

				cout << " controle ";
				for (int dc = 0; dc < dimC; dc++)
				    {
				    trajControlCoords[dc] =
					    controlCoords[bestCu][dc];
				    currentControl[dc] =
					    controlCoords[bestCu][dc];
				    cout << " " << currentControl[dc];
				    }
				trajC.push_back(trajControlCoords);

				}
			    }
			else
			    {
			    cout << "   Echec! Sortie de l'ensemble viable \n";
			    break;

			    }
			}
		    time += rho;
		    nbIter++;
		    }
		if (time >= finalTime)
		    {
		    succes = 1.0;
		    }
		}		//fin de else (reconstruction de trajectoire)

	    }
	else
	    {
	    printf(" Point initial hors de l'ensemble de contraintes. Arret\n");
	    succes = 0;
	    }
	}
    else
	{
	printf(" Point initial hors de grille de calcul. Arret\n");
	succes = 0;
	}

    // if(succes)
	{

	printf(" trajectoire trouvée. Enregistrement\n");

	FILE *fi;
	fi = fopen(fileName.c_str(), "w");
	if (fi == NULL)
	    {
	    printf("** error: impossible to open the file %s.\n",
		    fileName.c_str());

	    }
	else
	    {
	    list<valarray<double> >::iterator it = traj.begin();
	    list<valarray<double> >::iterator itc = trajC.end();
	    itc--;
	    trajC.push_back(*itc);
	    itc = trajC.begin();

	    while (it != traj.end())
		{

		for (int l1 = 0; l1 < dim; l1++)
		    {
		    fprintf(fi, "%15.8f ", (*it)[l1]);
		    //  cout<< " "<<  (*it)[l1];
		    }
		fprintf(fi, "%15.8f ", (*it)[dim]);
		//cout<< " "<<(*it)[dim];
		for (int dc = 0; dc < dimC; dc++)
		    {
		    fprintf(fi, "%15.8f ", (*itc)[dc]);
		    // cout<< " "<<(*itc)[dc]<<endl;
		    }
		fprintf(fi, "\n");
		it++;
		itc++;
		//   traj.pop_front();
		//   trajC.pop_front();
		}
	    fclose(fi);
	    }
	}
    return succes;

    }

int ViabiBitSet::findViabControl(double *currentPos, double &dt, int nbStepIter,
	double stepCoeff, double *resPos, int &nbViabVoisins, bool &succes)
    {

    int nbPointsCube = (int) pow(2.0, dim);			//pow(2.0, dim);
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

    /*
     * indices de déplacement pour parcourir les sommets d'une malle à partir
     * du sommet inf
     */
    long long int *indicesDecalCell = grid->getIndicesDecalCell();
    /*
     * coordonnées du points corent de la trajectoire
     */
    double xCoordsDouble[dim];
    /*
     * le pas de temps déterminé localement pour chaque point de la trajectoire
     */
    double rho;
    /*
     * numéros de mailles
     */
    int cellNum;

    int posTemp;

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;
    int cptOK = 0;

    for (int i = 0; i < dim; i++)
	{
	xCoordsDouble[i] = currentPos[i];
	}
    bool testviabInt = false;
    int maxnbViabPoints;
    int bestCu;
    rho = dt;
    dynsys->setRho(rho);
    // cout<< " rho= "<<rho<<endl;

    /*
     * on parcours tous les contrôles
     */
    maxnbViabPoints = 0;
    bestCu = 0;

    int iter = 0;

    testviabInt = false;
    testNonVide = false;
    double hMax = grid->maxStep;
    while (iter < nbStepIter && !testviabInt)
	{
	unsigned long long int cu = 0;
	testviabInt = false;
	testNonVide = false;
	dt = rho;
	while ((cu < nbCTotal) && !testviabInt)
	    {
	    /*
	     * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	     * au point en cours
	     */
	    if (dynsys->constraintsXU(xCoordsDouble,
		    controlCoords[cu])<PLUS_INF)
		{
		/*
		 * calcul du successeur  du point en cours
		 */
		rho = dynsys->calculRho_local(xCoordsDouble);
		(dynsys->*(dynsys->discretDynamics))(xCoordsDouble,
			controlCoords[cu], doubleVect1, rho);
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
			cptOK = 0;
			for (int ii = 0; ii < nbPointsCube; ii++)
			    {
			    posTemp = cellNum + indicesDecalCell[ii];
			    grid->numToIntAndDoubleCoords(posTemp, testI,
				    testV);
			    if (dynsys->dynConstraintsForTraj(xCoordsDouble,
				    testV) < PLUS_INF)// on v�rifie si la projection corrspond aux contraintes de dynamique
				{
				double dist = 0.0;
				for (int k = 0; k < dim; k++)
				    {
				    dist = max(dist,
					    abs(testV[k] - doubleVect1[k]));
				    }
				testviabInt = (grid->isInSet(testI)
					&& (dist <= hMax / 2.0));
				if (grid->isInSet(testI))
				    {
				    cptOK++;
				    }
				}

			    }
			//testviabInt=(cptOK==nbPointsCube);
			testNonVide = (cptOK > 0);
			if (cptOK >= maxnbViabPoints)
			    {
			    maxnbViabPoints = cptOK;
			    bestCu = cu;
			    }
			}
		    }
		}
	    cu++;
	    }				//fin de parcours de tous les contrôles
					// la boucle s'arête ici u premier contrôle
					// qui donne un successeur viable
	iter++;
	rho = rho * stepCoeff;
	dynsys->setRho(rho);
	//   cout<< " iteration  "<<iter<<" rho= "<<rho<<endl;
	}

    succes = testNonVide;
    nbViabVoisins = maxnbViabPoints;
    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[bestCu],
	    doubleVect1, rho);
    for (int i = 0; i < dim; i++)
	{
	resPos[i] = doubleVect1[i];
	}
    return bestCu;
    }

bool ViabiBitSet::findViabImagePoint_noControl(double *xCoordsDouble,
	bool print)
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
    long long int *indicesDecalCell = grid->getIndicesDecalCell();
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
	    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble,
		    controlCoords[cu], doubleVect1, rho);
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
			if (dynsys->dynConstraintsForTraj(xCoordsDouble,
				testV)<PLUS_INF)
			    {
			    double dist = 0.0;
			    for (int k = 0; k < dim; k++)
				{
				dist = max(dist,
					abs(testV[k] - doubleVect1[k]));
				}
			    testNonVide = grid->isInSet(testI)
				    && (dist <= hMax / 2.0);

			    }
			ii++;
			}
		    }
		}
	    else
		{
		//cout<< " sortie du domaine : autorisation : "<<grid->unboundedDomain<< "point "<< doubleVect1[0]<<" "<<doubleVect1[1]<<endl;
		testNonVide = grid->unboundedDomain
			&& grid->isPointInGridWithConstr(doubleVect1)
			&& (dynsys->constraintsX(doubleVect1) < PLUS_INF);
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
    long long int *indicesDecalCell = grid->getIndicesDecalCell();
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
	    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble,
		    controlCoords[cu], doubleVect1, rho);
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
			if (dynsys->dynConstraintsForTraj(xCoordsDouble,
				testV)<PLUS_INF)
			    {
			    double dist = 0.0;
			    for (int k = 0; k < dim; k++)
				{
				dist = max(dist,
					abs(testV[k] - doubleVect1[k]));
				}
			    testNonVide = grid->isInSet(testI)
				    && (dist <= hMax / 2.0);
			    ii++;
			    }
			}
		    }
		}
	    else
		{
		testNonVide = grid->unboundedDomain
			&& grid->isPointInGridWithConstr(doubleVect1)
			&& (dynsys->constraintsX(doubleVect1) < PLUS_INF);
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
		if (dynsys->constraintsXU(xCoordsDouble,
			controlCoords[cu])<PLUS_INF)
		    {
		    /*
		     * calcul du successeur  du point en cours
		     */
		    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble,
			    controlCoords[cu], doubleVect1, rho);
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
				grid->numToIntAndDoubleCoords(posTemp, testI,
					testV);
				if (dynsys->dynConstraintsForTraj(xCoordsDouble,
					testV) < PLUS_INF)
				    {
				    double dist = 0.0;
				    for (int k = 0; k < dim; k++)
					{
					dist = max(dist,
						abs(testV[k] - doubleVect1[k]));
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

int ViabiBitSet::computeViableTrajectorySetVal(double *initPosition,
	double finalTime, string fileName)
    {

    double succes = 0;
    int nbPointsCube = (int) pow(2.0, dim);			//pow(2.0, dim);
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
     * indices de déplacement pour parcourir les sommets d'une malle à partir
     * du sommet inf
     */
    long long int *indicesDecalCell = grid->getIndicesDecalCell();
    /*
     * coordonnées du points corent de la trajectoire
     */
    double xCoordsDouble[dim], imageVect[dim];
    /*
     * le pas de temps déterminé localement pour chaque point de la trajectoire
     */
    double rho;
    /*
     * numéros de mailles
     */
    int cellNum;
    /*
     * listes  pour contenir la trajectoire ( temps-position) et les contrôles
     */
    list<valarray<double> > traj, trajC;

    /*
     * structures accumulables dansune liste
     */
    valarray<double> newTrajPoint(dim + 1);
    valarray<double> trajControlCoords(dimC);

    int posTemp;

    cout << " calcul de traj a partir de coords \n";
    cout << " Postion initiale = ";

    for (int l1 = 0; l1 < dim; l1++)
	{
	cout << " " << initPosition[l1];
	}
    cout << " \n";
    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;
    int cptOK = 0;
    if (grid->isPointInGrid(initPosition))
	{
	if (dynsys->constraintsX(initPosition) < PLUS_INF)
	    {
	    cellNum = grid->localizePoint(initPosition);

	    for (int ii = 0; ii < nbPointsCube; ii++)
		{
		posTemp = cellNum + indicesDecalCell[ii];
		grid->numToIntAndDoubleCoords(posTemp, testI, testV);
		if (grid->isInSet(testI))
		    {
		    cptOK++;
		    }
		}

	    testNonVide = (cptOK > 0);

	    if (!testNonVide)
		{
		cout
			<< " La position initiale sélectionnée n'appartiant pas au noyau de viabilité\n";
		succes = 0;
		}
	    else
		{
		/*
		 * la position initiale se trouve dans le noyau de viabilité
		 * on initialise le temps à 0  et recopie la pos initiale
		 * dans le coordonnées temporaires du point en cours de la trajectoire
		 */
		double time = 0.0;
		for (int i = 0; i < dim; i++)
		    {
		    xCoordsDouble[i] = initPosition[i];
		    }
		int nbIter = 0;

		bool testviabInt = false;
		int maxnbViabPoints;
		int bestCu;
		while (time < finalTime && nbIter <= NB_MAX_TRAJ_ITER)
		    {
		    // cout<< " point en cours ";
		    for (int i = 0; i < dim; i++)
			{
			newTrajPoint[i] = xCoordsDouble[i];
			//  cout<< " "<<newTrajPoint[i];
			}
		    newTrajPoint[dim] = time;
		    // cout<< " temps= "<<newTrajPoint[dim]<<endl;
		    traj.push_back(newTrajPoint);

		    rho = dynsys->calculRho_local(xCoordsDouble);

		    rho = min(rho, finalTime - time);

		    //  cout<< " rho= "<<rho<<endl;

		    bestCu = this->findViabControl(xCoordsDouble, rho, 1, 1.0,
			    imageVect, maxnbViabPoints, testNonVide);
		    time += rho;
		    testviabInt = (maxnbViabPoints == nbPointsCube);

		    // la boucle s'arête ici u premier contrôle
		    // qui donne un successeur viable

		    //    cout<<   " Premiere recherche de controle viable  fini parcours de controles on a test interieur = "<<testviabInt<<
		    //        " test non vide "<<testNonVide<< " maxnbViabPoints =  "<<maxnbViabPoints<< " bes c u= "<<bestCu<<endl;
		    if (testNonVide)
			{
			// contrôle viable trouvé
			// on recopie ce contrôle dans la liste et
			// le successeur devient le point  courent
			if (testviabInt)
			    {

			    //   cout<<  " image interieure tourvee \n";

			    for (int dc = 0; dc < dimC; dc++)
				{
				trajControlCoords[dc] =
					controlCoords[bestCu][dc];
				}
			    trajC.push_back(trajControlCoords);
			    for (int i = 0; i < dim; i++)
				{
				xCoordsDouble[i] = imageVect[i];
				}
			    }
			else
			    {
			    //     cout<< " ======================= Recalage =======================\n";
			    grid->findNearestViabPointInCell(xCoordsDouble,
				    imageVect, xCoordsDouble,
				    dynsys->dynConstraintsForTraj);

			    for (int dc = 0; dc < dimC; dc++)
				{
				trajControlCoords[dc] =
					controlCoords[bestCu][dc];
				}
			    trajC.push_back(trajControlCoords);

			    }
			}
		    else
			{
			cout << "   Echec! Sortie de l'ensemble viable \n";
			break;

			}
		    nbIter++;
		    }
		if (time >= finalTime)
		    {
		    succes = 1.0;
		    }
		}		//fin de else (reconstruction de trajectoire)

	    }
	else
	    {
	    printf(" Point initial hors de l'ensemble de contraintes. Arret\n");
	    succes = 0;
	    }
	}
    else
	{
	printf(" Point initial hors de grille de calcul. Arret\n");
	succes = 0;
	}

    // if(succes)
	{

	printf(" trajectoire trouvée. Enregistrement\n");

	FILE *fi;
	fi = fopen(fileName.c_str(), "w");
	if (fi == NULL)
	    {
	    printf("** error: impossible to open the file %s.\n",
		    fileName.c_str());

	    }
	else
	    {
	    list<valarray<double> >::iterator it = traj.begin();
	    list<valarray<double> >::iterator itc = trajC.end();
	    itc--;
	    trajC.push_back(*itc);
	    itc = trajC.begin();

	    while (it != traj.end())
		{

		for (int l1 = 0; l1 < dim; l1++)
		    {
		    fprintf(fi, "%15.8f ", (*it)[l1]);
		    //  cout<< " "<<  (*it)[l1];
		    }
		fprintf(fi, "%15.8f ", (*it)[dim]);
		//cout<< " "<<(*it)[dim];
		for (int dc = 0; dc < dimC; dc++)
		    {
		    fprintf(fi, "%15.8f ", (*itc)[dc]);
		    // cout<< " "<<(*itc)[dc]<<endl;
		    }
		fprintf(fi, "\n");
		it++;
		itc++;
		//   traj.pop_front();
		//   trajC.pop_front();
		}
	    fclose(fi);
	    }
	}
    return succes;

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

    currentImagePointsList.maxNum = 0;
    currentImagePointsList.minNum = 0;
    currentImagePointsList.pointsList = list<unsigned long long int>();

    unsigned long long int pos;

    list<unsigned long long int>::iterator itStart =
	    currentImagePointsList.pointsList.begin(), itNew;

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
		addDataToPointsList(&itStart, pos, &itNew);
		itStart = itNew;

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

    cout << "  Compute First confexified  Image . Iteration num " << iter
	    << " points dans Cn = " << currentImagePointsList.pointsList.size()
	    << endl;

    unsigned long long int posX;

    list<unsigned long long int>::iterator itStart, itNew;

    currentImageList.cellsList.clear();
    currentImageList.maxNum = 0;
    currentImageList.minNum = PLUS_INF;
    unsigned long long int intPointCoords[dim];
    double doublePointCoords[dim];

    double rho;

    list<unsigned long long int>::iterator itPoint =
	    currentImagePointsList.pointsList.begin(), itLastPoint =
	    currentImagePointsList.pointsList.end(), itTemp;
    unsigned long long int tempImageCell;

    while (itPoint != itLastPoint)
	{
	posX = (*itPoint);
	//cout<< " posX="<<posX<<endl;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	//printVector(doublePointCoords, dim);
	////////system("pause");
	rho = dynsys->calculRho_local(doublePointCoords);
	//cout<< " rho= "<<rho;

	this->computeDiscreteImageOfPoint(doublePointCoords, intPointCoords);
	/*!
	 * L'image calcul�e est stock�e dans la structure pointDI sous forme de tableau ordonn�
	 *  de cellules
	 */
	unsigned long long int numCell;

	itStart = this->currentImageList.cellsList.begin();

	//cout<<  " analyse d'une image de point\n";
	for (unsigned long long int i = 0; i < pointDI.nbImageCells; i++)
	    {
	    numCell = pointDI.tabImageCells[i];
	    tempImageCell = numCell;

	    //cout<< " ajout  de cellule \n ";
	    addDataToCurrentImage(&itStart, tempImageCell, &itNew);

	    itStart = itNew;
	    addConvexCombinations(posX, numCell, &tempImageCell, rho, &itStart);
	    //cout<< " ares ajout  de l'image d'un point la taille de la liste est "<<currentImageList.cellsList.size()<<endl;

	    }

	itPoint++;

	}

    cout << " parcours de base termin�\n";

    }

void ViabiBitSet::addConvexCombinations(unsigned long long int posX,
	unsigned long long int numCell, unsigned long long int *tempImageCell,
	double rho, list<unsigned long long int>::iterator *itStart)
    {

    unsigned long long int intPointCoords[dim];
    double doublePointCoords[dim];

    //cout<< "add convex comnbin  posX= "<<posX<< " rho ="<<rho<<endl;
    grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
    //cout<< " point de d�part ";
    //    printVector(doublePointCoords, dim);
    //    cout<< " numCell = "<<numCell<<endl;
    grid->numToIntAndDoubleCoords(numCell, intPointCoords, doubleVect);
    //cout<< " point arrivee  ";
    //printVector(doubleVect, dim);
    list<unsigned long long int>::iterator itNew;
    double dist = 0.;
    /*!
     * Ensuite on calcule le vecteur diff�rence \f$ z=y-x\f$ et sa norme \f$ \|z\|_2 \f$ .
     */
    for (int i = 0; i < dim; i++)
	{
	doubleVect[i] = doubleVect[i] - doublePointCoords[i];
	dist += doubleVect[i] * doubleVect[i];
	}
    /*!
     * Le segment reliant \f$ x\f$ et \f$ y\simeq x-\rho F(x,u)\f$ pour un certain \f$ u\in U(x)\f$  a pour �quation
     * \f[
     * x+tz,\ t\in[0,1]
     * \f]
     * On d�termine alors un pas de progression \f[
     * \Delta t=\min(0.1, \frac{h_{max}}{\|z\|_2})
     * \f]
     * de fa�on � pouvoir passer d'une maille � l'autre en avan�ant avec  ce pas le long du segment
     * On parcourt ensuite le segment avec le pas calcul� et on ajoute � l'image les mailles crois�es
     * avec pour valeur une partie du pas de temps \f$\rho\f$ .
     */
    double deltat = min(0.1, grid->getMaxStep() / (2.0 * sqrt(dist)));
    //    cout<< " deltat= "<<deltat<<endl;
    //    cout<< "   norme de difference "<<sqrt(dist)<< " pas maxi "<<grid->getMaxStep()<<endl<< " difference vecteur ";
    //printVector(doubleVect, dim);
    double t = deltat;
    unsigned long long int newCellNum, lastVisitCellNum =
	    grid->getNbTotalCells() + 1;
    while (t < 1.0)
	{

	for (int i = 0; i < dim; i++)
	    {
	    doubleVect1[i] = doublePointCoords[i] + t * doubleVect[i];
	    }

	// cout<< "  point intermediaire  ";
	// printVector(doubleVect1, dim);
	if (grid->isPointInGrid(doubleVect1))
	    {
	    //cout<< " le point est das la grlle\n";
	    /*!
	     * Si l'image est  dans les limites de la grille on �tudie si elle v�rifie les contraintes
	     */

	    if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
		{
		// cout<< "  contraintes sur   X  ok\n";

		/*!
		 * Si l'image  est dans l'ensemble de contraintes sur l'�tat \f$ K \f$
		 * on calcule le num�ro de maille qui contient cette image
		 */
		newCellNum = grid->localizePoint(doubleVect1); // on enregistre le numero de maille
		if (newCellNum != lastVisitCellNum)
		    {
		    //cout<< " maille visee est "<<newCellNum<<endl;
		    (*tempImageCell) = newCellNum;
		    //itStart=this->currentImageList.cellsList.begin();
		    // cout<<  "  ajout de celule "<<newCellNum<< " avec val = "<<t*rho*(dynsys->lFunc(doubleVect1,cVect))<<endl;
		    this->addDataToCurrentImage(itStart, (*tempImageCell),
			    &itNew);
		    lastVisitCellNum = newCellNum;
		    (*itStart) = itNew;
		    //this->showCurrentImageList();
		    ////////system("pause");
		    }
		}
	    }
	t = t + deltat;

	}
    ////////system("pause");

    }
void ViabiBitSet::addDataToCurrentImage(
	list<unsigned long long int>::iterator *startIt,
	unsigned long long int newCell,
	list<unsigned long long int>::iterator *resIt)
    {
    /*!
     * Dans l'image  en construction de \f$ C_{n}\setminus C_{n-1}\f$ les mailles sont stock�es dans l'ordre croissant
     * de leur num�ros  et chaque maille garde la m�moire de la valeur minimale ainsi que de tous les ant�c�dents
     * de cette maille c'est � dire tous les couples viables (x,u) pour lesquels f(x,u) appartient � cette maille.
     * En plus  on distingue les couples (x,u) minimisant la valeur  et les autres, juste viables.
     *
     * Pour ajouter une nouvelle maille dans cette liste  on doit faire dans l'ordre les choses suivantes :
     *
     *    rechercher dans la liste la place nu num�ro de maille � inserrer
     *   deux cas de figure peuvent se pr�senter :
     *     la maille ayant le m�me num�ro  existe d�j�: on proc�de alors � la fusion des deux,  en d�terminant la valeur optimale et
     *    en  fusionnant les r�tro-actions
     *     la maille n'existe pas encore; alors on l'inserre simplement  dans l'ordre croissant.
     *
     */

    unsigned long long int numNewCell = newCell;

    list<unsigned long long int>::iterator itCell, itLast =
	    currentImageList.cellsList.end();

    if (currentImageList.cellsList.size() == 0)
	{
	currentImageList.cellsList.push_back(newCell);
	currentImageList.maxNum = numNewCell;
	currentImageList.minNum = numNewCell;
	(*resIt) = currentImageList.cellsList.end();
	(*resIt)--;
	}
    else
	{

	if (numNewCell > currentImageList.maxNum)
	    {
	    currentImageList.cellsList.push_back(newCell);
	    currentImageList.maxNum = numNewCell;
	    currentImageList.minNum = min(currentImageList.maxNum,
		    currentImageList.minNum);
	    (*resIt) = currentImageList.cellsList.end();
	    (*resIt)--;
	    }
	else
	    {
	    if (numNewCell < currentImageList.minNum)
		{
		currentImageList.cellsList.push_front(newCell);
		currentImageList.minNum = numNewCell;
		currentImageList.maxNum = max(currentImageList.maxNum,
			currentImageList.minNum);
		(*resIt) = currentImageList.cellsList.begin();
		}
	    else
		{
		itCell = *startIt;

		if ((numNewCell > (*itCell)))
		    {
		    while ((itCell != itLast) && (numNewCell > (*itCell)))
			{
			//cout<< " current cell num "<<(*itCell).cellNum<< " new cell num "<<numNewCell<<endl;
			itCell++;
			}
		    if (numNewCell < (*itCell))
			{
			currentImageList.cellsList.insert(itCell, newCell);
			}

		    (*resIt) = itCell;
		    }
		else
		    {
		    while ((numNewCell < (*itCell)))
			{
			//cout<< " current cell num "<<(*itCell).cellNum<< " new cell num "<<numNewCell<<endl;
			itCell--;
			}
		    if (numNewCell > (*itCell))
			{
			itCell++;
			currentImageList.cellsList.insert(itCell, newCell);
			}

		    (*resIt) = itCell;
		    }

		}
	    }

	}
    }

void ViabiBitSet::computeDiscreteImageOfPoint(double *doublePointCoords,
	unsigned long long int *intPointCoords)
    {

    //double testV[dim];
    //unsigned long long int testI[dim];

    //cout<< " calcul de l'image d'un point \n";
    double **controlCoords = dynsys->getControlCoords();
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int nbCellsTotal = grid->getNbTotalPoints();
    unsigned long long int cu;

    list<unsigned long long int> cellsList;
    double rho;

    //cout<< " calcul de l'image pour le point ";
    // printVector(doublePointCoords,dim);
    //double dv[dim];
    for (cu = 0; cu < nbCTotal; cu++)
	{

	// cout<< "  controle numero "<<cu<<endl;
	// printVector(controlCoords[cu],dimC);
	/*!
	 * on calcule les coordonn�es r�elles de l'image du point par la dynamique discrete
	 * elles sont stockes dans le tableau doubleVect
	 * si \f$ u\in U(x)\f$  (contr�le admissible) on calcule l'image r�elle par la dynamique
	 *  discr�tis�e  en temps  du point x avec le controle u
	 */
	if (dynsys->constraintsXU(doublePointCoords,
		controlCoords[cu])<PLUS_INF)
	    {
	    //////printf("  contraintes sur U et X  ok\n");
	    //  ////printf( " x= ");
	    // printVector(doublePointCoords,dim);
	    rho = dynsys->calculRho_local(doublePointCoords);

	    (dynsys->*(dynsys->discretDynamics))(doublePointCoords,
		    controlCoords[cu], doubleVect1, rho);

	    // cout<< " retrour dnamique discrete ";
	    // printVector(doubleVect1, dim);
	    //(dynsys->*discretDynamics)(doublePointCoords, controlCoords[cu], doubleVect);
	    if (grid->isPointInGrid(doubleVect1))
		{
		//  printf( "   le point est das la grlle\n " );

		//////printf(" le point est das la grlle\n");
		/*!
		 * Si l'image est  dans les limites de la grille on �tudie si elle v�rifie les contraintes
		 */
		if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
		    {
		    //     printf( "   contraintes Xu ok\n " );

		    //     ////printf("  contraintes sur   X  ok\n");

		    /*!
		     * Si l'image  est dans l'ensemble de contraintes sur l'�tat \f$ K \f$
		     * on calcule le num�ro de maille qui contient cett image
		     */
		    imageCells[cu] = grid->localizePoint(doubleVect1);

		    // on enregistre le numero de maille
		    cellsList.push_back(imageCells[cu]);
		    //////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,unsigned long long ints[cu]);
		    //   cout<< " les coordonn�es  de cin bas de la maille  sont ";
		    //   grid->numToIntAndDoubleCoords(unsigned long long ints[cu], intPointCoords, doubleVect1);
		    //   printVector(doubleVect1, dim);
		    }
		else
		    {
		    /*!
		     * Si l'image n'est pas dans  \f$ K \f$ on enregistre un num�ro de maille factice qui signifie que
		     * cette image est rejet�e
		     */
		    imageCells[cu] = nbCellsTotal + 10; // sinon on enregistre un nombre convenu reconnaissanble
		    cellsList.push_back(imageCells[cu]);
		    //    ////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,unsigned long long ints[cu]);
		    }
		}
	    else
		{
		/*!
		 * Si l'image n'est pas dans  la grille on enregistre un num�ro de maille factice qui signifie que
		 * cette image est rejet�e
		 */
		if (grid->isPointInGridWithConstr(doubleVect1))
		    {
		    imageCells[cu] = nbCellsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble
		    cellsList.push_back(imageCells[cu]);

		    }
		else
		    {
		    imageCells[cu] = nbCellsTotal + 10; // sinon on enregistre un nombre convenu reconnaissanble
		    cellsList.push_back(imageCells[cu]);
		    }
		//////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,unsigned long long ints[cu]);
		}

	    //cout<< " ajour d'une paire � la liste longueur "<<cellsList.size()<<"\n";

	    }
	else
	    {
	    /*!
	     * Si l'image n'est pas dans  la grille on enregistre un num�ro de maille factice qui signifie que
	     * cette image est rejet�e
	     */
	    /*!
	     * \todo pr�voir la gestion des autorisations de sortie  par axe pour les variables qui peuvent
	     * se trouver dans des intervalles non born�s
	     */
	    imageCells[cu] = nbCellsTotal + 10; // sinon on enregistre un nombre convenu reconnaissanble
	    cellsList.push_back(imageCells[cu]);
	    //////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,unsigned long long ints[cu]);
	    }

	}
    /*    cout<< " fini traitement parallele\n";
     for( cu=0;cu<nbCTotal;cu++)
     {
     cellsList.push_back(intPair(unsigned long long ints[cu],cu));
     cout<< " cu  = "<<cu<< " image cell = "<<unsigned long long ints[cu]<<endl;
     }*/
    ////////system("pause");
    /*!
     * Toutes les images sont calcul�es  et stock�es dans un tableau dans l'ordre
     *  de num�rotation des controles. La tache suivante consiste � trier ce tableau, �liminer les doublons
     *   et  compl�ter la structure  image de point
     */

    cellsList.sort();
    list<unsigned long long int>::iterator itCell, itDouble;

    itCell = cellsList.begin();
    itDouble = itCell;
    unsigned long long int currentCell;

    cu = 0;
    unsigned long long int iCellsTab = 0;
    //cout<< " juste avant while de parcours de la liste  first "<<(*itCell).first<< " nb cells total "<<(int)nbCellsTotal<<endl;
    while ((itCell != cellsList.end()) & ((*itCell) <= nbCellsTotal + 1))
	{

	currentCell = (*itCell);

	//    cout<< " current cell "<<currentCell<<endl;

	pointDI.tabImageCells[iCellsTab] = currentCell;

	iCellsTab++;
	//cout<< " itCellsTab="<<iCellsTab<<endl;
	while ((itDouble != cellsList.end()) & ((*itDouble) == currentCell))
	    {
	    itDouble++;
	    }
	/*!
	 * la boucle s'arrete au premier different
	 * � la fin de la boucle soit itDouble a atteint la fin de la liste
	 *  soit  il pointe sur une cellule diff�rente
	 */
	itCell = itDouble;
	}

    pointDI.nbImageCells = iCellsTab;
    //cout<< "  nb  de mailles differentes dans l'iage d'un point "<<iCellsTab<<endl;
    /*!
     * La derni�re valeur  du tableau des controles indique la fin  de la liste des controles associ�s
     * � la derni�re cellule
     */

    }

void ViabiBitSet::computeDiscreteImageOfPoint_noControl(
	double *doublePointCoords, unsigned long long int *intPointCoords)
    {

    //double testV[dim];
    //unsigned long long int testI[dim];

    //cout<< " calcul de l'image d'un point  SANS controle\n";

    double **controlCoords = dynsys->getControlCoords();
    unsigned long long int nbCellsTotal = grid->getNbTotalPoints();

    list<unsigned long long int> cellsList;

    //cout<< " calcul de l'image pour le point ";
    // printVector(doublePointCoords,dim);
    //double dv[dim];

    /*!
     * on calcule les coordonn�es r�elles de l'image du point par la dynamique discrete
     * elles sont stockes dans le tableau doubleVect
     * si \f$ u\in U(x)\f$  (contr�le admissible) on calcule l'image r�elle par la dynamique
     *  discr�tis�e  en temps  du point x avec le controle u
     */

    // printVector(doublePointCoords,dim);
    (dynsys->*(dynsys->discretDynamics))(doublePointCoords, controlCoords[0],
	    doubleVect1, 1.0);

    // cout<< " retrour dnamique discrete ";
    //printVector(doubleVect1, dim);
    //(dynsys->*discretDynamics)(doublePointCoords, controlCoords[cu], doubleVect);
    if (grid->isPointInGrid(doubleVect1))
	{
	//printf( "   le point est das la grlle\n " );

	//////printf(" le point est das la grlle\n");
	/*!
	 * Si l'image est  dans les limites de la grille on �tudie si elle v�rifie les contraintes
	 */
	if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
	    {
	    //     printf( "   contraintes Xu ok\n " );

	    //     ////printf("  contraintes sur   X  ok\n");

	    /*!
	     * Si l'image  est dans l'ensemble de contraintes sur l'�tat \f$ K \f$
	     * on calcule le num�ro de maille qui contient cett image
	     */
	    imageCells[0] = grid->localizePoint(doubleVect1);
	    /* cout<<  " num de cell image "<<imageCells[0]<<endl;
	     cout<< " int coords de d�part ";
	     for(int lm=0;lm<dim;lm++)
	     {
	     cout<< " "<<intPointCoords[lm];
	     }
	     cout<< endl;
	     cout<< "double coords de image ";
	     for(int lm=0;lm<dim;lm++)
	     {
	     cout<< " "<<doubleVect1[lm];
	     }*/
	    /*    grid->numToIntAndDoubleCoords(unsigned long long ints[cu],testI,testV);
	     if(testV[0]<doublePointCoords[0])
	     {
	     cout<< " num de d�parrt "<<num<< " num de cell image "<<unsigned long long ints[cu]<<endl;
	     cout<< " int coords de d�part ";
	     for(int lm=0;lm<dim;lm++)
	     {
	     cout<< " "<<intPointCoords[lm];
	     }
	     cout<< endl;
	     cout<< "double coords de d�part ";
	     for(int lm=0;lm<dim;lm++)
	     {
	     cout<< " "<<doublePointCoords[lm];
	     }
	     cout<< endl;
	     cout<< " double coords de result ";
	     for(int lm=0;lm<dim;lm++)
	     {
	     cout<< " "<<doubleVect1[lm];
	     }
	     cout<< endl;
	     cout<< " int coords de result projet� ";
	     for(int lm=0;lm<dim;lm++)
	     {
	     cout<< " "<<testI[lm];
	     }
	     cout<< endl;
	     cout<< "double coords de result projet�  ";
	     for(int lm=0;lm<dim;lm++)
	     {
	     cout<< " "<<testV[lm];
	     }
	     cout<< endl;

	     }*/
	    // on enregistre le numero de maille
	    cellsList.push_back(imageCells[0]);
	    //////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,unsigned long long ints[cu]);
	    //   cout<< " les coordonn�es  de cin bas de la maille  sont ";
	    //   grid->numToIntAndDoubleCoords(unsigned long long ints[cu], intPointCoords, doubleVect1);
	    //   printVector(doubleVect1, dim);
	    }
	else
	    {
	    /*!
	     * Si l'image n'est pas dans  \f$ K \f$ on enregistre un num�ro de maille factice qui signifie que
	     * cette image est rejet�e
	     */
	    imageCells[0] = nbCellsTotal + 10; // sinon on enregistre un nombre convenu reconnaissanble
	    cellsList.push_back(imageCells[0]);
	    //    ////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,unsigned long long ints[cu]);
	    }
	}
    else
	{
	/*!
	 * Si l'image n'est pas dans  la grille on enregistre un num�ro de maille factice qui signifie que
	 * cette image est rejet�e
	 */
	if (grid->isPointInGridWithConstr(doubleVect1))
	    {
	    imageCells[0] = nbCellsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble
	    cellsList.push_back(imageCells[0]);

	    }
	else
	    {
	    imageCells[0] = nbCellsTotal + 10; // sinon on enregistre un nombre convenu reconnaissanble
	    cellsList.push_back(imageCells[0]);
	    }
	//////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,unsigned long long ints[cu]);
	}

    // cout<< " ajour d'une paire � la liste longueur "<<cellsList.size()<<"\n";

    /*    cout<< " fini traitement parallele\n";
     for( cu=0;cu<nbCTotal;cu++)
     {
     cellsList.push_back(intPair(unsigned long long ints[cu],cu));
     cout<< " cu  = "<<cu<< " image cell = "<<unsigned long long ints[cu]<<endl;
     }*/
    ////////system("pause");
    /*!
     * Toutes les images sont calcul�es  et stock�es dans un tableau dans l'ordre
     *  de num�rotation des controles. La tache suivante consiste � trier ce tableau, �liminer les doublons
     *   et  compl�ter la structure  image de point
     */

    cellsList.sort();
    list<unsigned long long int>::iterator itCell, itDouble;

    itCell = cellsList.begin();
    itDouble = itCell;
    unsigned long long int currentCell;

    unsigned long long int iCellsTab = 0;
    //  cout<< " juste avant while de parcours de la liste taille de liste  "<<cellsList.size()<< " nb cells total "<<nbCellsTotal<<endl;
    while ((itCell != cellsList.end()) & ((*itCell) <= nbCellsTotal + 1))
	{

	currentCell = (*itCell);

	//   cout<< " current cell "<<currentCell<<endl;

	pointDI.tabImageCells[iCellsTab] = currentCell;

	iCellsTab++;
	// cout<< " itCellsTab="<<iCellsTab<<endl;
	while ((itDouble != cellsList.end()) & ((*itDouble) == currentCell))
	    {
	    itDouble++;
	    }
	/*!
	 * la boucle s'arrete au premier different
	 * � la fin de la boucle soit itDouble a atteint la fin de la liste
	 *  soit  il pointe sur une cellule diff�rente
	 */
	itCell = itDouble;
	// exit(1);
	}

    pointDI.nbImageCells = iCellsTab;
    // cout<< "  nb  de mailles differentes dans l'iage d'un point "<<iCellsTab<<endl;
    /*!
     * La derni�re valeur  du tableau des controles indique la fin  de la liste des controles associ�s
     * � la derni�re cellule
     */

    }

void ViabiBitSet::computeTrajectories()
    {
    computeViableTrajectories();
    }

void ViabiBitSet::computeViableTrajectories()
    {
    algoViabiParams *avp = modelParams->getAlgoParameters();
    int nbTrajs = avp->NB_TRAJS;
    int typeTraj = avp->TYPE_TRAJ;
    double T = modelParams->getSystemParameters()->maxTime;
    ostringstream os;
    string fileName;

    if (nbTrajs > 0)
	{
	for (int tr = 0; tr < nbTrajs; tr++)
	    {
	    if (typeTraj == VD)
		{
		os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1
			<< ".dat";
		fileName = os.str();
		os.str("");

		computeViableTrajectorySetVal(avp->INIT_POINTS + tr * dim, T,
			fileName);
		}
	    if (typeTraj == VL)
		{
		os << "../OUTPUT/" << filePrefix << "-traj-H-" << tr + 1
			<< ".dat";
		fileName = os.str();
		os.str("");

		computeViableTrajectoryHeavy(avp->INIT_POINTS + tr * dim,
			avp->INIT_CONTROLS + tr * dimC, T, fileName);
		}
	    }
	}
    }

void ViabiBitSet::CaptureBasin()
    {
    algoViabiParams *avp = modelParams->getAlgoParameters();
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
    os << "../OUTPUT/" << filePrefix << "-Capture.dat";
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
	os << "../OUTPUT/" << filePrefix << "-CaptureSlice.dat";
	fileName = os.str();
	os.str("");
	grid->saveCoupe(fileName);
	}
    if (avp->SAVE_SLICE_BOUND)
	{
	os << "../OUTPUT/" << filePrefix << "-CaptureSliceBound" << ".dat";
	fileName = os.str();
	os.str("");
	grid->saveCoupeBoundary(fileName);
	}
    if (avp->SAVE_BOUNDARY)
	{
	os << "../OUTPUT/" << filePrefix << "-Capture-bound.dat";
	fileName = os.str();
	os.str("");
	grid->saveBoundary(fileName);
	}

    if (avp->SAVE_PROJECTION)
	{
	os << "../OUTPUT/" << filePrefix << "-Capture-proj" << ".dat";
	fileName = os.str();
	os.str("");
	/*
	 *  calcul et sauvegarde  de la projection du  noyau
	 */
	grid->saveProjetion(fileName, avp->PROJECTION);
	}

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
    system("pause");
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
void ViabiBitSet::createPointsList()
    {

    currentImagePointsList.pointsList.clear();
    currentImagePointsList.maxNum = -1;
    currentImagePointsList.minNum = currentImageList.minNum;

    int nbPointsCube = (int) pow(2.0, dim);	//pow(2.0, dim);
    list<unsigned long long int>::iterator itCell =
	    currentImageList.cellsList.begin(), itLast =
	    currentImageList.cellsList.end();
    int i;

    long long int *indicesDecalCell = grid->getIndicesDecalCell();

    unsigned long long int numCell;
    unsigned long long int posX;

    imagePoint tempPoint;
    list<unsigned long long int>::iterator itStart, itNew;
    itStart = this->currentImagePointsList.pointsList.begin();

    while (itCell != itLast)
	{
	numCell = (*itCell);

	//cout<< "  maille num�ro "<<numCell<<endl;
	for (i = 0; i < nbPointsCube; i++)
	    {
	    posX = numCell + indicesDecalCell[i];

	    /*!
	     * \todo proc�der comme pour la cr�ation de la liste de cellules:
	     * cr�er un point avec le num�ro et les donn�es de la cellule.
	     * copier une seule fois les donn�es de la cellules sur la valeur et la r�tro-action
	     *  d'un point � l'autre de la m�me cellule  seul le num�ro du point change.
	     *
	     *  Puis inserrer le point dans la liste ordonn�e de points.
	     *   Exactement comme pour les celllules! Donc du copier coller de code.
	     *
	     */
	    //cout<< " on ajoute le point numero "<<posX<<endl;
	    this->addDataToPointsList(&itStart, posX, &itNew);
	    itStart = itNew;

	    }
	itCell++;
	}

    }

void ViabiBitSet::addDataToPointsList(
	list<unsigned long long int>::iterator *startIt,
	unsigned long long int newPoint,
	list<unsigned long long int>::iterator *resIt)
    {
    /*!
     * Dans l'image  en construction de \f$ C_{n}\setminus C_{n-1}\f$ les mailles sont stock�es dans l'ordre croissant
     * de leur num�ros  et chaque maille garde la m�moire de la valeur minimale ainsi que de tous les ant�c�dents
     * de cette maille c'est � dire tous les couples viables (x,u) pour lesquels f(x,u) appartient � cette maille.
     * En plus  on distingue les couples (x,u) minimisant la valeur  et les autres, juste viables.
     *
     * Pour ajouter une nouvelle maille dans cette liste  on doit faire dans l'ordre les choses suivantes :
     *
     *    rechercher dans la liste la place nu num�ro de maille � inserrer
     *   deux cas de figure peuvent se pr�senter :
     *     la maille ayant le m�me num�ro  existe d�j�: on proc�de alors � la fusion des deux,  en d�terminant la valeur optimale et
     *    en  fusionnant les r�tro-actions
     *     la maille n'existe pas encore; alors on l'inserre simplement  dans l'ordre croissant.
     *
     */
    //cout<< " l'ajout du point itStart pointe sur "<<(*(*startIt)).PointNum <<endl;
    unsigned long long int numnewPoint = newPoint;

    list<unsigned long long int>::iterator itCell;

    if (currentImagePointsList.pointsList.size() == 0)
	{
	currentImagePointsList.pointsList.push_back(newPoint);
	currentImagePointsList.maxNum = numnewPoint;
	currentImagePointsList.minNum = numnewPoint;
	(*resIt) = currentImagePointsList.pointsList.end();
	(*resIt)--;
	//cout<< " on ajoute le premier element nouveax min et max  de la liste sont "<<currentImagePointsList.minNum<< " "<<currentImagePointsList.maxNum<<endl;
	}
    else
	{
	if (numnewPoint
		> (unsigned long long int) currentImagePointsList.maxNum)
	    {
	    currentImagePointsList.pointsList.push_back(newPoint);
	    currentImagePointsList.maxNum = numnewPoint;
	    (*resIt) = currentImagePointsList.pointsList.end();
	    (*resIt)--;
	    //cout<< " on ajoute � la fin  nouveax min et max  de la liste sont "<<currentImagePointsList.minNum<< " "<<currentImagePointsList.maxNum<<endl;
	    }
	else
	    {
	    if (numnewPoint
		    < (unsigned long long int) currentImagePointsList.minNum)
		{
		currentImagePointsList.pointsList.push_front(newPoint);
		currentImagePointsList.minNum = numnewPoint;
		(*resIt) = currentImagePointsList.pointsList.begin();
		//    cout<< " on ajoute au debut nouveax min et max  de la liste sont "<<currentImagePointsList.minNum<< " "<<currentImagePointsList.maxNum<<endl;

		}
	    else
		{
		itCell = *startIt;

		if (numnewPoint < (*itCell))
		    {

		    while ((numnewPoint < (*itCell)))
			{
			//cout<< " current cell num "<<(*itCell).PointNum<< " new cell num "<<numnewPoint<<endl;
			itCell--;
			}
		    if (numnewPoint > (*itCell))
			{
			itCell++;
			currentImagePointsList.pointsList.insert(itCell,
				newPoint);
			}

		    (*resIt) = itCell;

		    }
		else
		    {
		    while ((numnewPoint > (*itCell)))
			{
			//      cout<< " current point num "<<(*itCell).PointNum<< " new point  num "<<numnewPoint<<endl;
			itCell++;
			}
		    if (numnewPoint < (*itCell))
			{
			currentImagePointsList.pointsList.insert(itCell,
				newPoint);
			}

		    (*resIt) = itCell;

		    }
		}
	    }
	}
    //cout<< " l'ajout du point est termin�  on a ajout� juste avant "<<(*(*resIt)).PointNum<<endl;
    ////////system("pause");
    }

void ViabiBitSet::computeCurrentImage(int iter)
    {

    cout << "  Compute Current Image . Iteration num " << iter
	    << " points dans Cn = " << currentImagePointsList.pointsList.size()
	    << endl;

    double t1, t2, elapsed_time;

    timeval tim;
    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

    int posX;

    unsigned long long int intPointCoords[dim];
    double doublePointCoords[dim];

    list<unsigned long long int>::iterator itStart, itNew;

    /*!
     * \todo  introduire � ce niveau un pointeur sur la liste  des cellules images
     *  ce pointeur devra suivre l'insertion des cellules
     *  puisque elles sont dans l'ordre croissant on
     *  recherchera la suivant � partir du pointeur sur la derbi�re  qui a �t� ajout�e
     *
     *
     *  Aussi il faut faire l'insertion  de tout le bloc  des donn�es
     *  pour la m�me cellule
     *  donc former tout de m�me une cellule  et apr�s l'ajouter dans la liste
     *  ok
     */

    currentImageList.cellsList.clear();
    currentImageList.maxNum = -1;
    currentImageList.minNum = PLUS_INF;

    list<unsigned long long int>::iterator itPoint =
	    currentImagePointsList.pointsList.begin(), itLastPoint =
	    currentImagePointsList.pointsList.end(), itTemp;
    itStart = this->currentImageList.cellsList.begin();

    while (itPoint != itLastPoint)
	{
	posX = (*itPoint);
	//cout<< " posX="<<posX<<endl;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);

	this->computeDiscreteImageOfPoint(doublePointCoords, intPointCoords);
	/*!
	 * L'image calcul�e est stock�e dans la structure pointDI sous forme de tableau ordonn�
	 *  de cellules
	 */
	for (unsigned long long int i = 0; i < pointDI.nbImageCells; i++)
	    {
	    addDataToCurrentImage(&itStart, pointDI.tabImageCells[i], &itNew);
	    //cout<< " ares ajout  de l'image d'un point la taille de la liste est "<<currentImageList.cellsList.size()<<endl;
	    itStart = itNew;
	    }

	itPoint++;

	}

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));

    cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;
    }
int ViabiBitSet::addNewPoints()
    {
    /*!
     * Cette fonction  r�alise l'ajout  dans la base de donn�es  de la nouvelle couche
     * \f$ \Phi(C_{n}\setminus C_{n-1})\f$  dans la base de donn�es.
     * Elle doit �galement aouter aux bases de donn�es de r�troaction optimale et viable
     * les r�troactions des points calcul�s
     *
     * Les points sont stock�s sous forme de liste de structures imagePoint ordonn�e par num�ro de point x.
     *  Chaque structure contient le num�ro du  point, la fonction valeur et
     *  les listes de triplets \f$ (x,u,\rho)\f$ repr�sentant les r�trocations du point, optimale et viable.
     *
     *  Pour chaque point  de la liste on v�rifie  d'abord s'il est d�j�  dans la base de donn�es. S'il n'y est pas,
     *  on l'ajoute � l'ensemble en construction et en m�me temps on ajoute ses r�troactions aux deux bases de r�troaction
     *  correspondantes.
     *
     *  Si le point existe d�j� dans l'ensemble, il poss�de d�j� une r�tro-action �gaement. Dans ce cas, on
     *  v�rifie les fonctions valeurs. Si celle du nuveau point est inf�rieure, alors on modifie le point existant
     *  et la r�troaction optimale. Sinon,  on modifie  selement les r�troactons viables.
     *
     *  Important!  Il est � noter ici que c'est cette m�me liste de points qui doit �tre ajout�e  � la base
     *  par cette fonction qui servira ensuite � la construction de la couche suivante.
     *  Ainsi, si un point de cette liste exste d�j� dans la base et que sa fonction valeur de la base n'est pas modifi�e,
     *  il n'est pas consid�r� comme nouveau, il n'appartient pas � \f$ C_{n+1}\setminus C_n\f$ . On doit donc le supprimer
     *   de la liste pour ne pas  calculer son image plutard.
     *
     */

    cout << " ajout de nouveaux points\n";

    double t1, t2, elapsed_time;

    timeval tim;
    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    unsigned long long int intCoords[dim];

    int nbNewPoints = 0;

    list<unsigned long long int>::iterator itPoint =
	    currentImagePointsList.pointsList.begin(), itLastPoint =
	    currentImagePointsList.pointsList.end(), itTemp;
    list<triple>::iterator itLastR, itR;
    unsigned long long int pointNum;

    while (itPoint != itLastPoint)
	{

	pointNum = (*itPoint);
	grid->numToIntCoords(pointNum, intCoords);
	if (!grid->isInSet(intCoords))
	    {

	    //            cout << " le point est nouveau on l'ajout dans la base  de fonc val \n";
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
	    itTemp = itPoint;
	    itPoint++;
	    currentImagePointsList.pointsList.erase(itTemp);

	    }

	}

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));

    cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;
    return nbNewPoints;
    }

void ViabiBitSet::noyauViabi_sansControle(bool sortieOK, int nbArret)
    {

    int dirTramage = grid->getDirTram();

    unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    bool testF;

    //    cout<<"masque points enleves cree"<<masquePointsEnleves;
    unsigned long long int longTrame = grid->getLongTrame();
    boost::dynamic_bitset<> masque;
    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(
	    longTrame, 0);
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
			xCoordsDouble[j] = limInf[j]
				+ indice[j - 1] * gridStep[j];
			}
		    ///cout<< " Analyse en cours "<<endl;
		    for (unsigned long long int k = 0; k < longTrame; k++)
			{
			if (masque[k])
			    {
			    xCoordsDouble[dirTramage] = limInf[dirTramage]
				    + k * gridStep[dirTramage];
			    testNonVide = this->findViabImagePoint_noControl(
				    xCoordsDouble, false);
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
		    if (masquePointsEnleves->count()
			    < (unsigned long long int) longTrame)
			{

			*gridTab[posX] &= (*masquePointsEnleves);
			}

		    }

		}			//fin de if la trame n'est pas vide

	    }				// fin de for de parcours de la trame

	cout << "Itération " << nbIter
		<< " terminée. Nombre de points  points enlevés: "
		<< comptEnleves << "\n";

	nbIter++;

	}
    cout << "fini nbIter=" << nbIter;
    //foncCarNoyau->printTrame();

    }

void ViabiBitSet::noyauViabi_sansControle_omp(bool sortieOK, int nbArret)
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

		unsigned long long int *nbPointsSub =
			grid->getNbPointsSubGrid();

		unsigned long long int longTrame = grid->getLongTrame();

		double *limInf = grid->limInf;
		double *gridStep = grid->step;

		bool testNonVide;
		//printf("thread numero %d nouvelle boucle while\n", tid);
		boost::dynamic_bitset<> masque;
		boost::dynamic_bitset<> *masquePointsEnleves =
			new boost::dynamic_bitset<>(longTrame, 0);

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
			    xCoordsDouble[j] = limInf[j]
				    + indice[j] * gridStep[j];
			    }

			for (int j = dirTramage + 1; j < dim; j++)
			    {
			    xCoordsDouble[j] = limInf[j]
				    + indice[j - 1] * gridStep[j];
			    }
			for (unsigned long long int k = 0; k < longTrame; k++)
			    {
			    if (masque[k])
				{
				xCoordsDouble[dirTramage] = limInf[dirTramage]
					+ k * gridStep[dirTramage];
				testNonVide =
					this->findViabImagePoint_noControl(
						xCoordsDouble, false);

				if (!testNonVide)
				    {

				    masquePointsEnleves->set(k, false);
				    comptEnleves++;
				    }
				}			// fin de if masque[k]
			    }		// fin de for  de parcours de masque

			if (masquePointsEnleves->count()
				< (unsigned long long int) longTrame)
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

	cout << "Itération " << nbIter
		<< " terminée. Nombre de points  points enlevés: "
		<< comptEnleves << "\n";
	nbIter++;

	}
    cout << "fini nbIter=" << nbIter;
    }

void ViabiBitSet::GarantedViabilityKernel(bool sortieOK, int nbArret)
    {

    }

SysDyn* ViabiBitSet::GetSysDynForViabProblem()
    {
    return this->dynsys;
    }

void ViabiBitSet::loadViableSets()
    {

    algoViabiParams *avp = modelParams->getAlgoParameters();
    int refine = avp->GRID_REFINMENTS_NUMBER;

    ostringstream os;
    string fileName;
    // on charge dans la mémoire l'ensemble calculé et enregistré
    // correspondant au dernier raffinement
    os << "../OUTPUT/" << filePrefix << "-viab-" << refine << ".dat";
    fileName = os.str();
    os.str("");
    grid->loadSet(fileName);

    if (avp->SAVE_SLICE)
	{
	os << "../OUTPUT/" << filePrefix << "-viab-" << refine << "-Slice"
		<< ".dat";
	fileName = os.str();
	os.str("");
	grid->saveCoupe(fileName);
	}
    if (avp->SAVE_SLICE_BOUND)
	{
	os << "../OUTPUT/" << filePrefix << "-viab-" << refine << "-SliceBound"
		<< ".dat";
	fileName = os.str();
	os.str("");
	grid->saveCoupeBoundary(fileName);
	}

    if (avp->SAVE_PROJECTION)
	{
	os << "../OUTPUT/" << filePrefix << "-viab-" << refine << "-proj"
		<< ".dat";
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
    ostringstream os;
    string fileName;
    algoViabiParams *avp = modelParams->getAlgoParameters();
    int refine = avp->GRID_REFINMENTS_NUMBER;
    os << "../OUTPUT/" << filePrefix << "-viab-" << refine << ".dat";
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
	os << "../OUTPUT/" << filePrefix << "-viab-" << refine << "-Slice"
		<< ".dat";
	fileName = os.str();
	os.str("");
	grid->saveCoupe(fileName);
	}
    if (avp->SAVE_SLICE_BOUND)
	{
	os << "../OUTPUT/" << filePrefix << "-viab-" << refine << "-SliceBound"
		<< ".dat";
	fileName = os.str();
	os.str("");
	grid->saveCoupeBoundary(fileName);
	}
    if (avp->SAVE_BOUNDARY)
	{
	os << "../OUTPUT/" << filePrefix << "-viab-" << refine << "-bound.dat";
	fileName = os.str();
	os.str("");
	grid->saveBoundary(fileName);
	}

    if (avp->SAVE_PROJECTION)
	{
	os << "../OUTPUT/" << filePrefix << "-viab-" << refine << "-proj"
		<< ".dat";
	fileName = os.str();
	os.str("");
	/*
	 *  calcul et sauvegarde  de la projection du  noyau
	 */
	grid->saveProjetion(fileName, avp->PROJECTION);
	}
    }

void ViabiBitSet::ViabilityKernel(bool sortieOK, int nbArret)
    {
    algoViabiParams *avp = modelParams->getAlgoParameters();
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
	cout
		<< "=============================Debut viab kernel =============================="
		<< endl;
	ViabilityKernelSimple(1, seuilArret);
	cout
		<< "=============================================================================="
		<< endl;
	gettimeofday(&tim, NULL);
	t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	elapsed_time = (double) ((t2 - t1));
	cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;
	//   cout<< "=============================================================================="<<endl;
	/*  * sauvegardes
	 */
	if (avp->INTERMEDIATE_SAVINGS)
	    {
	    os << "../OUTPUT/" << filePrefix << "-viab-" << refIter + 1
		    << ".dat";
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
		os << "../OUTPUT/" << filePrefix << "-viab-" << refIter + 1
			<< "-Slice" << ".dat";
		fileName = os.str();
		os.str("");
		grid->saveCoupe(fileName);
		}
	    if (avp->SAVE_SLICE_BOUND)
		{
		os << "../OUTPUT/" << filePrefix << "-viab-" << refIter + 1
			<< "-SliceBound" << ".dat";
		fileName = os.str();
		os.str("");
		grid->saveCoupeBoundary(fileName);
		}
	    if (avp->SAVE_BOUNDARY)
		{
		os << "../OUTPUT/" << filePrefix << "-viab-" << refIter + 1
			<< "-bound.dat";
		fileName = os.str();
		os.str("");
		grid->saveBoundary(fileName);
		}

	    if (avp->SAVE_PROJECTION)
		{
		os << "../OUTPUT/" << filePrefix << "-viab-" << refIter + 1
			<< "-proj" << ".dat";
		fileName = os.str();
		os.str("");
		/*
		 *  calcul et sauvegarde  de la projection du  noyau
		 */
		grid->saveProjetion(fileName, avp->PROJECTION);
		}
	    }

	refIter++;
	seuilArret *= 2;
	/*
	 * raffinement du maillage
	 */

	if (refIter < refine)
	    {
	    cout << "refine\n";
	    grid->refine();

	    os << "../OUTPUT/" << filePrefix << "-viab-Refined" << refIter
		    << ".dat";
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

	    if (avp->SAVE_BOUNDARY)
		{
		os << "../OUTPUT/" << filePrefix << "-viab-Refined" << refIter
			<< "-bound.dat";
		fileName = os.str();
		os.str("");
		grid->saveBoundary(fileName);
		}
	    }
	}
    if (!avp->INTERMEDIATE_SAVINGS)
	{
	os << "../OUTPUT/" << filePrefix << "-viab-" << refIter << ".dat";
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
	    os << "../OUTPUT/" << filePrefix << "-viab-" << refine << "-Slice"
		    << ".dat";
	    fileName = os.str();
	    os.str("");
	    grid->saveCoupe(fileName);
	    }
	if (avp->SAVE_SLICE_BOUND)
	    {
	    os << "../OUTPUT/" << filePrefix << "-viab-" << refine
		    << "-SliceBound" << ".dat";
	    fileName = os.str();
	    os.str("");
	    grid->saveCoupeBoundary(fileName);
	    }
	if (avp->SAVE_BOUNDARY)
	    {
	    os << "../OUTPUT/" << filePrefix << "-viab-" << refine
		    << "-bound.dat";
	    fileName = os.str();
	    os.str("");
	    grid->saveBoundary(fileName);
	    }

	if (avp->SAVE_PROJECTION)
	    {
	    os << "../OUTPUT/" << filePrefix << "-viab-" << refine << "-proj"
		    << ".dat";
	    fileName = os.str();
	    os.str("");
	    /*
	     *  calcul et sauvegarde  de la projection du  noyau
	     */
	    grid->saveProjetion(fileName, avp->PROJECTION);
	    }
	}

    }
void ViabiBitSet::ViabilityKernelSimple(bool sortieOK, int nbArret)
    {
    if (dynsys->getDimC() == 0)
	{
	if (nbOMPThreads > 1)
	    {
	    noyauViabi_sansControle_omp(sortieOK, nbArret);
	    }
	else
	    {
	    cout << "viab kernel sans controle" << endl;
	    noyauViabi_sansControle(sortieOK, nbArret);
	    }
	}
    else
	{
	if ((dynsys->getDynType() == CC) || (dynsys->getDynType() == DC))
	    {
	    if (nbOMPThreads > 1)
		{
		noyauViabi_omp(sortieOK, nbArret);
		}
	    else
		{
		noyauViabi(sortieOK, nbArret);
		}
	    }
	else
	    {
	    if ((dynsys->getDynType() == DD))
		{
		noyauViabi_FD(sortieOK, nbArret);
		}
	    }
	}
    }

void ViabiBitSet::noyauViabi(bool sortieOK, int nbArret)
    {

    int dirTramage = grid->getDirTram();

    unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    bool testF;

    unsigned long long int longTrame = grid->getLongTrame();
    boost::dynamic_bitset<> masque;
    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(
	    longTrame, 0);
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
			xCoordsDouble[j] = limInf[j]
				+ indice[j - 1] * gridStep[j];
			}

		    for (unsigned long long int k = 0; k < longTrame; k++)
			{
			if (masque[k])
			    {
			    xCoordsDouble[dirTramage] = limInf[dirTramage]
				    + k * gridStep[dirTramage];
			    testNonVide = this->findViabImagePoint(
				    xCoordsDouble, false);

			    if (!testNonVide)
				{
				masquePointsEnleves->set(k, false);
				comptEnleves++;
				}
			    }	// fin de if masque[k]
			}	// fin de for  de parcours de masque

		    if (masquePointsEnleves->count()
			    < (unsigned long long int) longTrame)
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

	cout << "Itération " << nbIter
		<< " terminée. Nombre de points  points enlevés: "
		<< comptEnleves << "\n";

	nbIter++;

	}
    cout << "fini nbIter=" << nbIter;
    //foncCarNoyau->printTrame();

    }

void ViabiBitSet::noyauViabi_omp(bool sortieOK, int nbArret)
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

		unsigned long long int *nbPointsSub =
			grid->getNbPointsSubGrid();

		unsigned long long int longTrame = grid->getLongTrame();

		double *limInf = grid->limInf;
		double *gridStep = grid->step;
		// printf("thread numero %d nouvelle boucle while\n", nbTh);

		if (!gridTab[posX]->none())
		    {
		    // cout<<  (*gridTab[posX])<<endl;

		    boost::dynamic_bitset<> masque =
			    grid->analyseTrameMasqueBis(posX, false);

		    //   cout<<"masque d'analyse  "<<masque<<endl;

		    if (!masque.none())
			{
			bool testNonVide;

			boost::dynamic_bitset<> *masquePointsEnleves =
				new boost::dynamic_bitset<>(longTrame, 0);
			numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

			masquePointsEnleves->set();

			for (int j = 0; j < dirTramage; j++)
			    {
			    xCoordsDouble[j] = limInf[j]
				    + indice[j] * gridStep[j];
			    }

			for (int j = dirTramage + 1; j < dim; j++)
			    {
			    xCoordsDouble[j] = limInf[j]
				    + indice[j - 1] * gridStep[j];
			    }

			for (unsigned long long int k = 0; k < longTrame; k++)
			    {

			    if (masque[k])
				{
				xCoordsDouble[dirTramage] = limInf[dirTramage]
					+ k * gridStep[dirTramage];
				testNonVide = this->findViabImagePoint(
					xCoordsDouble, false);

				//cout<< " fini\n";
				if (!testNonVide)
				    {
				    masquePointsEnleves->set(k, false);
				    comptEnleves++;
				    }
				}			// fin de if masque[k]
			    }		// fin de for  de parcours de masque

			if (masquePointsEnleves->count()
				< (unsigned long long int) longTrame)
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

	cout << "Itération " << nbIter
		<< " terminée. Nombre de points  points enlevés: "
		<< comptEnleves << "\n";

	nbIter++;

	}
    cout << "Calcul de noyau fini. Nb total d'itérations: " << nbIter;

    }

void ViabiBitSet::noyauViabiGaranti_FD(bool sortieOK, int nbArret)
    {
    int dirTramage = grid->getDirTram();

    unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

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
    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(
	    longTrame, 0);
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
	int dirParcours = -1;
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

		    if (!((masque.count() == 2)
			    && (masque[0] && masque[longTrame - 1])))
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
				    if (dynsys->constraintsXU_fd(coordDiscretes,
					    controlIntCoords[compteComm])<PLUS_INF)
					{
					// cout<< "  controle "<<controlIntCoords[compteComm][0]<<" "<<controlIntCoords[compteComm][1]<< " est autorise \n";

					testAllTych = true;
					for (compteTych = 0; compteTych < nbTy;
						compteTych++)
					    {
					    //   cout<< "  tychet "<<tychIntCoords[compteTych][0]<<endl;
					    dynsys->dynamics_tych_fd(
						    coordDiscretes,
						    controlIntCoords[compteComm],
						    tychIntCoords[compteTych],
						    imageXU);
					    // cout<<  " image ";
					    //                                             printVector(imageXU,dim);

					    if (grid->isPointInGrid_fd(imageXU))
						{
						testAllTych &= grid->isInSet(
							imageXU);

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

			if (masquePointsEnleves->count()
				< (unsigned long long int) longTrame)
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
			    if (dynsys->constraintsXU_fd(coordDiscretes,
				    controlIntCoords[compteComm]) < PLUS_INF)
				{
				testAllTych = true;
				for (compteTych = 0; compteTych < nbTy;
					compteTych++)
				    {
				    dynsys->dynamics_tych_fd(coordDiscretes,
					    controlIntCoords[compteComm],
					    tychIntCoords[compteTych], imageXU);
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
			    if (dynsys->constraintsXU_fd(coordDiscretes,
				    controlIntCoords[compteComm]) < PLUS_INF)
				{
				testAllTych = true;
				for (compteTych = 0; compteTych < nbTy;
					compteTych++)
				    {
				    dynsys->dynamics_tych_fd(coordDiscretes,
					    controlIntCoords[compteComm],
					    tychIntCoords[compteTych], imageXU);
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

			if (masquePointsEnleves->count()
				< (unsigned long long int) longTrame)
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

	unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

	unsigned long long int tailleTrame = grid->getNbPointsTotalSubGrid();
	//	cout<<"masque points enleves cree"<<masquePointsEnleves;
	unsigned long long int longTrame = grid->getLongTrame();

	unsigned long long int nbC = dynsys->getTotalNbPointsC();
	boost::dynamic_bitset<> masque;
	boost::dynamic_bitset<> **gridTab = grid->getGridTab();

	unsigned long long int indice[dim - 1];

	unsigned long long int **controlIntCoords =
		dynsys->getControlIntCoords();

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
			grid->intCoordsToNum(coordDiscretes, &posIm);
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
			    if (dynsys->constraintsXU_fd(coordDiscretes,
				    controlIntCoords[compteComm]) < PLUS_INF)
				{
				//cout<< "  controle "<<controlIntCoords[compteComm][0]<<" "<<controlIntCoords[compteComm][1]<< " est autorise \n";

				dynsys->dynamics_fd(coordDiscretes,
					controlIntCoords[compteComm], imageXU);
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
				grid->intCoordsToNum(imageXU, &posIm);
				//	 	cout<< " numero image "<<posIm<<endl;
				viabControls.push_back(compteComm);
				images.push_back(posIm);
				}
			    compteComm++;
			    }
			// system("pause");
			fichierB << images.size() << " ";
			for (unsigned long long int jj = 0; jj < images.size();
				jj++)
			    {
			    fichierB << viabControls.at(jj) << " "
				    << images.at(jj) << " ";
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

	unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

	int tailleTrame = grid->getNbPointsTotalSubGrid();
	//	cout<<"masque points enleves cree"<<masquePointsEnleves;
	unsigned long long int longTrame = grid->getLongTrame();

	unsigned long long int nbC = dynsys->getTotalNbPointsC();
	unsigned long long int nbTy = dynsys->getTotalNbPointsTy();

	boost::dynamic_bitset<> masque;
	boost::dynamic_bitset<> **gridTab = grid->getGridTab();
	unsigned long long int indice[dim - 1];

	unsigned long long int **controlIntCoords =
		dynsys->getControlIntCoords();
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
			grid->intCoordsToNum(coordDiscretes, &posIm);
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
			    if (dynsys->constraintsXU_fd(coordDiscretes,
				    controlIntCoords[compteComm]) < PLUS_INF)
				{
				//cout<< "  controle "<<controlIntCoords[compteComm][0]<<" "<<controlIntCoords[compteComm][1]<< " est autorise \n";

				testAllTych = true;
				for (compteTych = 0; compteTych < nbTy;
					compteTych++)
				    {
				    dynsys->dynamics_tych_fd(coordDiscretes,
					    controlIntCoords[compteComm],
					    tychIntCoords[compteTych], imageXU);
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
				for (compteTych = 0; compteTych < nbTy;
					compteTych++)
				    {
				    dynsys->dynamics_tych_fd(coordDiscretes,
					    controlIntCoords[compteComm],
					    tychIntCoords[compteTych], imageXU);
				    grid->intCoordsToNum(imageXU, &posIm);
				    images.push_back(posIm);
				    }

				}
			    compteComm++;
			    }
			// system("pause");
			fichierB << viabControls.size() << " ";
			int cptIm = 0;
			for (unsigned long long int jj = 0;
				jj < viabControls.size(); jj++)
			    {
			    fichierB << viabControls.at(jj) << " ";
			    for (compteTych = 0; compteTych < nbTy;
				    compteTych++)
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

