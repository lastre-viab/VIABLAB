/*
 * ViabiHJB.cpp
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
 *  Created on: 9 december 2013
 *      Author: Anna DESILLES
 */

#include "../include/ViabiMicroMacro.h"

ViabiMicroMacro::ViabiMicroMacro(ParametersManager *pm)
    {
    modelParams = pm;

    systemParams sp = *(pm->getSystemParameters());

    algoViabiParams avp = *(pm->getAlgoParameters());

    controlParams cp = *(pm->getControlParameters());
    gridParams gp = *(pm->getGridParameters());

    dim = gp.DIM;
    grid = new GridMicroMacro(gp);
    grid->printGrid();

    dynsys = new SysDyn(sp, dim, cp, grid);
    InitViabiMicroMacro(avp);
    trajectoryHelper = new ViabiMicroMacroTrajectoryHelper(grid, dynsys, avp.TYPE_TRAJ);
    }

SysDyn*
ViabiMicroMacro::GetSysDynForViabProblem()
    {
    return this->dynsys;
    }

GridMicroMacro*
ViabiMicroMacro::GetGridForViabProblem()
    {
    return grid;
    }

void ViabiMicroMacro::InitViabiMicroMacro(algoViabiParams avp)
    {

    ostringstream os;

    vTab = grid->getGridPtr();

    dim = grid->dim;
    dimC = dynsys->getDimC();

    nbOMPThreads = avp.NB_OMP_THREADS;

    spdlog::info("Initialization of Viability algorithm : Micro-Macro model");

    /*!
     * On initialise la structure  qui servira au stockage temporaire de l'image discrete d'un point de grille
     * Les tableaux sont initialisés  avec leur taille maximale possible :
     * égale au nb  de controles.
     * C'est la valeur effective  de nbImageCells , calculé é chaque evaluation de l'image discrete
     * d'un point qui  servira  é lire et remplir  correctement
     * la bonne partie de ces tableaux
     */
    pointDI.tabCellEntrees = new unsigned long long int[dynsys->getTotalNbPointsC() + 1];
    pointDI.tabImageCells = new unsigned long long int[dynsys->getTotalNbPointsC()];
    pointDI.tabImageControls = new unsigned long long int[dynsys->getTotalNbPointsC()];

    /*
     * Initialisation de tableaux servant de variables globales pour certaines fonctions
     * Cela évte de multiples allocations/destructions de mémoire
     */
    intPointCoords = new unsigned long long int[dim];
    intVect1 = new unsigned long long int[dim];
    doublePointCoords = new double[dim];
    doubleVect = new double[dim];
    doubleVect1 = new double[dim];
    intControlCoords = new unsigned long long int[dimC];
    doubleControlCoords = new double[dimC];
    imageCells = new unsigned long long int[dynsys->getTotalNbPointsC()];

    filePrefix = avp.FILE_PREFIX;

    targ_or_dep = avp.TARGET_OR_DEPARTURE;
    computeTmin = avp.COMPUTE_TMIN;

    ViabiMicroMacro::addNewPointsToSet = &ViabiMicroMacro::addNewPoints;
    ViabiMicroMacro::addDataToCurrentCell = &ViabiMicroMacro::addDataToCell;
    ViabiMicroMacro::addDataToCurrentPoint = &ViabiMicroMacro::addDataToPoint;


	ViabiMicroMacro::createCurrentPointsList = &ViabiMicroMacro::createPointsList;

    if (avp.COMPUTE_TMIN)
	{
	ViabiMicroMacro::computeFirstConvexifiedImage = &ViabiMicroMacro::computeConvexifiedImage_tmin;
	ViabiMicroMacro::computeCurrentImage = &ViabiMicroMacro::computeCurrIm_tmin;
	//  ViabiMicroMacro::computeFirstConvexifiedImage_omp=&ViabiMicroMacro::computeConvexifiedImage_tmin_omp;
	}
    else
	{

	    ViabiMicroMacro::computeCurrentImage = &ViabiMicroMacro::computeCurrIm_Lmin;
	    ViabiMicroMacro::computeFirstConvexifiedImage = &ViabiMicroMacro::computeConvexifiedImage_Lmin;

	}

    tempPointsList1 = new list<imagePoint>();
    tempPointsList2 = new list<imagePoint>();
    whichPointListToUse = 1;
    currentImagePointsList.pointsList = tempPointsList1;
    }

void ViabiMicroMacro::computeDiscreteImageOfPoint(unsigned long long int num)
    {
    double **controlCoords = dynsys->getControlCoords();
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int nbCellsTotal = grid->getNbTotalPoints();
    grid->numToIntAndDoubleCoords(num, intPointCoords, doublePointCoords);
    double rho;
    try
    {
	rho = dynsys->calculRho_local(doublePointCoords);
    }
    catch (const std::bad_alloc &e)
	{
	std::cout << "Allocation failed dans calculRho: " << e.what() << '\n';
	cout << " doublePointCoords = ";
	printVector(doublePointCoords, dim);
	exit(1);
	}
    unsigned long long int cu;

    list<intPair> cellsList = list<intPair>();
    for (cu = 0; cu < nbCTotal; cu++)
	{

	/*!
	 * on calcule les coordonnées réelles de l'image du point par la dynamique discrete
	 * elles sont stockes dans le tableau doubleVect
	 * si \f$ u\in U(x)\f$  (contréle admissible) on calcule l'image réelle par la dynamique
	 *  discrétisée  en temps  du point x avec le controle u
	 */
	if (dynsys->constraintsXU(doublePointCoords, controlCoords[cu]) < PLUS_INF)
	    {
	    try
	    {
		(dynsys->*(dynsys->discretDynamics))(doublePointCoords, controlCoords[cu], doubleVect1, rho);
	    }
	    catch (const std::bad_alloc &e)
		{
		std::cout << "Allocation failed dans discrete dyn: " << e.what() << '\n';
		cout << " doublePointCoords = ";
		printVector(doublePointCoords, dim);
		exit(1);
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
		     * on calcule le numéro de maille qui contient cett image
		     */
		    imageCells[cu] = grid->localizePoint(doubleVect1);
		    // on enregistre le numero de maille
		    cellsList.push_back(intPair(imageCells[cu], cu));
		    }
		else
		    {

		    /*!
		     * Si l'image n'est pas dans  \f$ K \f$ on enregistre un numéro de maille factice qui signifie que
		     * cette image est rejetée
		     */
		    imageCells[cu] = nbCellsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble
		    cellsList.push_back(intPair(imageCells[cu], cu));
		    }
		}
	    else
		{
		/*!
		 * Si l'image n'est pas dans  la grille on enregistre un numéro de maille factice qui signifie que
		 * cette image est rejetée
		 */
		/*!
		 * \todo prévoir la gestion des autorisations de sortie  par axe pour les variables qui peuvent
		 * se trouver dans des intervalles non bornés
		 */
		if (grid->unboundedDomain && grid->isPointInGridWithConstr(doubleVect1) && (dynsys->constraintsX(doubleVect1) < PLUS_INF))
		    {
		    imageCells[cu] = nbCellsTotal; // sinon on enregistre un nombre convenu reconnaissanble
		    cellsList.push_back(intPair(imageCells[cu], cu));
		    }
		else
		    {
		    imageCells[cu] = nbCellsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble
		    cellsList.push_back(intPair(imageCells[cu], cu));
		    }
		}
	    }
	else
	    {
	    /*!
	     * Si l'image n'est pas dans  la grille on enregistre un numéro de maille factice qui signifie que
	     * cette image est rejetée
	     */
	    /*!
	     * \todo prévoir la gestion des autorisations de sortie  par axe pour les variables qui peuvent
	     * se trouver dans des intervalles non bornés
	     */
	    imageCells[cu] = nbCellsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble
	    cellsList.push_back(intPair(imageCells[cu], cu));
	    }
	}
    for (cu = 0; cu < nbCTotal; cu++)
	{
	cellsList.push_back(intPair(imageCells[cu], cu));
	}
    /*!
     * Toutes les images sont calculées  et stockées dans un tableau dans l'ordre
     *  de numérotation des controles. La tache suivante consiste é trier ce tableau, éliminer les doublons
     *   et  compléter la structure  image de point
     */

    cellsList.sort(pairCompare);
    list<intPair>::iterator itCell, itDouble;

    itCell = cellsList.begin();
    itDouble = itCell;
    int currentCell;

    cu = 0;
    unsigned long long int iControlTab = 0, iCellsTab = 0, iEntreesTab = 0;
    while ((itCell != cellsList.end()) & ((*itCell).first <= (int) nbCellsTotal))
	{

	currentCell = (*itCell).first;
	if (this->testConstraintesForCell(currentCell))
	    {
	    pointDI.tabCellEntrees[iEntreesTab] = iControlTab;
	    pointDI.tabImageCells[iCellsTab] = currentCell;
	    iEntreesTab++;
	    iCellsTab++;
	    while ((itDouble != cellsList.end()) & ((*itDouble).first == currentCell))
		{
		pointDI.tabImageControls[iControlTab] = (*itDouble).second;
		iControlTab++;
		itDouble++;

		}
	    /*!
	     * la boucle s'arrete au premier different
	     * é la fin de la boucle soit itDouble a atteint la fin de la liste
	     *  soit  il pointe sur une cellule différente
	     */

	    itCell = itDouble;
	    }
	else
	    {
	    while ((itDouble != cellsList.end()) & ((*itDouble).first == currentCell))
		{

		itDouble++;
		}
	    /*!
	     * la boucle s'arrete au premier different
	     * é la fin de la boucle soit itDouble a atteint la fin de la liste
	     *  soit  il pointe sur une cellule différente
	     */

	    itCell = itDouble;
	    }

	}

    pointDI.nbImageCells = iCellsTab;

    /*!
     * La derniére valeur  du tableau des controles indique la fin  de la liste des controles associés
     * é la derniére cellule
     */
    pointDI.tabCellEntrees[iCellsTab] = iControlTab;
    }


void ViabiMicroMacro::initialiseTarget()
    {
    /*!
     *  cette fonction initialise la base de données pour la dynamique. Elle
     *  ajoute dnas la base les premiers points pour lesquels la fonction valeur
     *  est réelle. Au début de l'algorithme de bassin de capture
     *    seuls les points de la  cible  ont une fonction valeur réelle
     *
     */

	initialiseTargetHJB();
    }

void ViabiMicroMacro::computeTrajectories()
    {

    dynsys->setDynamicsForward();

    algoViabiParams *avp = modelParams->getAlgoParameters();
    if (avp->TYPE_TRAJ == OP)
	computeOptimalTrajectories();
    else
	computeViableTrajectories();

    }

void ViabiMicroMacro::computeOptimalTrajectories()
    {
    algoViabiParams *avp = modelParams->getAlgoParameters();
    int nbTrajs = avp->NB_TRAJS;
    ostringstream os;
    string fileName;
    if (nbTrajs > 0)
	{
	bool success[nbTrajs];
	for (int tr = 0; tr < nbTrajs; tr++)
	    {
	    os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1 << ".dat";
	    fileName = os.str();
	    os.str("");

	    computeOptimalCaptTrajectory(avp->INIT_POINTS + tr * dim, fileName, success[tr]);
	    }
	}
    }

void ViabiMicroMacro::computeViableTrajectories()
    {
    algoViabiParams *avp = modelParams->getAlgoParameters();
    int nbTrajs = avp->NB_TRAJS;
    int typeTraj = avp->TYPE_TRAJ;
    ostringstream os;
    string endOfFileName = ".dat";
    switch (typeTraj)
	{
    case VD:
	{
	endOfFileName = "-viabDefault.dat";
	break;
	}
    case VDI:
	{
	endOfFileName = "-viabDiffControls.dat";
	break;
	}
    case VMM:
	{
	endOfFileName = "-viabMinValue.dat";
	break;
	}
    case VG:
	{
	endOfFileName = "-viabGaranti.dat";
	break;
	}
    default:
	{
	endOfFileName = "-viabDefault.dat";
	break;
	}
	}

    string fileName;
    if (nbTrajs > 0)
	{
	spdlog::debug(" Number of trajectories {}", nbTrajs);
	bool success[nbTrajs];
	for (int tr = 0; tr < nbTrajs; tr++)
	    {
	    os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1 << endOfFileName;
	    fileName = os.str();
	    os.str("");


		trajectoryHelper->computeViableTrajectory(avp->INIT_POINTS + tr * dim, avp->INIT_VALUES[tr], fileName, success[tr]);
	    }
	}
    }

void ViabiMicroMacro::initialiseTargetHJB()
    {
    /*
     *  cette fonction initialise la base de données pour . Elle
     *  ajoute dnas la base les premiers points pour lesquels la fonction valeur
     *
     *   est réelle. Au début de l'algorithme de bassin de capture
     *    seuls le spoints de la  cible  ont une fonction valeur réelle
     *
     */

    spdlog::info("Target set initialization");
    double c;

    unsigned long long int dim = grid->dim;

    unsigned long long int *x = new unsigned long long int[dim];

    double *xReel = new double[dim];

    int totalPointsC = 0;  // nombre de points  de l'espace des commandes

    int totalPointsX = grid->getNbTotalPoints();

    imagePoint currentPoint;

    currentImagePointsList.maxNum = 0;
    currentImagePointsList.minNum = 0;
    currentImagePointsList.pointsList = tempPointsList1;
    ;

    unsigned long long int pos;

    list<imagePoint>::iterator itStart = currentImagePointsList.pointsList->begin(), itNew;
    int cpt = 0;
    /*
     *  on parcourt  tous les points de l'espace discret  fini
     *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
     */
    for (pos = 0; pos < (unsigned long long int) totalPointsX; pos++)
	{

	/*!
	 * le compteur pos  ets l'unique numéro entier du point en cours
	 * dans la numérotation alphabétique : on parcourt axe par axe
	 */

	/*!
	 * on restitue les  coordonnées netiéres  du point é partir de son numéro
	 * ainsi que ses coordonnées réelles
	 */
	grid->numToIntAndDoubleCoords(pos, x, xReel);

	c = max((*(dynsys->target))(xReel), (*(dynsys->constraintsX))(xReel));

	vTab[pos] = c;
	if (c < PLUS_INF)
	    {
	    totalPointsC++;
	    currentPoint.minVal = c;
	    currentPoint.PointNum = pos;

	    addDataToPointsList(&itStart, currentPoint, &itNew);
	    itStart = itNew;
	    cpt++;
	    }
	}
    spdlog::info("Target set initialized, number of grid points in target set : {0:d}", cpt);
    }


void ViabiMicroMacro::loadViableSets()
    {

    spdlog::info("Loading viable sets for the problem with prefix {}", filePrefix);
    algoViabiParams *avp = modelParams->getAlgoParameters();

    ostringstream os;
    string fileName;
    // on charge dans la mémoire l'ensemble calculé et enregistré
    // correspondant au dernier raffinement
    spdlog::info("Value function file {}", fileName);
    os << "../OUTPUT/" << filePrefix << "-valFunc.dat";
    fileName = os.str();
    os.str("");
    grid->loadSet(fileName);

    if (avp->SAVE_PROJECTION)
	{
	os << "../OUTPUT/" << filePrefix << "-proj" << ".dat";
	fileName = os.str();
	os.str("");
	/*
	 *  calcul et sauvegarde  de la projection du  noyau
	 */
	grid->saveProjetion(fileName, avp->PROJECTION);
	}
    }

void ViabiMicroMacro::saveViableSets()
    {
    saveValFunctions();
    }
void ViabiMicroMacro::saveValFunctions()
    {
    algoViabiParams *avp = modelParams->getAlgoParameters();

    ostringstream os;
    string fileName;

	if (avp->SAVE_SUBLEVEL)
	    {
	    os << "../OUTPUT/" << filePrefix << "-subLevel.dat";
	    fileName = os.str();
	    os.str("");
	    grid->saveSubLevelset(avp->LEVEL, fileName);
	    }

	if (avp->SAVE_BOUNDARY)
	    {
	    os << "../OUTPUT/" << filePrefix << "-valFunc.dat";
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
	    }

	if (avp->SAVE_PROJECTION)
	    {
	    os << "../OUTPUT/" << filePrefix << "-proj" << ".dat";
	    fileName = os.str();
	    os.str("");
	    /*
	     *  calcul et sauvegarde  de la projection du  noyau
	     */
	    grid->saveProjetion(fileName, avp->PROJECTION);
	    }

	if (avp->SAVE_SLICE_BOUND)
	    {
	    os << "../OUTPUT/" << filePrefix << "-slice" << ".dat";
	    fileName = os.str();
	    os.str("");
	    /*
	     *  calcul et sauvegarde  de la projection du  noyau
	     */
	    grid->saveCoupeBoundary(fileName);
	    }


    }

void ViabiMicroMacro::viabKerValFunc(unsigned long long int nbArret)
    {
    unsigned long long int iCoords[dim];
    unsigned long long int cptChanged = 20000000, nbIter = 0;
    unsigned long long int pos;
    unsigned long long int nbCellsTotal = grid->getNbTotalPoints();
    unsigned long long int compteComm, cellNum, numControl;

    double rCoords[dim];
    spdlog::info("Viable kernel computation fo continuous dynamical system with Micro-Macro model");
    double rho;

    uintPair image;

    int totalPointsX = grid->getNbTotalPoints();

    long long int *indicesDecalCell = grid->getIndicesDecalCell();
    int nbPointsCube = (int) pow(2.0, dim);

    double **controlCoords = dynsys->getControlCoords();
    double minValCell, valAtPos;

    while ((cptChanged > nbArret) & (nbIter < 15000))
	{
	cptChanged = 0;
	/*
	 *  on parcourt  tous les points de l'espace discret  fini
	 *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
	 */
	for (pos = 0; pos < (unsigned long long int) totalPointsX; pos++)
	    {
	    double tempVal, tempL, tempM;
	    valAtPos = vTab[pos];
	    if (valAtPos < PLUS_INF)
		{
		grid->numToIntAndDoubleCoords(pos, iCoords, rCoords);
		try
		{
		    rho = dynsys->calculRho_local(rCoords);
		}
		catch (const std::bad_alloc &e)
		    {
		    std::cout << "Allocation failed dans calculRho: " << e.what() << '\n';
		    cout << " pos = " << pos << " rcoords = ";
		    printVector(rCoords, dim);
		    exit(1);
		    }

		bool testNonVide = false;

		try
		{
		    this->computeDiscreteImageOfPoint(pos);

		}
		catch (const std::bad_alloc &e)
		    {
		    std::cout << "Allocation failed dans computeDiscrete inmage: " << e.what() << '\n';
		    exit(1);
		    }

		compteComm = 0;
		minValCell = PLUS_INF;

		while (compteComm < pointDI.nbImageCells && !testNonVide)
		    {
		    cellNum = pointDI.tabImageCells[compteComm];
		    if (cellNum < nbCellsTotal)
			{

			for (int iCell = 0; iCell < nbPointsCube - 1; iCell++)
			    {

			    double imageVal = vTab[cellNum + indicesDecalCell[iCell]];
			    for (unsigned long long int j = pointDI.tabCellEntrees[compteComm];
				    j < pointDI.tabCellEntrees[compteComm + 1] && !testNonVide; j++)
				{
				numControl = pointDI.tabImageControls[j];

				tempL = dynsys->lFunc(rCoords, controlCoords[numControl]);

				tempM = dynsys->mFunc(rCoords, controlCoords[numControl]);
				if (tempL < PLUS_INF && tempM < PLUS_INF)
				    {
				    tempVal = (imageVal + rho * tempL) / (1 - rho * tempM);
				    testNonVide = tempVal < valAtPos;
				    minValCell = min(minValCell, tempVal);

				    }
				if (testNonVide)
				    break;
				}
			    if (testNonVide)
				break;
			    }
			}
		    else if (cellNum == nbCellsTotal)
			{
			for (unsigned long long int j = pointDI.tabCellEntrees[compteComm]; j < pointDI.tabCellEntrees[compteComm + 1]; j++)
			    {
			    numControl = pointDI.tabImageControls[j];
			    (dynsys->*(dynsys->discretDynamics))(rCoords, controlCoords[numControl], doubleVect1, 1.0);
			    printVector(doubleVect1, dim);
			    tempL = dynsys->lFunc(rCoords, controlCoords[numControl]);
			    tempM = dynsys->mFunc(rCoords, controlCoords[numControl]);

			    double tempV = dynsys->constraintsX(doubleVect1);
			    tempVal = (tempV + rho * tempL) / (1 - rho * tempM);
			    testNonVide = tempVal < valAtPos;
			    minValCell = min(minValCell, tempVal);
			    if (testNonVide)
				break;
			    }
			}

		    compteComm++;
		    }

		if (vTab[pos] < minValCell)
		    {
		    cptChanged++;

		    }
		vTab[pos] = max(vTab[pos], minValCell);
		}
	    }

	nbIter++;
	spdlog::info("Iteration {} finished. Number of points with updated value funtion : {}", nbIter, cptChanged);
	}
    spdlog::info("Viability kernel computation finished. next step : saving results");
    saveValFunctions();
    }

void ViabiMicroMacro::viabKerValFunc_omp(unsigned long long int nbArret)
    {

    unsigned long long int cptChanged = 20000000, nbIter = 0;
    spdlog::info("This is a OMP parallelized viability algorithm for micro-macro model. Number of threads : {}", nbOMPThreads);

    int totalPointsX = grid->getNbTotalPoints();

    double *gridTab = grid->getGridPtr();
    double *gridTabNew = grid->getGridPtr_tmp();
    grid->copyGrid(gridTab, gridTabNew);

    while ((cptChanged > nbArret) & (nbIter < 15000))
	{
	cptChanged = 0;
	/*
	 *  on parcourt  tous les points de l'espace discret  fini
	 *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
	 */
	unsigned long long int pos = 0;
#pragma omp parallel for num_threads(nbOMPThreads)  reduction(+:cptChanged) private(pos)  shared( gridTab, gridTabNew, totalPointsX) default(none)
	for (pos = 0; pos < (unsigned long long int) totalPointsX; pos++)
	    {

	    double **controlCoords = dynsys->getControlCoords();
	    double minValCell, valAtPos;

	    double rCoords[dim];
	    long long int *indicesDecalCell = grid->getIndicesDecalCell();
	    int nbPointsCube = (int) pow(2.0, dim);
	    unsigned long long int nbCellsTotal = grid->getNbTotalPoints();
	    unsigned long long int iCoords[dim];
	    unsigned long long int compteComm, cellNum, numControl;
	    double rho;
	    double *doubleVect1 = new double[dim];
	    uintPair image;
	    double tempVal, tempL, tempM;
	    valAtPos = gridTab[pos];
	    if (valAtPos < PLUS_INF)
		{
		grid->numToIntAndDoubleCoords(pos, iCoords, rCoords);
		try
		{
		    rho = dynsys->calculRho_local(rCoords);
		}
		catch (const std::bad_alloc &e)
		    {

		    printVector(rCoords, dim);
		    exit(1);
		    }

		bool testNonVide = false;
		try
		{
		    this->computeDiscreteImageOfPoint(pos);
		}
		catch (const std::bad_alloc &e)
		    {

		    printVector(rCoords, dim);
		    exit(1);
		    }

		compteComm = 0;
		minValCell = PLUS_INF;

		while (compteComm < pointDI.nbImageCells && !testNonVide)
		    {
		    cellNum = pointDI.tabImageCells[compteComm];
		    if (cellNum < nbCellsTotal)
			{

			for (int iCell = 0; iCell < nbPointsCube - 1; iCell++)
			    {
			    double imageVal = gridTab[cellNum + indicesDecalCell[iCell]];

			    for (unsigned long long int j = pointDI.tabCellEntrees[compteComm];
				    j < pointDI.tabCellEntrees[compteComm + 1] && !testNonVide; j++)
				{
				numControl = pointDI.tabImageControls[j];
				tempL = dynsys->lFunc(rCoords, controlCoords[numControl]);
				tempM = dynsys->mFunc(rCoords, controlCoords[numControl]);
				if (tempL < PLUS_INF && tempM < PLUS_INF)
				    {
				    tempVal = (imageVal + rho * tempL) / (1 - rho * tempM);
				    testNonVide = tempVal < valAtPos;
				    minValCell = min(minValCell, tempVal);

				    }
				if (testNonVide)
				    break;
				}
			    if (testNonVide)
				break;
			    }
			}
		    else if (cellNum == nbCellsTotal)
			{
			for (unsigned long long int j = pointDI.tabCellEntrees[compteComm]; j < pointDI.tabCellEntrees[compteComm + 1]; j++)
			    {
			    numControl = pointDI.tabImageControls[j];
			    (dynsys->*(dynsys->discretDynamics))(rCoords, controlCoords[numControl], doubleVect1, 1.0);
			    tempL = dynsys->lFunc(rCoords, controlCoords[numControl]);
			    tempM = dynsys->mFunc(rCoords, controlCoords[numControl]);

			    double tempV = dynsys->constraintsX(doubleVect1);
			    tempVal = (tempV + rho * tempL) / (1 - rho * tempM);
			    testNonVide = tempVal < valAtPos;
			    minValCell = min(minValCell, tempVal);
			    if (testNonVide)
				break;
			    }
			}

		    compteComm++;
		    }

		if (gridTab[pos] < minValCell)
		    {
		    cptChanged++;
		    }
		gridTabNew[pos] = max(gridTab[pos], minValCell);
		}
	    delete[] doubleVect1;
	    }
	grid->copyGrid(gridTabNew, gridTab);
	nbIter++;
	spdlog::info("Iteration {} finished. Number of points with updated value funtion : {}", nbIter, cptChanged);
	}
    spdlog::info("Viability kernel computation finished. next step : saving results");
    saveValFunctions();
    }

void ViabiMicroMacro::CaptureBasin()
    {
    spdlog::info("Capture bassin algorithm : start");

    dynsys->setDynamicsBackward();
    /*!
     * \var nbNewPoints : nb de nouveaux points ajoutés é l'étape n
     */
    int nbNewPoints = 1;

    int iter = 0;

    /*!
     * On calcule la premiére itération, en tanant compte de la convexification de la dynamique sur la cible
     * On appelle pour cela la fonction computeConvexifiedImage().
     */
    (this->*computeFirstConvexifiedImage)(iter);

    /*!
     * -A partir de la liste des mailles on crée une liste ordonnée de points représentant l'image é l'aide de la fonction
     * createPointsList().
     */
    (this->*createCurrentPointsList)();
    /*!
     * - A partir de la liste des points représentant l'image  trois bases  de données sont alimentées: La base
     * principale, associée é la grille  et représentant la fonction valeur, la base de rétroaction optimale et la base de
     *  rétro-action viable. On appelle pour cela la fonction addNewPoints().
     */
    if (!(dynsys->getDynType() == DD))
	{
	nbNewPoints = (this->*addNewPointsToSet)();
	}

    iter++;

    /*!
     * Tant qu'il y a de nouveaux points ajoutés on répéte les opérations suivantes.
     */
    while ((nbNewPoints > 0))
	{
	spdlog::info("Iteration {} : new points in the set {}", iter, nbNewPoints);

	/*!
	 * -Calcul de l'image de \f^C_{n}\setminus C_n\f$  par la fonction computeCurrentImageLocalRho(): l'image est enregistrée en émoire vive sous forme de liste
	 * de références de mailles dans lesquelles arrive au moins une évolution. Chaque référence de maille contient
	 * des informations sur tous les antécédants de cette maille ainsi  que la valeur minimale de
	 * temps. Dans cette version oé \f$ \rho\f$ est global, la fonction valeur prend la méme valeur
	 * é chaque étape : \f$ \rho \cdot n \f$.
	 */
	(this->*computeCurrentImage)(iter);
	//this->showCurrentImageList();

	/*!
	 * -A partir de la liste des mailles on crée une liste ordonnée de points représentant l'image en appelant la fonction createPointsList().
	 *  Comme pour le mailles, chaque référence de point regroupe les informatons (regroupées é partir de différentes
	 * mailles dont est vertex)  sur les antécédents  de ce point. Le but de la création de cette liste est d'éliminer
	 * les doublons afin de minimiser les accés é la base de données
	 */
	(this->*createCurrentPointsList)();
	//this->showCurrentImagePointsList();
	/*!
	 * - A partir de la liste des points représentant l'image  trois bases  de données sont alimentées: La base
	 * principale, associée é la grille  et représentant la fonction valeur, la base de rétroaction optimale et la base de
	 *  rétro-action viable. On appelle ici la fonction addNewPoints().
	 */
	nbNewPoints = (this->*addNewPointsToSet)();
	iter++;
	}

    saveValFunctions();
    }

double ViabiMicroMacro::computeOptimalCaptTrajectory(double *initPosition, string fileName, bool &succes)
    {
    return trajectoryHelper->computeOptimalTrajectory(initPosition, fileName, succes);
    /* if (computeTmin)
    		return trajectoryHelper->computeOptimalTrajectory(initPosition,
    				fileName, succes);
    	else
    		return trajectoryHelper->computeOptimalTrajectory_Lmin(initPosition,
    				fileName, succes);*/
    }

void ViabiMicroMacro::computeCurrIm_tmin(int iter)
    {
    double t1, t2, elapsed_time;

    timeval tim;
    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

    int posX;

    list<imageCell>::iterator itStart, itNew;

    /*!
     * \todo  introduire à ce niveau un pointeur sur la liste  des cellules images
     *  ce pointeur devra suivre l'insertion des cellules
     *  puisque elles sont dans l'ordre croissant on
     *  recherchera la suivant à partir du pointeur sur la derbière  qui a été ajoutée
     *
     *
     *  Aussi il faut faire l'insertion  de tout le bloc  des données
     *  pour la méme cellule
     *  donc former tout de méme une cellule  et aprés l'ajouter dans la liste
     *  ok
     */

    currentImageList.cellsList.clear();
    currentImageList.maxNum = -1;
    currentImageList.minNum = 1000 + grid->nbTotalCells;
    double rho;

    list<imagePoint>::iterator itPoint = currentImagePointsList.pointsList->begin(), itTemp;
    imageCell tempImageCell;

    itStart = this->currentImageList.cellsList.begin();

    while (!currentImagePointsList.pointsList->empty())	//(itPoint!=itLastPoint)
	{

	posX = (*itPoint).PointNum;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	rho = dynsys->calculRho_local(doublePointCoords);
	{
	    /*!
	     * On calcule l'image discrète  du point
	     */
	    tempImageCell.minVal = (*itPoint).minVal + rho;

	    this->computeDiscreteImageOfPoint(posX);
	    /*!
	     * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
	     *  de cellules
	     */

	    for (unsigned long long int i = 0; i < pointDI.nbImageCells; i++)
		{
		tempImageCell.cellNum = pointDI.tabImageCells[i];
		addDataToCurrentImage(&itStart, tempImageCell, &itNew);
		itStart = itNew;
		}
	}
	itPoint++;
	currentImagePointsList.pointsList->pop_front();
	}
    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));
    spdlog::debug("[Min Time problem] : Elapsed time to compute current image {} sec", elapsed_time);
    }

void ViabiMicroMacro::computeCurrIm_Lmin(int iter)
    {

    double t1, t2, elapsed_time;

    timeval tim;
    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

    int posX;

    unsigned long long int nbCellsTotal = grid->getNbTotalPoints();

    double **controlCoords = dynsys->getControlCoords();
    list<imageCell>::iterator itStart, itNew;

    /*!
     * \todo  introduire à ce niveau un pointeur sur la liste  des cellules images
     *  ce pointeur devra suivre l'insertion des cellules
     *  puisque elles sont dans l'ordre croissant on
     *  recherchera la suivant à partir du pointeur sur la derbière  qui a été ajoutée
     *
     *
     *  Aussi il faut faire l'insertion  de tout le bloc  des données
     *  pour la méme cellule
     *  donc former tout de méme une cellule  et aprés l'ajouter dans la liste
     *  ok
     */

    currentImageList.cellsList.clear();
    currentImageList.maxNum = -1;
    currentImageList.minNum = 1000 + grid->nbTotalCells;
    double rho, tempL;

    list<imagePoint>::iterator itPoint = currentImagePointsList.pointsList->begin(), itTemp;
    imageCell tempImageCell;

    itStart = this->currentImageList.cellsList.begin();
    double imageCoords[dim];

    while (!currentImagePointsList.pointsList->empty())	//(itPoint!=itLastPoint)
	{

	posX = (*itPoint).PointNum;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	rho = dynsys->calculRho_local(doublePointCoords);

	/*!
	 * On calcule l'image discrète  du point
	 */

	this->computeDiscreteImageOfPoint(posX);
	/*!
	 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
	 *  de cellules
	 */
	unsigned long long int numCell, numControl;
	for (unsigned long long int i = 0; i < pointDI.nbImageCells; i++)
	    {
	    numCell = pointDI.tabImageCells[i];
	    tempImageCell.cellNum = numCell;

	    if (numCell < nbCellsTotal)
		{
		tempImageCell.minVal = PLUS_INF;
		for (unsigned long long int j = pointDI.tabCellEntrees[i]; j < pointDI.tabCellEntrees[i + 1]; j++)
		    {
		    numControl = pointDI.tabImageControls[j];
		    tempL = dynsys->lFunc(doublePointCoords, controlCoords[numControl]);
		    (dynsys->*(dynsys->discretDynamics))(doublePointCoords, controlCoords[numControl], imageCoords, 1.0);

		    double tempL1 = dynsys->lFunc(imageCoords, controlCoords[numControl]);
		    tempImageCell.minVal = min(tempImageCell.minVal, (*itPoint).minVal + rho * 0.5 * (tempL + tempL1));
		    }
		addDataToCurrentImage(&itStart, tempImageCell, &itNew);
		itStart = itNew;
		}
	    }

	itPoint++;
	currentImagePointsList.pointsList->pop_front();
	}

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));

    spdlog::debug("[Integral cost problem] : Elapsed time to compute current image {} sec", elapsed_time);

    }


void ViabiMicroMacro::computeConvexifiedImage_Lmin_omp(int iter)
    {
    int posX;
    list<imageCell>::iterator itStart, itNew;

    currentImageList.cellsList.clear();
    currentImageList.maxNum = -1;
    currentImageList.minNum = 1000 + grid->nbTotalCells;
    double rho, tempL;

    list<imagePoint>::iterator itPoint = currentImagePointsList.pointsList->begin(), itLastPoint = currentImagePointsList.pointsList->end(), itTemp;
    imageCell tempImageCell;

    unsigned long long int nbCellsTotal = grid->getNbTotalPoints();
    double **controlCoords = dynsys->getControlCoords();

    itStart = this->currentImageList.cellsList.begin();

    while (itPoint != itLastPoint)
	{
	posX = (*itPoint).PointNum;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	rho = dynsys->calculRho_local(doublePointCoords);
	/*!
	 * On calcule l'image discrète  du point
	 */
	this->computeDiscreteImageOfPoint(posX);
	/*!
	 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
	 *  de cellules
	 */
	unsigned long long int numCell, numControl;
	for (unsigned long long int i = 0; i < pointDI.nbImageCells; i++)
	    {
	    numCell = pointDI.tabImageCells[i];
	    tempImageCell.cellNum = numCell;

	    if (numCell < nbCellsTotal)
		{
		tempImageCell.minVal = PLUS_INF;
		for (unsigned long long int j = pointDI.tabCellEntrees[i]; j < pointDI.tabCellEntrees[i + 1]; j++)
		    {
		    numControl = pointDI.tabImageControls[j];
		    tempL = dynsys->lFunc(doublePointCoords, controlCoords[numControl]);
		    tempImageCell.minVal = min(tempImageCell.minVal, (*itPoint).minVal + rho * tempL);
		    }
		addDataToCurrentImage(&itStart, tempImageCell, &itNew);
		itStart = itNew;
		addConvexCombinations(itPoint, numCell, &tempImageCell, rho, &itStart);
		}
	    }
	itPoint++;
	}
    }

void ViabiMicroMacro::computeConvexifiedImage_tmin(int iter)
    {
    spdlog::info("[Min Time problem] : Computing convexified image of teh target set. number of points {}",
	    currentImagePointsList.pointsList->size());
    int posX;

    list<imageCell>::iterator itStart, itNew;

    currentImageList.cellsList.clear();
    currentImageList.maxNum = -1;
    currentImageList.minNum = 1000 + grid->nbTotalCells;
    double rho;

    list<imagePoint>::iterator itPoint = currentImagePointsList.pointsList->begin(), itLastPoint = currentImagePointsList.pointsList->end(), itTemp;
    imageCell tempImageCell;

    itStart = this->currentImageList.cellsList.begin();

    while (itPoint != itLastPoint)
	{
	posX = (*itPoint).PointNum;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	rho = dynsys->calculRho_local(doublePointCoords);
	{
	    /*!
	     * On calcule l'image discrète  du point
	     */
	    tempImageCell.minVal = (*itPoint).minVal + rho;
	    this->computeDiscreteImageOfPoint(posX);
	    /*!
	     * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
	     *  de cellules
	     */
	    unsigned long long int numCell;

	    for (unsigned long long int i = 0; i < pointDI.nbImageCells; i++)
		{
		numCell = pointDI.tabImageCells[i];
		tempImageCell.cellNum = numCell;
		addDataToCurrentImage(&itStart, tempImageCell, &itNew);
		itStart = itNew;
		addConvexCombinations(itPoint, numCell, &tempImageCell, rho, &itStart);

		tempImageCell.minVal = (*itPoint).minVal + rho;
		}
	}
	itPoint++;

	}
    }


void ViabiMicroMacro::computeConvexifiedImage_Lmin(int iter)
    {
    spdlog::info("[Integral cost problem] : Computing convexified image of the target set. number of points {}",
	    currentImagePointsList.pointsList->size());

    int posX;
    list<imageCell>::iterator itStart, itNew;

    currentImageList.cellsList.clear();
    currentImageList.maxNum = -1;
    currentImageList.minNum = 1000 + grid->nbTotalCells;
    double rho, tempL;

    list<imagePoint>::iterator itPoint = currentImagePointsList.pointsList->begin(), itLastPoint = currentImagePointsList.pointsList->end(), itTemp;
    imageCell tempImageCell;

    unsigned long long int nbCellsTotal = grid->getNbTotalPoints();

    double **controlCoords = dynsys->getControlCoords();

    itStart = this->currentImageList.cellsList.begin();
    double imageCoords[dim];
    while (itPoint != itLastPoint)
	{
	posX = (*itPoint).PointNum;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	rho = dynsys->calculRho_local(doublePointCoords);
	/*!
	 * On calcule l'image discrète  du point
	 */
	this->computeDiscreteImageOfPoint(posX);
	/*!
	 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
	 *  de cellules
	 */
	unsigned long long int numCell, numControl;

	for (unsigned long long int i = 0; i < pointDI.nbImageCells; i++)
	    {
	    numCell = pointDI.tabImageCells[i];
	    tempImageCell.cellNum = numCell;

	    if (numCell < nbCellsTotal)
		{
		tempImageCell.minVal = PLUS_INF;
		for (unsigned long long int j = pointDI.tabCellEntrees[i]; j < pointDI.tabCellEntrees[i + 1]; j++)
		    {
		    numControl = pointDI.tabImageControls[j];
		    tempL = dynsys->lFunc(doublePointCoords, controlCoords[numControl]);
		    (dynsys->*(dynsys->discretDynamics))(doublePointCoords, controlCoords[numControl], imageCoords, rho);
		    double tempL1 = dynsys->lFunc(imageCoords, controlCoords[numControl]);
		    tempImageCell.minVal = min(tempImageCell.minVal, (*itPoint).minVal + rho * 0.5 * (tempL + tempL1));
		    }
		addDataToCurrentImage(&itStart, tempImageCell, &itNew);
		itStart = itNew;
		addConvexCombinations(itPoint, numCell, &tempImageCell, rho, &itStart);
		}
	    }
	itPoint++;
	}
    }

void ViabiMicroMacro::addConvexCombinations(list<imagePoint>::iterator itPoint, unsigned long long int numCell, imageCell *tempImageCell, double rho,
	list<imageCell>::iterator *itStart)
    {

    /*!
     * On appelle d'abord deux fois la fonction numToIntAndDoubleCoords()
     * pour calculer les coordonnées réelles du point de départ x  et  du coin inférieur de la maille
     *  appartenant é l'image \f$ \Phi(x)\f$ y.
     */

    unsigned long long int posX = (*itPoint).PointNum;
    double pointVal = (*itPoint).minVal;
    double newCellVal = (*tempImageCell).minVal;

    double LVal = (newCellVal - pointVal) / rho;

    grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
    grid->numToIntAndDoubleCoords(numCell, intPointCoords, doubleVect);
    list<imageCell>::iterator itNew;
    double dist = 0.;
    /*!
     * Ensuite on calcule le vecteur différence \f$ z=y-x\f$ et sa norme \f$ \|z\|_2 \f$ .
     */
    for (int i = 0; i < dim; i++)
	{
	doubleVect[i] = doubleVect[i] - doublePointCoords[i];
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
		if (newCellNum != lastVisitCellNum)
		    {

		    (*tempImageCell).cellNum = newCellNum;
		    (*tempImageCell).minVal = pointVal + t * rho * LVal;
		    this->addDataToCurrentImage(itStart, (*tempImageCell), &itNew);
		    lastVisitCellNum = newCellNum;
		    (*itStart) = itNew;
		    }
		}
	    }
	t = t + deltat;
	}
    }
void ViabiMicroMacro::addDataToPointsList(list<imagePoint>::iterator *startIt, imagePoint newPoint, list<imagePoint>::iterator *resIt)
    {
    /*!
     * Dans l'image  en construction de \f$ C_{n}\setminus C_{n-1}\f$ les mailles sont stockées dans l'ordre croissant
     * de leur numéros  et chaque maille garde la mémoire de la valeur minimale ainsi que de tous les antécédents
     * de cette maille c'est é dire tous les couples viables (x,u) pour lesquels f(x,u) appartient é cette maille.
     * En plus  on distingue les couples (x,u) minimisant la valeur  et les autres, juste viables.
     *
     * Pour ajouter une nouvelle maille dans cette liste  on doit faire dans l'ordre les choses suivantes :
     *
     *    rechercher dans la liste la place nu numéro de maille é inserrer
     *   deux cas de figure peuvent se présenter :
     *     la maille ayant le méme numéro  existe déjé: on procéde alors é la fusion des deux,  en déterminant la valeur optimale et
     *    en  fusionnant les rétro-actions
     *     la maille n'existe pas encore; alors on l'inserre simplement  dans l'ordre croissant.
     *
     */
    unsigned long long int numnewPoint = newPoint.PointNum;
    list<imagePoint>::iterator itCell;

    if (currentImagePointsList.pointsList->size() == 0)
	{
	currentImagePointsList.pointsList->push_back(newPoint);
	currentImagePointsList.maxNum = numnewPoint;
	currentImagePointsList.minNum = numnewPoint;
	(*resIt) = currentImagePointsList.pointsList->end();
	(*resIt)--;
	}
    else
	{
	if (numnewPoint > currentImagePointsList.maxNum)
	    {
	    currentImagePointsList.pointsList->push_back(newPoint);
	    currentImagePointsList.maxNum = numnewPoint;
	    (*resIt) = currentImagePointsList.pointsList->end();
	    (*resIt)--;
	    }
	else
	    {
	    if (numnewPoint < currentImagePointsList.minNum)
		{
		currentImagePointsList.pointsList->push_front(newPoint);
		currentImagePointsList.minNum = numnewPoint;
		(*resIt) = currentImagePointsList.pointsList->begin();
		}
	    else
		{
		itCell = *startIt;

		if (numnewPoint < (*itCell).PointNum)
		    {

		    while ((numnewPoint < (*itCell).PointNum))
			{
			itCell--;
			}
		    if (numnewPoint > (*itCell).PointNum)
			{
			itCell++;
			currentImagePointsList.pointsList->insert(itCell, newPoint);
			}
		    else
			{
			(this->*addDataToCurrentPoint)(itCell, newPoint);
			}
		    (*resIt) = itCell;

		    }
		else
		    {
		    while ((numnewPoint > (*itCell).PointNum))
			{
			itCell++;
			}
		    if (numnewPoint < (*itCell).PointNum)
			{
			currentImagePointsList.pointsList->insert(itCell, newPoint);
			}
		    else
			{
			(this->*addDataToCurrentPoint)(itCell, newPoint);
			}
		    (*resIt) = itCell;

		    }
		}
	    }
	}
    }

void ViabiMicroMacro::addDataToGivenPointsList(imagePointsList *tempImagePointsList, list<imagePoint>::iterator *startIt, imagePoint newPoint,
	list<imagePoint>::iterator *resIt)
    {
    /*!
     * Dans l'image  en construction de \f$ C_{n}\setminus C_{n-1}\f$ les mailles sont stockées dans l'ordre croissant
     * de leur numéros  et chaque maille garde la mémoire de la valeur minimale ainsi que de tous les antécédents
     * de cette maille c'est é dire tous les couples viables (x,u) pour lesquels f(x,u) appartient é cette maille.
     * En plus  on distingue les couples (x,u) minimisant la valeur  et les autres, juste viables.
     *
     * Pour ajouter une nouvelle maille dans cette liste  on doit faire dans l'ordre les choses suivantes :
     *
     *    rechercher dans la liste la place nu numéro de maille é inserrer
     *   deux cas de figure peuvent se présenter :
     *     la maille ayant le méme numéro  existe déjé: on procéde alors é la fusion des deux,  en déterminant la valeur optimale et
     *    en  fusionnant les rétro-actions
     *     la maille n'existe pas encore; alors on l'inserre simplement  dans l'ordre croissant.
     *
     */
    unsigned long long int numnewPoint = newPoint.PointNum;
    list<imagePoint>::iterator itCell;

    if (tempImagePointsList->pointsList->size() == 0)
	{
	tempImagePointsList->pointsList->push_back(newPoint);
	tempImagePointsList->maxNum = numnewPoint;
	tempImagePointsList->minNum = numnewPoint;
	(*resIt) = tempImagePointsList->pointsList->end();
	(*resIt)--;
	}
    else
	{
	if (numnewPoint > tempImagePointsList->maxNum)
	    {
	    tempImagePointsList->pointsList->push_back(newPoint);
	    tempImagePointsList->maxNum = numnewPoint;
	    (*resIt) = tempImagePointsList->pointsList->end();
	    (*resIt)--;
	    }
	else
	    {
	    if (numnewPoint < tempImagePointsList->minNum)
		{
		tempImagePointsList->pointsList->push_front(newPoint);
		tempImagePointsList->minNum = numnewPoint;
		(*resIt) = tempImagePointsList->pointsList->begin();

		}
	    else
		{
		itCell = *startIt;

		if (numnewPoint < (*itCell).PointNum)
		    {

		    while ((numnewPoint < (*itCell).PointNum))
			{
			itCell--;
			}
		    if (numnewPoint > (*itCell).PointNum)
			{
			itCell++;
			tempImagePointsList->pointsList->insert(itCell, newPoint);
			}
		    else
			{
			(this->*addDataToCurrentPoint)(itCell, newPoint);
			}
		    (*resIt) = itCell;

		    }
		else
		    {
		    while ((numnewPoint > (*itCell).PointNum))
			{
			itCell++;
			}
		    if (numnewPoint < (*itCell).PointNum)
			{
			tempImagePointsList->pointsList->insert(itCell, newPoint);
			}
		    else
			{
			(this->*addDataToCurrentPoint)(itCell, newPoint);
			}
		    (*resIt) = itCell;

		    }
		}
	    }
	}
    }

void ViabiMicroMacro::addDataToPoint(list<imagePoint>::iterator itCell, imagePoint newPoint)
    {
    (*itCell).minVal = min((*itCell).minVal, newPoint.minVal);
    }


void ViabiMicroMacro::createPointsList()
    {

    cout << " create point list opti new\n";
    list<triple>::iterator itR;
    double t1, t2, elapsed_time;

    timeval tim;
    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    unsigned long long int posX;
    double testV[dim];
    unsigned long long int testI[dim];
    currentImagePointsList.pointsList->clear();
    currentImagePointsList.maxNum = -1;
    currentImagePointsList.minNum = currentImageList.minNum;

    int nbPointsCube = (int) pow(2.0, dim);	//pow(2.0, dim);
    list<imageCell>::iterator itCell = currentImageList.cellsList.begin(), itLast = currentImageList.cellsList.end();
    int i;

    long long int *indicesDecalCell = grid->getIndicesDecalCell();

    imagePoint tempPoint;
    list<imagePoint>::iterator itStart, itNew;
    itStart = this->currentImagePointsList.pointsList->begin();
    for (i = 0; i < nbPointsCube - 1; i++)
	{
	while (itCell != itLast)
	    {
	    tempPoint.minVal = (*itCell).minVal;

	    posX = (*itCell).cellNum + indicesDecalCell[i];
	    grid->numToIntAndDoubleCoords(posX, testI, testV);
	    if (dynsys->constraintsX(testV) < PLUS_INF)
		{

		tempPoint.PointNum = (*itCell).cellNum + indicesDecalCell[i];

		this->addDataToPointsList(&itStart, tempPoint, &itNew);
		itStart = itNew;

		}
	    itCell++;

	    }
	itCell = currentImageList.cellsList.begin();
	}

    while (!currentImageList.cellsList.empty())	//(itCell!=itLast)
	{
	tempPoint.minVal = (*itCell).minVal;

	/*!
	 * \todo procéder comme pour la création de la liste de cellules:
	 * créer un point avec le numéro et les données de la cellule.
	 * copier une seule fois les données de la cellules sur la valeur et la rétro-action
	 *  d'un point é l'autre de la méme cellule  seul le numéro du point change.
	 *
	 *  Puis inserrer le point dans la liste ordonnée de points.
	 *   Exactement comme pour les celllules! Donc du copier coller de code.
	 *
	 */
	tempPoint.PointNum = (*itCell).cellNum + indicesDecalCell[nbPointsCube - 1];
	this->addDataToPointsList(&itStart, tempPoint, &itNew);
	itStart = itNew;

	itCell++;
	currentImageList.cellsList.pop_front();
	}

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));

    cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;

    }

void ViabiMicroMacro::addDataToCurrentImage(list<imageCell>::iterator *startIt, imageCell newCell, list<imageCell>::iterator *resIt)
    {
    /*!
     * Dans l'image  en construction de \f$ C_{n}\setminus C_{n-1}\f$ les mailles sont stockées dans l'ordre croissant
     * de leur numéros  et chaque maille garde la mémoire de la valeur minimale ainsi que de tous les antécédents
     * de cette maille c'est é dire tous les couples viables (x,u) pour lesquels f(x,u) appartient é cette maille.
     * En plus  on distingue les couples (x,u) minimisant la valeur  et les autres, juste viables.
     *
     * Pour ajouter une nouvelle maille dans cette liste  on doit faire dans l'ordre les choses suivantes :
     *
     *    rechercher dans la liste la place nu numéro de maille é inserer
     *   deux cas de figure peuvent se présenter :
     *     la maille ayant le méme numéro  existe déjé: on procéde alors é la fusion des deux,  en déterminant la valeur optimale et
     *    en  fusionnant les rétro-actions
     *     la maille n'existe pas encore; alors on l'inserre simplement  dans l'ordre croissant.
     *
     */

    unsigned long long int numNewCell = newCell.cellNum;

    list<imageCell>::iterator itCell, itLast = currentImageList.cellsList.end();

    if (currentImageList.cellsList.size() == 0)
	{
	//cout<< " c1"<<endl;
	currentImageList.cellsList.push_back(newCell);
	currentImageList.maxNum = numNewCell;
	currentImageList.minNum = numNewCell;
	(*resIt) = currentImageList.cellsList.end();
	(*resIt)--;
	}
    else
	{

	if ((int) numNewCell > currentImageList.maxNum)
	    {
	    currentImageList.cellsList.push_back(newCell);
	    currentImageList.maxNum = numNewCell;
	    currentImageList.minNum = min(currentImageList.maxNum, currentImageList.minNum);
	    (*resIt) = currentImageList.cellsList.end();
	    (*resIt)--;
	    }
	else
	    {
	    if ((int) numNewCell < currentImageList.minNum)
		{
		currentImageList.cellsList.push_front(newCell);
		currentImageList.minNum = numNewCell;
		currentImageList.maxNum = max(currentImageList.maxNum, currentImageList.minNum);
		(*resIt) = currentImageList.cellsList.begin();
		}
	    else
		{
		itCell = *startIt;

		if ((numNewCell > (*itCell).cellNum))
		    {
		    while ((itCell != itLast) && (numNewCell > (*itCell).cellNum))
			{
			itCell++;
			}
		    if (numNewCell < (*itCell).cellNum)
			{
			currentImageList.cellsList.insert(itCell, newCell);
			}
		    else
			{
			(this->*addDataToCurrentCell)(itCell, newCell);
			}
		    (*resIt) = itCell;
		    }
		else
		    {
		    while ((numNewCell < (*itCell).cellNum))
			{
			itCell--;
			}
		    if (numNewCell > (*itCell).cellNum)
			{
			itCell++;
			currentImageList.cellsList.insert(itCell, newCell);
			}
		    else
			{
			(this->*addDataToCurrentCell)(itCell, newCell);
			}
		    (*resIt) = itCell;
		    }
		}
	    }
	}
    }

void ViabiMicroMacro::addDataToCell(list<imageCell>::iterator itCell, imageCell newCell)
    {

    (*itCell).minVal = min((*itCell).minVal, newCell.minVal);

    }

/*!
 * Le destructeur nettoie la mémoire réservée pour certaines variables globales é la classe et ferme
 * les bases de données de rétro-actions
 */
ViabiMicroMacro::~ViabiMicroMacro()
    {

    cout << " desturction de classe viabi hjb \n";
    // TODO Auto-generated destructor stub
    delete[] pointDI.tabCellEntrees;
    delete[] pointDI.tabImageCells;
    delete[] pointDI.tabImageControls;

    delete[] intPointCoords;
    delete[] doublePointCoords;
    delete[] doubleVect;
    delete[] doubleVect1;
    delete[] intControlCoords;
    delete[] doubleControlCoords;
    delete[] imageCells;

    cout << " divers tableaux OK\n";

    }

void ViabiMicroMacro::printViabiInfo()
    {
    grid->printGrid();
    }

void ViabiMicroMacro::showCurrentImageList()
    {
    imageCell c;
    list<imageCell>::iterator itCell = currentImageList.cellsList.begin(), itLast = currentImageList.cellsList.end();
    list<triple>::iterator it;
    cout << " image  de CN calculée  est la suivante \n";
    cout << "*******************************************************\n";
    while ((itCell != itLast))
	{
	c = (*itCell);

	cout << " maille  num " << c.cellNum;
	cout << " value= " << c.minVal << endl;

	itCell++;
	}
    cout << "*******************************************************\n";
    }
void ViabiMicroMacro::showCurrentImagePointsList()
    {
    imagePoint c;
    list<imagePoint>::iterator itCell = currentImagePointsList.pointsList->begin(), itLast = currentImagePointsList.pointsList->end();
    cout << " image  de CN calculée  est la suivante : liste de POINTS \n";
    cout << "*******************************************************\n";
    while ((itCell != itLast))
	{
	c = (*itCell);

	cout << " point  num " << c.PointNum;
	cout << " value= " << c.minVal;
	itCell++;
	}
    cout << "*******************************************************\n";
    }

int ViabiMicroMacro::addNewPoints()
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

    list<imagePoint>::iterator itPoint = currentImagePointsList.pointsList->begin(), itLastPoint = currentImagePointsList.pointsList->end(), itTemp;

    imagePoint tempPoint;
    while (itPoint != itLastPoint)
	{

	tempPoint = (*itPoint);
	if ((*itPoint).minVal < vTab[(*itPoint).PointNum])
	    {
	    nbNewPoints++;
	    vTab[(*itPoint).PointNum] = min(vTab[(*itPoint).PointNum], (*itPoint).minVal);
	    itPoint++;
	    }
	else
	    {
	    itTemp = itPoint;
	    itPoint++;
	    currentImagePointsList.pointsList->erase(itTemp);
	    }
	}

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));

    cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;
    return nbNewPoints;
    }

bool ViabiMicroMacro::testConstraintesForCell(unsigned long long int numCell)
    {
    bool res = true;
    long long int *indicesDecalCell = grid->getIndicesDecalCell();
    double *doublePointCoords = new double[dim];

    int nbPointsCube = (int) pow(2.0, dim);	//pow(2.0, dim);
    unsigned long long int posX;
    int i = 0;

    while (res & (i < nbPointsCube))
	{
	posX = numCell + indicesDecalCell[i];
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	res = (dynsys->constraintsX(doublePointCoords) < PLUS_INF);

	i++;

	}

    delete[] doublePointCoords;

    return res;
    }

void ViabiMicroMacro::initialiseConstraints()
    {
    spdlog::info("Initialization on the contrats set K0");


	initialiseConstraints_CC();

    }

void ViabiMicroMacro::initialiseConstraints_CC()
    {
    spdlog::info("Initialization on the contrats set K0");

    unsigned long long int iCoords[dim];
    double xCoords[dim];
    unsigned long long int totalPointsX = grid->getNbTotalPoints();
    for (unsigned long long int pos = 0; pos < totalPointsX; pos++)
	{
	grid->numToIntAndDoubleCoords(pos, iCoords, xCoords);

	vTab[pos] = dynsys->constraintsX(xCoords);
	}
    }


void ViabiMicroMacro::ViabilityKernel(bool sortieOK, int nbArret)
    {
    dynsys->setDynamicsForward();

	if (nbOMPThreads > 1)
	    {
	    viabKerValFunc_omp(nbArret);
	    }
	else
	    {
	    viabKerValFunc(nbArret);
	    }


    }

void ViabiMicroMacro::GarantedViabilityKernel(bool sortieOK, int nbArret)
    {
    dynsys->setDynamicsForward();

    }

