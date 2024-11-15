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

    currentPointsImage = map<unsigned long long int, double>();
    currentCellsImage = map<unsigned long long int, double>();

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

    ViabiMicroMacro::addNewPointsToSet = &ViabiMicroMacro::addNewPoints_new;

    ViabiMicroMacro::createCurrentPointsList = &ViabiMicroMacro::createPointsList_new;

    if (avp.COMPUTE_TMIN)
	{
	ViabiMicroMacro::computeFirstConvexifiedImage = &ViabiMicroMacro::computeConvexifiedImage_tmin_new;
	ViabiMicroMacro::computeCurrentImage = &ViabiMicroMacro::computeCurrIm_tmin_new;
	//  ViabiMicroMacro::computeFirstConvexifiedImage_omp=&ViabiMicroMacro::computeConvexifiedImage_tmin_omp;
	}
    else
	{

	ViabiMicroMacro::computeCurrentImage = &ViabiMicroMacro::computeCurrIm_Lmin_new;
	ViabiMicroMacro::computeFirstConvexifiedImage = &ViabiMicroMacro::computeConvexifiedImage_Lmin_new;

	}


    whichPointListToUse = 1;
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

    initialiseTargetHJB_new();
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

void ViabiMicroMacro::initialiseTargetHJB_new()
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
    currentPointsImage.clear();
    currentCellsImage.clear();

    unsigned long long int pos;

    map<unsigned long long int, double>::iterator itStart = currentPointsImage.begin(), itNew;
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

	    currentPointsImage[pos] = c;
	    }
	}
    spdlog::info("Target set initialized, number of grid points in target set : {0:d}", totalPointsC);
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

void ViabiMicroMacro::viabKerValFunc_new(unsigned long long int nbArret)
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


    DiscretPointImage *tempPointImage = new DiscretPointImage((unsigned long long int) 1, dim, dynsys, grid);
    map<unsigned long long int, list<unsigned long long int> > *cellsWithControls = tempPointImage->GetImageCellsWithControls();
    map<unsigned long long int, list<unsigned long long int> >::iterator itCell, itLastCell, itCellsTemp;
    list<unsigned long long int> retroControls;
    list<unsigned long long int>::iterator ic, icLast;
    double imageCoords[dim];
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
		catch (...)
		    {
		    std::cout << "Allocation failed dans calculRho: ";
		    cout.flush();
		    throw;
		    }

		bool testNonVide = false;

		tempPointImage->setPointNumAndRebuild(pos);

		minValCell = PLUS_INF;
		itCell = cellsWithControls->begin();
		itLastCell = cellsWithControls->end();
		while (itCell != itLastCell && !testNonVide)
		    {
		    cellNum = itCell->first;

		    if (cellNum < nbCellsTotal)
			{
			retroControls = itCell->second;
			for (int iCell = 0; iCell < nbPointsCube - 1; iCell++)
			    {
			    double imageVal = vTab[cellNum + indicesDecalCell[iCell]];

			    ic = retroControls.begin();
			    icLast = retroControls.end();
			    tempVal = PLUS_INF;
			    while(ic != icLast && !testNonVide)
				{
				numControl = (*ic);
				tempL = dynsys->lFunc(rCoords, controlCoords[numControl]);

				tempM = dynsys->mFunc(rCoords, controlCoords[numControl]);
				if (tempL < PLUS_INF && tempM < PLUS_INF)
				    {
				    tempVal = (imageVal + rho * tempL) / (1 - rho * tempM);
				    testNonVide = tempVal < valAtPos;
				    minValCell = min(minValCell, tempVal);
				    }
				ic++;
				}

			    }
			}
		    else if (cellNum == nbCellsTotal)
			{
			retroControls = itCell->second;
			for (int iCell = 0; iCell < nbPointsCube - 1; iCell++)
			    {
			    double imageVal = vTab[cellNum + indicesDecalCell[iCell]];

			    ic = retroControls.begin();
			    icLast = retroControls.end();
			    tempVal = PLUS_INF;
			    while(ic != icLast && !testNonVide)
				{
				numControl = (*ic);
				(dynsys->*(dynsys->discretDynamics))(rCoords, controlCoords[numControl], imageCoords, rho);
				printVector(doubleVect1, dim);
				tempL = dynsys->lFunc(rCoords, controlCoords[numControl]);
				tempM = dynsys->mFunc(rCoords, controlCoords[numControl]);

				double tempV = dynsys->constraintsX(imageCoords);
				tempVal = (tempV + rho * tempL) / (1 - rho * tempM);
				testNonVide = tempVal < valAtPos;
				minValCell = min(minValCell, tempVal);
				ic++;
				}

			    }
			}

		    itCell++;
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
    try
    {
	spdlog::info("STEP 1: compute first convex image ");
	(this->*computeFirstConvexifiedImage)(iter);    //adds convex combinations
    }
    catch (...)
	{
	cout << " exception while computing first convexified image " << std::flush;
	cout.flush();
	throw;
	}

    /*!
     * -A partir de la liste des mailles on crée une liste ordonnée de points représentant l'image é l'aide de la fonction
     * createPointsList().
     */

    try
    {
	spdlog::info("STEP 2: Gets current points list");
	(this->*createCurrentPointsList)();    //ok
    }
    catch (...)
	{
	cout << " exception while computing currentpoints list " << std::flush;
	cout.flush();
	throw;
	}
    /*!
     * - A partir de la liste des points représentant l'image  trois bases  de données sont alimentées: La base
     * principale, associée é la grille  et représentant la fonction valeur, la base de rétroaction optimale et la base de
     *  rétro-action viable. On appelle pour cela la fonction addNewPoints().
     */
    if (!(dynsys->getDynType() == DD))
	{
	try
	{
	    spdlog::info("STEP 3: Update of value function ");
	    cout << " add new points to set " << std::flush;
	    nbNewPoints = (this->*addNewPointsToSet)();
	}
	catch (...)
	    {
	    cout << " exception while computing addPointsToSet " << std::flush;
	    cout.flush();
	    throw;
	    }
	}

    iter++;

    /*!
     * Tant qu'il y a de nouveaux points ajoutés on répéte les opérations suivantes.
     */
    while ((nbNewPoints > 0))
	{
	//spdlog::info("Iteration {} : new points in the set {}", iter,
	//		nbNewPoints);

	/*!
	 * -Calcul de l'image de \f^C_{n}\setminus C_n\f$  par la fonction computeCurrentImageLocalRho(): l'image est enregistrée en émoire vive sous forme de liste
	 * de références de mailles dans lesquelles arrive au moins une évolution. Chaque référence de maille contient
	 * des informations sur tous les antécédants de cette maille ainsi  que la valeur minimale de
	 * temps. Dans cette version oé \f$ \rho\f$ est global, la fonction valeur prend la méme valeur
	 * é chaque étape : \f$ \rho \cdot n \f$.
	 */
	try
	{
	    cout << " compute current image " << std::flush;
	    (this->*computeCurrentImage)(iter);
	}
	catch (...)
	    {
	    cout << " exception while computing currentImage " << std::flush;
	    cout.flush();
	    throw;
	    }
	//this->showCurrentImageList();

	/*!
	 * -A partir de la liste des mailles on crée une liste ordonnée de points représentant l'image en appelant la fonction createPointsList().
	 *  Comme pour le mailles, chaque référence de point regroupe les informatons (regroupées é partir de différentes
	 * mailles dont est vertex)  sur les antécédents  de ce point. Le but de la création de cette liste est d'éliminer
	 * les doublons afin de minimiser les accés é la base de données
	 */
	try
	{
	    cout << " create point list " << std::flush;
	    (this->*createCurrentPointsList)();
	}
	catch (...)
	    {
	    cout << " exception while computing CreateCurrentPointsList " << std::flush;
	    cout.flush();
	    throw;
	    }
	//this->showCurrentImagePointsList();
	/*!
	 * - A partir de la liste des points représentant l'image  trois bases  de données sont alimentées: La base
	 * principale, associée é la grille  et représentant la fonction valeur, la base de rétroaction optimale et la base de
	 *  rétro-action viable. On appelle ici la fonction addNewPoints().
	 */
	try
	{
	    cout << " add new points " << std::flush;
	    nbNewPoints = (this->*addNewPointsToSet)();
	}
	catch (...)
	    {
	    cout << " exception while computing addPointsToSet " << std::flush;
	    cout.flush();
	    throw;
	    }
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


void ViabiMicroMacro::computeCurrIm_tmin_new(int iter)
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


void ViabiMicroMacro::computeCurrIm_Lmin_new(int iter)
    {
    spdlog::info("[Integral cost problem] : Computing current image. number of points {}", currentCellsImage.size());
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
    double imageCoords[dim];
    double tempVal, pointVal;
    unsigned long long int numControl;
    double tempL, tempL1;
    double **controlCoords = dynsys->getControlCoords();
    unsigned long long int cellNum;
    list<unsigned long long int> retroControls;
    list<unsigned long long int>::iterator ic, icLast;

    while (itPoint != itLastPoint)
	{
	posX = itPoint->first;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	rho = dynsys->calculRho_local(doublePointCoords);
	/*!
	 * On calcule l'image discrète  du point
	 */

	tempPointImage->setPointNumAndRebuild(posX);
	pointVal = itPoint->second;
	itCell = cellsWithControls->begin();
	itLastCell = cellsWithControls->end();
	while (itCell != itLastCell)
	    {

	    cellNum = itCell->first;
	    if(cellNum < grid->getNbTotalPoints())
		{
		retroControls = itCell->second;
		ic = retroControls.begin();
		icLast = retroControls.end();
		tempVal = PLUS_INF;
		while(ic != icLast)
		    {
		    numControl = (*ic);
		    tempL = dynsys->lFunc(doublePointCoords, controlCoords[numControl]);
		    (dynsys->*(dynsys->discretDynamics))(doublePointCoords, controlCoords[numControl], imageCoords, 1.0);

		    double tempL1 = dynsys->lFunc(imageCoords, controlCoords[numControl]);
		    tempVal = min(tempVal, pointVal + rho * 0.5 * (tempL + tempL1));
		    ic++;
		    }

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
    currentPointsImage.clear();
    spdlog::info("[Integral cost problem] : Computing current image finished. number of points {}", currentCellsImage.size());

    }

void ViabiMicroMacro::computeConvexifiedImage_tmin_new(int iter)
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
		addConvexCombinations_new(posX, itPoint->second, tempVal, itCell->first, rho);
		}
	    itCell++;
	    }
	itPoint++;
	}
    spdlog::info("[Min Time problem] : Computing convexified image finished. number of points {}", currentPointsImage.size());
    }


void ViabiMicroMacro::computeConvexifiedImage_Lmin_new(int iter)
    {
    spdlog::info("[Integral cost problem] : Computing current image. number of points {}", currentCellsImage.size());
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
    double imageCoords[dim];
    double tempVal, pointVal;
    unsigned long long int numControl;
    double tempL, tempL1;
    double **controlCoords = dynsys->getControlCoords();
    unsigned long long int cellNum;
    list<unsigned long long int> retroControls;
    list<unsigned long long int>::iterator ic, icLast;

    while (itPoint != itLastPoint)
	{
	posX = itPoint->first;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	rho = dynsys->calculRho_local(doublePointCoords);
	/*!
	 * On calcule l'image discrète  du point
	 */

	tempPointImage->setPointNumAndRebuild(posX);
	pointVal = itPoint->second;
	itCell = cellsWithControls->begin();
	itLastCell = cellsWithControls->end();
	while (itCell != itLastCell)
	    {

	    cellNum = itCell->first;
	    if(cellNum < grid->getNbTotalPoints())
		{
		retroControls = itCell->second;
		ic = retroControls.begin();
		icLast = retroControls.end();
		tempVal = PLUS_INF;
		while(ic != icLast)
		    {
		    numControl = (*ic);
		    tempL = dynsys->lFunc(doublePointCoords, controlCoords[numControl]);
		    (dynsys->*(dynsys->discretDynamics))(doublePointCoords, controlCoords[numControl], imageCoords, 1.0);

		    double tempL1 = dynsys->lFunc(imageCoords, controlCoords[numControl]);
		    tempVal = min(tempVal, pointVal + rho * 0.5 * (tempL + tempL1));
		    ic++;
		    }

		if (auto result = currentCellsImage.find(cellNum); result != currentCellsImage.end())
		    {
		    result->second = min((double) result->second, tempVal);
		    }
		else
		    {
		    currentCellsImage[cellNum] = tempVal;
		    }
		addConvexCombinations_new(posX, itPoint->second, tempVal, itCell->first, rho);
		}
	    itCell++;
	    }

	itPoint++;
	}
    currentPointsImage.clear();
    spdlog::info("[Integral cost problem] : Computing current image finished. number of points {}", currentCellsImage.size());

    }


void ViabiMicroMacro::addConvexCombinations_new(unsigned long long int posX, double pointVal, double newCellVal, unsigned long long int numCell,
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


void ViabiMicroMacro::createPointsList_new()
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

    long long int *indicesDecalCell = grid->getIndicesDecalCell();

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

int ViabiMicroMacro::addNewPoints_new()
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

    map<unsigned long long int, double>::iterator itPoint = currentPointsImage.begin(), itLastPoint = currentPointsImage.end(), itTemp;

    imagePoint tempPoint;
    while (itPoint != itLastPoint)
	{
	unsigned long long int pointNum = itPoint->first;
	double value = itPoint->second;

	if (value < vTab[pointNum])
	    {
	    nbNewPoints++;
	    vTab[pointNum] = value;
	    itPoint++;
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
    viabKerValFunc_new(nbArret);
    /*
    if (nbOMPThreads > 1)
	{
	viabKerValFunc_omp(nbArret);
	}
    else
	{
	viabKerValFunc_new(nbArret);
	}
     */
    }

void ViabiMicroMacro::GarantedViabilityKernel(bool sortieOK, int nbArret)
    {
    dynsys->setDynamicsForward();

    }

