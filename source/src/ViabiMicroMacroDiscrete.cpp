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

#include "../include/ViabiMicroMacroDiscrete.h"

ViabiMicroMacroDiscrete::ViabiMicroMacroDiscrete(ParametersManager *pm)
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
    InitViabiMicroMacroDiscrete(avp);
    trajectoryHelper = new ViabiMicroMacroTrajectoryHelper(grid, dynsys, avp.TYPE_TRAJ);
    }

SysDyn* ViabiMicroMacroDiscrete::GetSysDynForViabProblem()
    {
    return this->dynsys;
    }

GridMicroMacro* ViabiMicroMacroDiscrete::GetGridForViabProblem()
    {
    return grid;
    }

void ViabiMicroMacroDiscrete::InitViabiMicroMacroDiscrete(algoViabiParams avp)
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
    pointDI.tabPointsEntrees = new unsigned long long int[dynsys->getTotalNbPointsC() + 1];
    pointDI.tabImagePoints = new unsigned long long int[dynsys->getTotalNbPointsC()];
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

    ViabiMicroMacroDiscrete::addNewPointsToSet = &ViabiMicroMacroDiscrete::addNewPoints;
    ViabiMicroMacroDiscrete::addDataToCurrentCell = &ViabiMicroMacroDiscrete::addDataToCell;
    ViabiMicroMacroDiscrete::addDataToCurrentPoint = &ViabiMicroMacroDiscrete::addDataToPoint;

    ViabiMicroMacroDiscrete::createCurrentPointsList = &ViabiMicroMacroDiscrete::createPointsList;

    ViabiMicroMacroDiscrete::computeCurrentImage = &ViabiMicroMacroDiscrete::computeCurrIm;
    tempPointsList1 = new list<imagePoint>();
    tempPointsList2 = new list<imagePoint>();
    ;
    whichPointListToUse = 1;
    currentImagePointsList.pointsList = tempPointsList1;
    }

void ViabiMicroMacroDiscrete::computeDiscreteImageOfPoint(unsigned long long int num)
    {

    unsigned long long int **controlCoords = dynsys->getControlIntCoords();
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int nbPointsTotal = grid->getNbTotalPoints();

    unsigned long long int imagePos;

    grid->numToIntAndDoubleCoords(num, intPointCoords, doublePointCoords);
    unsigned long long int cu;

    list<uintPair> pointsList;
    for (cu = 0; cu < nbCTotal; cu++)
	{
	/*!
	 * on calcule les coordonnées réelles de l'image du point par la dynamique discrete
	 * elles sont stockes dans le tableau doubleVect
	 * si \f$ u\in U(x)\f$  (contréle admissible) on calcule l'image réelle par la dynamique
	 *  discrétisée  en temps  du point x avec le controle u
	 */
	if (dynsys->constraintsXU_fd(intPointCoords, controlCoords[cu]) < PLUS_INF)
	    {

	    dynsys->dynamics_fd(intPointCoords, controlCoords[cu], intVect1);
	    if (grid->isPointInGrid_fd(intVect1))
		{
		grid->intCoordsToNum(intVect1, &imagePos);
		/*!
		 * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
		 */
		if (vTab[imagePos] < PLUS_INF)
		    {
		    /*!
		     * Si l'image  est dans l'ensemble de contraintes sur l'état \f$ K \f$
		     * on calcule le numéro de maille qui contient cett image
		     */
		    imageCells[cu] = imagePos;
		    // on enregistre le numero de maille
		    pointsList.push_back(uintPair(imageCells[cu], cu));
		    }
		else
		    {
		    /*!
		     * Si l'image n'est pas dans  \f$ K \f$ on enregistre un numéro de maille factice qui signifie que
		     * cette image est rejetée
		     */
		    imageCells[cu] = nbPointsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble
		    try
			{
			pointsList.push_back(uintPair(imageCells[cu], cu));
			}
		    catch (const std::bad_alloc &e)
			{
			std::cout << "Allocation failed: in discretImage points push_back 1" << e.what() << '\n';
			}
		    }
		}
	    else
		{
		/*!
		 * \todo prévoir la gestion des autorisations de sortie  par axe pour les variables qui peuvent
		 * se trouver dans des intervalles non bornés
		 */
		imageCells[cu] = nbPointsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble
		try
		    {
		    pointsList.push_back(uintPair(imageCells[cu], cu));
		    }
		catch (const std::bad_alloc &e)
		    {
		    std::cout << "Allocation failed: in discretImage points push_back 2" << e.what() << '\n';
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
	    imageCells[cu] = nbPointsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble
	    try
		{
		pointsList.push_back(uintPair(imageCells[cu], cu));
		}
	    catch (const std::bad_alloc &e)
		{
		std::cout << "Allocation failed: in discretImage points push_back 3" << e.what() << '\n';
		}
	    }

	}
    /*!
     * Toutes les images sont calculées  et stockées dans un tableau dans l'ordre
     *  de numérotation des controles. La tache suivante consiste é trier ce tableau, éliminer les doublons
     *   et  compléter la structure  image de point
     */

    pointsList.sort(upairCompare);
    list<uintPair>::iterator itCell, itDouble;

    itCell = pointsList.begin();
    itDouble = itCell;
    unsigned long long int currentPoint;

    cu = 0;
    unsigned long long int iControlTab = 0, iCellsTab = 0, iEntreesTab = 0;
    while ((itCell != pointsList.end()) & ((*itCell).first < nbPointsTotal))
	{

	currentPoint = (*itCell).first;

	pointDI.tabPointsEntrees[iEntreesTab] = iControlTab;
	pointDI.tabImagePoints[iCellsTab] = currentPoint;
	iEntreesTab++;
	iCellsTab++;

	while ((itDouble != pointsList.end()) & ((*itDouble).first == currentPoint))
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

    pointDI.nbImagePoints = iCellsTab;
    /*!
     * La derniére valeur  du tableau des controles indique la fin  de la liste des controles associés
     * é la derniére cellule
     */
    pointDI.tabPointsEntrees[iCellsTab] = iControlTab;
    }

void ViabiMicroMacroDiscrete::initialiseTarget()
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

void ViabiMicroMacroDiscrete::computeTrajectories()
    {

    dynsys->setDynamicsForward();

    algoViabiParams *avp = modelParams->getAlgoParameters();
    if (avp->TYPE_TRAJ == OP)
	computeOptimalTrajectories();
    else
	computeViableTrajectories();

    }

void ViabiMicroMacroDiscrete::computeOptimalTrajectories()
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

void ViabiMicroMacroDiscrete::computeViableTrajectories()
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

	    logVector("Initial point ", avp->INIT_POINTS_FD + tr * dim, dim);
	    if (typeTraj == VG)
		{
		trajectoryHelper->computeViableTrajectory_tych_DD(avp->INIT_POINTS_FD + tr * dim, avp->INIT_VALUES_FD[tr], fileName, success[tr]);
		}
	    else
		{
		trajectoryHelper->computeViableTrajectory_DD(avp->INIT_POINTS_FD + tr * dim, avp->INIT_VALUES_FD[tr], fileName, success[tr]);
		}
	    }

	}
    }

void ViabiMicroMacroDiscrete::initialiseTargetHJB()
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

	c = max((*(dynsys->target_fd))(x), (*(dynsys->constraintsX_fd))(x));

	vTab[pos] = c;
	if (c < PLUS_INF)
	    {
	    totalPointsC++;
	    currentPoint.minVal = c;
	    currentPoint.PointNum = pos;

	    addDataToPointsList(&itStart, currentPoint, &itNew);
	    //cout<<" fini ajout point  \n";
	    // printVector(x,dim);
	    //////system("pause");
	    itStart = itNew;
	    cpt++;
	    }
	}
    spdlog::info("Target set initialized, number of grid points in target set : {0:d}", cpt);
    }

void ViabiMicroMacroDiscrete::loadViableSets()
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

void ViabiMicroMacroDiscrete::saveViableSets()
    {
    saveValFunctions();
    }
void ViabiMicroMacroDiscrete::saveValFunctions()
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

void ViabiMicroMacroDiscrete::viabKerValFunc(unsigned long long int nbArret)
    {
    unsigned long long int cptChanged = 1000000, nbIter = 0;
    unsigned long long int pos;
    unsigned long long int nbPointsTotal = grid->getNbTotalPoints();
    spdlog::info("Viability kernel algorithm for discrete systems");
    unsigned long long int **controlCoords = dynsys->getControlIntCoords();

    uintPair image;

    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int imagePos;
    double minValCell;

    while ((cptChanged > nbArret) & (nbIter < 15000))
	{
	cptChanged = 0;
	/*
	 *  on parcourt  tous les points de l'espace discret  fini
	 *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
	 */
	for (pos = 0; pos < (unsigned long long int) nbPointsTotal; pos++)
	    {
	    double currentVal = vTab[pos];
	    grid->numToIntAndDoubleCoords(pos, intPointCoords, doublePointCoords);
	    if (currentVal < PLUS_INF)
		{
		double tempVal, tempL, tempM;
		minValCell = PLUS_INF;
		unsigned long long int cu = 0;
		bool testNonVide = false;
		while (!testNonVide && (cu < nbCTotal))
		    {
		    if (dynsys->constraintsXU_fd(intPointCoords, controlCoords[cu]) < PLUS_INF)
			{
			dynsys->dynamics_fd(intPointCoords, controlCoords[cu], intVect1);
			if (grid->isPointInGrid_fd(intVect1))
			    {
			    grid->intCoordsToNum(intVect1, &imagePos);
			    /*!
			     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
			     */
			    if (vTab[imagePos] < PLUS_INF)
				{
				tempL = dynsys->lFunc_fd(intPointCoords, controlCoords[cu]);

				tempM = dynsys->muFunc_fd(intPointCoords, controlCoords[cu]);

				tempVal = max((vTab[imagePos] - tempL), tempM);

				minValCell = min(minValCell, tempVal);
				testNonVide = (tempVal <= currentVal);
				}
			    }
			}
		    cu++;
		    }
		}
	    if (currentVal < minValCell)
		{
		cptChanged++;
		}
	    vTab[pos] = max(vTab[pos], minValCell);
	    }

	nbIter++;
	spdlog::info("Iteration {} finished. Number of grid points with updated value funtion : {}", nbIter, cptChanged);
	}

    saveValFunctions();
    }

void ViabiMicroMacroDiscrete::viabKerGarantiValFunc(unsigned long long int nbArret)
    {
    unsigned long long int cptChanged = 1000000, nbIter = 0;
    unsigned long long int pos;
    unsigned long long int nbPointsTotal = grid->getNbTotalPoints();
    spdlog::info("Guaranted viability kernel algorithm for discrete systems with tychastic controls");
    unsigned long long int **controlCoords = dynsys->getControlIntCoords();
    unsigned long long int **tychIntCoords = dynsys->getTychIntCoords();

    int dimTy = dynsys->getDimTy();

    uintPair image;

    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int imagePos;

    double minValCell;

    while ((cptChanged > nbArret) & (nbIter < 15000))
	{
	cptChanged = 0;
	/*
	 *  on parcourt  tous les points de l'espace discret  fini
	 *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
	 */
	for (pos = 0; pos < (unsigned long long int) nbPointsTotal; pos++)
	    {
	    double currentVal = vTab[pos];
	    grid->numToIntAndDoubleCoords(pos, intPointCoords, doublePointCoords);
	    if (currentVal < PLUS_INF)
		{
		double tempVal, tempL, tempM;
		minValCell = PLUS_INF;
		unsigned long long int cu = 0;
		bool testNonVide = false;
		while (!testNonVide && (cu < nbCTotal))
		    {
		    if (dynsys->constraintsXU_fd(intPointCoords, controlCoords[cu]) < PLUS_INF)
			{
			dynsys->dynamics_tych_fd(intPointCoords, controlCoords[cu], tychIntCoords[0], intVect1);
			if (grid->isPointInGrid_fd(intVect1))
			    {
			    grid->intCoordsToNum(intVect1, &imagePos);
			    /*!
			     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
			     */
			    if (vTab[imagePos] < PLUS_INF)
				{
				double tempValTy = -PLUS_INF;
				tempM = dynsys->muFunc_fd(intPointCoords, controlCoords[cu]);
				for (int jj = 0; jj < dimTy; jj++)
				    {
				    tempL = dynsys->lFunc_tych_fd(intPointCoords, controlCoords[cu], tychIntCoords[jj]);

				    tempValTy = max(tempValTy, max((vTab[imagePos] - tempL), tempM));
				    }

				tempVal = tempValTy;
				minValCell = min(minValCell, tempVal);
				testNonVide = (tempVal <= currentVal);
				}
			    }
			}
		    cu++;
		    }
		}
	    if (currentVal < minValCell)
		{
		cptChanged++;
		}
	    vTab[pos] = max(vTab[pos], minValCell);
	    }

	nbIter++;
	spdlog::info("Iteration {} finished. Number of grid points with updated value funtion : {}", nbIter, cptChanged);
	}

    saveValFunctions();
    }

void ViabiMicroMacroDiscrete::CaptureBasin()
    {
    spdlog::info("Capture bassin algorithm : start");

    dynsys->setDynamicsBackward();
    /*!
     * \var nbNewPoints : nb de nouveaux points ajoutés é l'étape n
     */
    int nbNewPoints = 1;

    int iter = 0;

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

    nbNewPoints = (this->*addNewPointsToSet)();

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

	/*!
	 * -A partir de la liste des mailles on crée une liste ordonnée de points représentant l'image en appelant la fonction createPointsList().
	 *  Comme pour le mailles, chaque référence de point regroupe les informatons (regroupées é partir de différentes
	 * mailles dont est vertex)  sur les antécédents  de ce point. Le but de la création de cette liste est d'éliminer
	 * les doublons afin de minimiser les accés é la base de données
	 */
	(this->*createCurrentPointsList)();
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

double ViabiMicroMacroDiscrete::computeOptimalCaptTrajectory(double *initPosition, string fileName, bool &succes)
    {
    return trajectoryHelper->computeOptimalTrajectory(initPosition, fileName, succes);
    }

void ViabiMicroMacroDiscrete::computeCurrIm(int iter)
    {

    cout << " debut compute current image DD" << endl;
    cout << " currentpoints list size is " << currentImagePointsList.pointsList->size() << endl;
    double t1, t2, elapsed_time;

    timeval tim;
    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

    unsigned long long int posX;

    unsigned long long int nbPointsTotal = grid->getNbTotalPoints();

    unsigned long long int **controlCoords = dynsys->getControlIntCoords();

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

    imagePointsList tempPointsList;
    tempPointsList.maxNum = 0;
    tempPointsList.minNum = 0;
    if (this->whichPointListToUse == 1)
	{
	tempPointsList.pointsList = tempPointsList2;
	tempPointsList.pointsList->clear();
	}
    else
	{
	tempPointsList.pointsList = tempPointsList1;
	tempPointsList.pointsList->clear();
	}

    double tempL, tempMu;

    list<imagePoint>::iterator itPoint = currentImagePointsList.pointsList->begin(), itTemp;
    imagePoint tempImagePoint;

    list<imagePoint>::iterator itStart, itNew;
    itStart = tempPointsList.pointsList->begin();
    double imageCoords[dim];

    while (!currentImagePointsList.pointsList->empty())    //(itPoint!=itLastPoint)
	{

	posX = (*itPoint).PointNum;
	//cout<< " posX="<<posX<<endl;
	grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
	//  cout<< " rho= "<<rho;

	/*!
	 * On calcule l'image discrète  du point
	 */

	this->computeDiscreteImageOfPoint(posX);
	/*!
	 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
	 *  de cellules
	 */
	unsigned long long int numPoint, numControl;

	//cout<< "  nb cells dans l'image "<<pointDI_DD.nbImagePoints<<endl;

	// cout<<  " analyse d'une image de point\n";
	for (int i = 0; i < pointDI.nbImagePoints; i++)
	    {
	    numPoint = pointDI.tabImagePoints[i];
	    tempImagePoint.PointNum = numPoint;
	    // cout<< " i= "<<i<<  "numPoint= "<<numPoint<<endl;
	    //cout<< " intervalle "<<pointDI_DD.tabPointsEntrees[i]<< " "<< pointDI_DD.tabPointsEntrees[i+1]<<endl;
	    if (numPoint < nbPointsTotal)
		{
		tempImagePoint.minVal = PLUS_INF;
		for (int j = pointDI.tabPointsEntrees[i]; j < pointDI.tabPointsEntrees[i + 1]; j++)
		    {
		    numControl = pointDI.tabImageControls[j];
		    //cout<< " numControl = "<<numControl<<endl;
		    tempL = dynsys->lFunc_fd(intPointCoords, controlCoords[numControl]);
		    tempMu = dynsys->muFunc_fd(intPointCoords, controlCoords[numControl]);
		    // cout<< " tempL = "<<tempL<< " tempmu = "<<tempMu<<endl;
		    tempImagePoint.minVal = min(tempImagePoint.minVal, max(tempMu, (*itPoint).minVal + tempL));
		    }
		// cout<< " ajout  de cellule avec valeur "<< tempImageCell.minVal<<"\n ";
		addDataToGivenPointsList(&tempPointsList, &itStart, tempImagePoint, &itNew);
		itStart = itNew;
		// cout<< " ares ajout  de l'image d'un point la taille de la liste est "<<currentImageList.cellsList.size()<<endl;
		}
	    }

	//        cout<< " \n";

	itPoint++;
	currentImagePointsList.pointsList->pop_front();
	}

    cout << " parcours de base terminé\n";
    currentImagePointsList.pointsList = tempPointsList.pointsList;

    if (this->whichPointListToUse == 1)
	{
	whichPointListToUse = 2;
	}
    else
	{
	whichPointListToUse = 1;
	}

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));

    cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;

    }

void ViabiMicroMacroDiscrete::addDataToPointsList(list<imagePoint>::iterator *startIt, imagePoint newPoint, list<imagePoint>::iterator *resIt)
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

void ViabiMicroMacroDiscrete::addDataToGivenPointsList(imagePointsList *tempImagePointsList, list<imagePoint>::iterator *startIt, imagePoint newPoint,
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

void ViabiMicroMacroDiscrete::addDataToPoint(list<imagePoint>::iterator itCell, imagePoint newPoint)
    {
    (*itCell).minVal = min((*itCell).minVal, newPoint.minVal);
    }

void ViabiMicroMacroDiscrete::createPointsList()
    {

    }

void ViabiMicroMacroDiscrete::addDataToCurrentImage(list<imageCell>::iterator *startIt, imageCell newCell, list<imageCell>::iterator *resIt)
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

void ViabiMicroMacroDiscrete::addDataToCell(list<imageCell>::iterator itCell, imageCell newCell)
    {

    (*itCell).minVal = min((*itCell).minVal, newCell.minVal);

    }

/*!
 * Le destructeur nettoie la mémoire réservée pour certaines variables globales é la classe et ferme
 * les bases de données de rétro-actions
 */
ViabiMicroMacroDiscrete::~ViabiMicroMacroDiscrete()
    {

    cout << " desturction de classe viabi hjb \n";
    // TODO Auto-generated destructor stub

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

void ViabiMicroMacroDiscrete::printViabiInfo()
    {
    grid->printGrid();
    }

void ViabiMicroMacroDiscrete::showCurrentImageList()
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
void ViabiMicroMacroDiscrete::showCurrentImagePointsList()
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

int ViabiMicroMacroDiscrete::addNewPoints()
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

bool ViabiMicroMacroDiscrete::testConstraintesForCell(unsigned long long int numCell)
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

void ViabiMicroMacroDiscrete::initialiseConstraints()
    {
    spdlog::info("Initialization on the contrats set K0");

    unsigned long long int iCoords[dim];
    double xCoords[dim];
    unsigned long long int totalPointsX = grid->getNbTotalPoints();
    for (unsigned long long int pos = 0; pos < totalPointsX; pos++)
	{
	grid->numToIntAndDoubleCoords(pos, iCoords, xCoords);
	vTab[pos] = dynsys->constraintsX_fd(iCoords);
	}
    }

void ViabiMicroMacroDiscrete::ViabilityKernel(bool sortieOK, int nbArret)
    {
    dynsys->setDynamicsForward();

    viabKerValFunc(nbArret);
    }

void ViabiMicroMacroDiscrete::GarantedViabilityKernel(bool sortieOK, int nbArret)
    {
    dynsys->setDynamicsForward();

    viabKerGarantiValFunc(nbArret);
    }

