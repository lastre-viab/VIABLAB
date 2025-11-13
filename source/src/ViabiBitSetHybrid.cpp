/*
 * ViabiBitSetHybrid.cpp
 *
 *  Created on: 21 août 2025
 *      Author: adesi
 */

#include "../include/ViabiBitSetHybrid.h"

ViabiBitSetHybrid::ViabiBitSetHybrid()
    {
    // TODO Auto-generated constructor stub

    }

ViabiBitSetHybrid::ViabiBitSetHybrid(ParametersManager *pm)
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
    dim_c = gp->DIM_HC;
    dim_d = gp->DIM_HD;
    spdlog::info("[Viabi] : ViabiBitSet Hybrid grid initialization");
    grid = new GridBitSetHybrid(*gp);
    dynsys = new SysDyn(*sp, dim_c, dim_d, *cp, grid);

    InitViabiBitSetHybrid(*avp);
    }

void ViabiBitSetHybrid::InitViabiBitSetHybrid(const algoViabiParams &avbp)
    {

    spdlog::info("[Viabi] : ViabiBitSet Hybrid startInit");
    dimC_c = dynsys->getDimC();
    dimC_d = dynsys->getDimHybrid();
    nbPointsCube_c = (int) pow(2.0, dim_c);
    filePrefix = avbp.FILE_PREFIX;
    nbOMPThreads = avbp.NB_OMP_THREADS;
    indicesDecalCell_c = grid->GetContinuousStateShifts();

    /*!
     * On initialise la structure  qui servira au stockage temporaire de l'image discrete d'un point de grille
     * Les tableaux sont initialis�s  avec leur taille maximale possible :
     * �gale au nb  de controles.
     * C'est la valeur effective  de nbImageCells , calcul� � chaque evaluation de l'image discrete
     * d'un point qui  servira  � lire et remplir  correctement
     * la bonne partie de ces tableaux
     */

    pointDI.tabImageCells = new unsigned long long int[dynsys->getTotalNbPointsC()];

    imageCells = new unsigned long long int[dynsys->getTotalNbPointsC()];
    }

void ViabiBitSetHybrid::printViabiInfo()
    {
    grid->printGrid();
    }

void ViabiBitSetHybrid::ViabilityKernel(int nbArret)
    {

    dynsys->setDynamicsForward();
    const algoViabiParams *avp = modelParams->getAlgoParameters();

    ostringstream os;
    string fileName;
    double t1, t2, elapsed_time, t1_glob, t2_glob, elapsed_time_glob;
    timeval tim, tim_glob;

    gettimeofday(&tim, NULL);              //mesure le temps d'execution
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    //  Calcul du noyau de viabilité
    ViabilityKernelSimple(nbArret);
    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));
    cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;

    this->saveViableSets();
    }

void ViabiBitSetHybrid::ViabilityKernelSimple(int nbArret)
    {

    noyauViabi(nbArret);

    }

void ViabiBitSetHybrid::noyauViabi(int nbArret)
    {

    int dirTramage = grid->getDirTram();

    const unsigned long long int *nbPointsSub = grid->getNbPointsSubGrid();

    bool testF;

    unsigned long long int longTrame = grid->getLongTrame();
    boost::dynamic_bitset<> masque;
    boost::dynamic_bitset<> *masquePointsEnleves = new boost::dynamic_bitset<>(longTrame, 0);
    boost::dynamic_bitset<> **gridTab = grid->getGridTab();

    double *limInf = grid->limInf;
    double *gridStep = grid->step;

    unsigned long long int indice[dim - 1];

    double xCoordsDouble[dim];
    unsigned long long int xCoordsInt[dim];
    bool testNonVide;
    int nbIter = 0;
    unsigned long long int comptEtats = 0, comptEnleves = nbArret + 1;
    unsigned long long int subGridSize = grid->getNbPointsTotalSubGrid();

    unsigned long long int * currentDiscreteState = new unsigned long long int[dim_d];
    double * currentContinuousState = new double[dim_c];

    testK0();

    while (comptEnleves > (unsigned long long int) nbArret)
	{
	cout << "nouvelle boucle while\n";

	comptEnleves = 0;

	comptEtats = 0;
	for (unsigned long long int posX = 0; posX < subGridSize; posX++)
	    {

	    if (!gridTab[posX]->none())
		{
		comptEtats++;
		testF = false;
		numToIntCoords_gen(posX, dim - 1, nbPointsSub, indice);

		masque = grid->analyseTrameMasque(posX);

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
			xCoordsInt[j] = indice[j];
			}

		    for (int j = dirTramage + 1; j < dim; j++)
			{
			xCoordsDouble[j] = limInf[j] + indice[j - 1] * gridStep[j];
			xCoordsInt[j] = indice[j - 1];
			}

		    for (unsigned long long int k = 0; k < longTrame; k++)
			{
			if (masque[k])
			    {
			    xCoordsDouble[dirTramage] = limInf[dirTramage] + k * gridStep[dirTramage];
			    xCoordsInt[dirTramage] = k;
//			    logVector("[ViabiHybrid] : coordsDouble: ", xCoordsDouble, dim);
//			    logVector("[ViabiHybrid] : coordsint: ", xCoordsInt, dim);

			    for(int j = 0; j< dim_c; j++)
				{
				currentContinuousState[j] = xCoordsDouble[j];
				}
			    for(int i = 0; i < dim_d; i++)
				{
				currentDiscreteState[i] = xCoordsInt[dim_c + i];
				}

//			    logVector("[ViabiHybrid] : ===========cont state: ", currentContinuousState, 2);
//			   logVector("[ViabiHybrid] : ============discrete state: ",currentDiscreteState, 2);
			    testNonVide = this->findViabImagePoint(currentContinuousState, currentDiscreteState, false);

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

bool ViabiBitSetHybrid::findViabImagePoint(double *xCoordsDouble, unsigned long long int *xCoordsInt, bool print)
    {
    print = false;

    /*
     * Contrinuous controls grid
     */
    double **controlCoords = dynsys->getControlCoords();
    //Discrete controls grid
    unsigned long long int **hybridControlCoords = dynsys->getHybridIntCoords();

    /*
     * Numbers of control points : continuous and discrete
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int nbCHybridTotal = dynsys->getTotalNbPointsHybrid();

    bool testNonVide = false;

    double *doubleImage = new double[dim_c];
    double *tempReset = new double[dim_c];
    unsigned long long int *tempResetd = new unsigned long long int[dim_d];
    //TODO : calcul de rho hybrid
    double rho = dynsys->calculRho_local_hybrid(xCoordsDouble, xCoordsInt);
    unsigned long long int *intImage = new unsigned long long int[dim_d];
    unsigned long long int *gridCoords = new unsigned long long int[dim];

    unsigned long long int cu_c = 0;
    unsigned long long int cu_d = 0;

    /*
     * On recherche la plus grande puissance de 2 plus grande  que le nombre totale de
     * contrôles à parcourir
     */

        //logVector("[ViabiHybrid] : xDouble: ", xCoordsDouble, 2);
        //    logVector("[ViabiHybrid] : xint: ",xCoordsInt, 2);
    while ((cu_c < nbCTotal) && !testNonVide)
	{
	cu_d = 0;
	while ((cu_d < nbCHybridTotal) && !testNonVide)
	    {

	    testNonVide = this->CheckViability(xCoordsDouble, xCoordsInt, controlCoords[cu_c], hybridControlCoords[cu_d], doubleImage, intImage,
		    gridCoords, tempReset, tempResetd, rho);
	    cu_d++;
	    }
	cu_c ++;

	}			//fin de parcours de tous les contrôles
    delete[] doubleImage;
    delete[] intImage;
    delete[] gridCoords;
    delete[] tempReset;

    return testNonVide;

    }

bool ViabiBitSetHybrid::CheckViability(double *xCoordsDouble, unsigned long long int *xCoordsInt, double *uc, unsigned long long int *ud,
	double *doubleImage, unsigned long long int *intImage, unsigned long long int *gridCoords, double * tempResetc, unsigned long long int * tempResetd, double rho)
    {
    bool testNonVide = false;

    unsigned long long int cellNum, posTemp;
    if (dynsys->constraintsXU_hybrid(xCoordsDouble, xCoordsInt, uc, ud) < PLUS_INF)
	{

	(dynsys->*(dynsys->discretDynamics_hybrid))(xCoordsDouble, xCoordsInt, uc, ud, doubleImage, intImage, tempResetc, tempResetd, rho);
//	logVector("[ViabiHybrid] :Image Double: ", doubleImage, 2);
//	logVector("[ViabiHybrid] : Image int: ",intImage, 2);
	if (grid->isHybridPointInGrid(doubleImage, intImage))
	    {
	    /*
	     *  The sucessor is in the grid
	     */
	    //Only continuous state is checked for state constraints,
	    //as the discrete state is supposed to be in a finite set
	    if (dynsys->constraintsX(doubleImage) < PLUS_INF)
		{

		cellNum = grid->LocalizeHybridPoint(doubleImage, intImage);
		//cout<< " num cellule "<<cellNum<<endl;

		/*
		 * On parcours les sommets de la maille
		 * autour du sucesseur pour voir s'il y a des
		 * points viables
		 */
		int ii = 0;
		while (ii < nbPointsCube_c && !testNonVide)
		    {
		    posTemp = cellNum + indicesDecalCell_c[ii];
		    grid->numToIntCoords(posTemp, gridCoords);
		   // logVector("[ViabiHybrid] :coords du voisin: ", gridCoords, dim);
		    testNonVide = grid->isInSet(gridCoords);
		    ii++;
		    }
		}
	    }
	else
	    {
	    testNonVide = grid->unboundedDomain && grid->isPointInGridWithConstr(doubleImage) && (dynsys->constraintsX(doubleImage) < PLUS_INF);
	    }

	}
    return testNonVide;
    }


void ViabiBitSetHybrid::initialiseConstraints()
    {
    setK0();
    }


void ViabiBitSetHybrid::setK0()
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


void ViabiBitSetHybrid::testK0()
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

void ViabiBitSetHybrid::CaptureBasin()
    {

    }

void ViabiBitSetHybrid::GuarantedViabilityKernelSimple(int nbArret)
    {
    noyauViabiGuaranti(nbArret);
    }



void ViabiBitSetHybrid::noyauViabiGuaranti(int nbArret)
    {


    }


string ViabiBitSetHybrid::getSetName(SetType type)
    {
    string setName = toString(type);
    transform(setName.begin(), setName.end(), setName.begin(), [](char c) {
	return tolower(c);
    });
    // On veut que VIABG renvoie "viabG" (et non pas "viabg")
    if (type == VIABG) setName[4] = 'G';
    return setName;
    }

void ViabiBitSetHybrid::loadViableSets()
    {

    const algoViabiParams *avp = modelParams->getAlgoParameters();

    ostringstream os;
    string fileName;
    string setName = getSetName(avp->SET_TYPE);

    // on charge dans la mémoire l'ensemble calculé et enregistré
    // correspondant au dernier raffinement
    os << "../OUTPUT/" << filePrefix << "-" << setName << ".dat";
    fileName = os.str();
    os.str("");
    grid->loadSetHybrid(fileName);

    }

void ViabiBitSetHybrid::saveViableSets(const string &baseFilenameSuffix) {
    ostringstream os;
    string fileName;
    const algoViabiParams *avp = modelParams->getAlgoParameters();
    string setName = getSetName(avp->SET_TYPE);

    os << "../OUTPUT/" << filePrefix << "-" << setName << baseFilenameSuffix << ".dat";
    fileName = os.str();
    os.str("");
    if (avp->SAVE_VIAB_LIGHT)
	{
	grid->saveValOnGridLightHybrid(fileName);
	}
    else
	{
	grid->saveValOnGridHybrid(fileName);
	}
    if (avp->SAVE_SLICE)
	{
	os << "../OUTPUT/" << filePrefix << "-" << setName << baseFilenameSuffix << "-Slice.dat";
	fileName = os.str();
	os.str("");
	grid->saveHybridCoupe(fileName);
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


void ViabiBitSetHybrid::computeTrajectories()
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

void ViabiBitSetHybrid::computeViableTrajectories()
    {
    //    ostringstream os;
    //    string fileName;
    //
    //    const int nbTrajs = modelParams->getNbTrajectories();
    //
    //    for (int tr = 0; tr < nbTrajs; tr++)
    //    {
    //        TrajectoryParametersManager tpm{modelParams, tr};
    //        const trajectoryParams *tp = tpm.getTrajectoryParameters();
    //        TypeTraj typeTraj = tp->TRAJECTORY_TYPE;
    //        double T = tp->maxTime;
    //        ViabiBitSetTrajectoryHelper trajHelper(grid, dynsys, &tpm);
    //
    //        if (isViableDefault(typeTraj))
    //        {
    //            os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1 << ".dat";
    //            fileName = os.str();
    //            os.str("");
    //
    //            trajHelper.computeViableTrajectory(tp->INIT_POINT, T, fileName);
    //        }
    //        else if (typeTraj == VL)
    //        {
    //            os << "../OUTPUT/" << filePrefix << "-traj-H-" << tr + 1 << ".dat";
    //            fileName = os.str();
    //            os.str("");
    //
    //            trajHelper.computeViableTrajectoryHeavy(tp->INIT_POINT, tp->INIT_CONTROL, T, fileName);
    //        }
    //        else if (typeTraj == CAUTIOUS || typeTraj == WEIGHTED_CONTROLS_CAUTIOUS) {
    //            os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1 << ".dat";
    //            fileName = os.str();
    //            os.str("");
    //
    //            trajHelper.computeCautiousTrajectory(tp->INIT_POINT, T, fileName);
    //        }
    //        else if (typeTraj == CAUTIOUS_HEAVY || typeTraj == WEIGHTED_CONTROLS_CAUTIOUS_HEAVY) {
    //            os << "../OUTPUT/" << filePrefix << "-traj-H-" << tr + 1 << ".dat";
    //            fileName = os.str();
    //            os.str("");
    //
    //            trajHelper.computeCautiousTrajectoryHeavy(tp->INIT_POINT, tp->INIT_CONTROL, T, fileName);
    //        }
    //        else if (typeTraj == STRATEGY_LIST) {
    //            os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1 << ".dat";
    //            fileName = os.str();
    //            os.str("");
    //
    //            trajHelper.computeStrategyTrajectory(tp->INIT_POINT, T, fileName);
    //        }
    //        else {
    //            spdlog::warn("Unimplemeted trajectory type. No trajectory calculation.");
    //        }
    //	    }
    }

void ViabiBitSetHybrid::computeTychasticViableTrajectories() {
    //    ostringstream os;
    //    string fileName;
    //
    //    const int nbTrajs = modelParams->getNbTrajectories();
    //
    //    for (int tr = 0; tr < nbTrajs; tr++)
    //    {
    //        TrajectoryParametersManager tpm{modelParams, tr};
    //        const trajectoryParams *tp = tpm.getTrajectoryParameters();
    //        TypeTraj typeTraj = tp->TRAJECTORY_TYPE;
    //        double T = tp->maxTime;
    //        ViabiBitSetTrajectoryHelper trajHelper(grid, dynsys, &tpm);
    //        TychePicker tychePicker(&tpm, dynsys);
    //
    //        if (typeTraj == STRATEGY_LIST) {
    //            os << "../OUTPUT/" << filePrefix << "-traj-" << tr + 1 << ".dat";
    //            fileName = os.str();
    //            os.str("");
    //
    //            trajHelper.computeTychasticStrategyTrajectory(tp->INIT_POINT, T, fileName, tychePicker);
    //        }
    //        else {
    //            spdlog::error("Unimplemented trajectory type. No trajectory calculation");
    //        }
    //    }
}

void ViabiBitSetHybrid::initialiseTarget()
    {

    /*!
     *  cette fonction initialise la base de donn�es pour la dynamique. Elle
     *  ajoute dnas la base les premiers points pour lesquels la fonction valeur
     *  est r�elle. Au d�but de l'algorithme de bassin de capture
     *    seuls les points de la  cible  ont une fonction valeur r�elle
     *
     */
    // initialiseTargetPointList();
    }

void ViabiBitSetHybrid::GarantedViabilityKernel(int nbArret)
    {
    dynsys->setDynamicsForward();
    const algoViabiParams *avp = modelParams->getAlgoParameters();

    ostringstream os;
    string fileName;
    double t1, t2, elapsed_time, t1_glob, t2_glob, elapsed_time_glob;
    timeval tim, tim_glob;

    gettimeofday(&tim, NULL);              //mesure le temps d'execution
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    //  Calcul du noyau de viabilité
    cout << "=============================Debut viab kernel ==============================" << endl;
    GuarantedViabilityKernelSimple(nbArret);
    cout << "==============================================================================" << endl;
    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    elapsed_time = (double) ((t2 - t1));
    cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;


    this->saveViableSets();
    }

SysDyn* ViabiBitSetHybrid::GetSysDynForViabProblem()
    {
    return this->dynsys;
    }

void ViabiBitSetHybrid::saveViableSets()
    {
    saveViableSets("");
    }
ViabiBitSetHybrid::~ViabiBitSetHybrid()
    {
    // TODO Auto-generated destructor stub
    }

