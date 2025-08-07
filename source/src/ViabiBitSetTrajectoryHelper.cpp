#include <numeric> // pour std::iota

#include "../include/ViabiBitSetTrajectoryHelper.h"

#include "../include/TrajectoryStorage.h"
#include "../include/TychasticTrajectoryStorage.h"
#include "../include/TrajectoryPointsStorage.h"
#include "../include/ControlPickerBitSet.h"
#include "../include/TychasticControlPickerBitSet.h"

ViabiBitSetTrajectoryHelper::ViabiBitSetTrajectoryHelper(Grid_BitSet *gr, SysDyn *ds, TrajectoryParametersManager *tpm) :
    grid(gr),
    tpm(tpm),
    dynsys(ds) {
    const trajectoryParams *tp = tpm->getTrajectoryParameters();
    TypeTraj typeTraj = tp->TRAJECTORY_TYPE;

    dim = grid->dim;
    dimC = dynsys->getDimC();

    unsigned long long int nbTotalC = dynsys->getTotalNbPointsC();    
    preferedControlIndexes = new int[nbTotalC];
    std::iota(preferedControlIndexes, preferedControlIndexes+nbTotalC, 0);
    controlWeight = tp->CONTROL_WEIGHT;
    sortIndexes = tp->SORT_INDEXES;
    trajIndex = tpm->getTrajectoryIndex();
}

ViabiBitSetTrajectoryHelper::~ViabiBitSetTrajectoryHelper() {
    delete [] preferedControlIndexes;
}

int *ViabiBitSetTrajectoryHelper::getPreferedControlIndexes() {
    return preferedControlIndexes;
}

int *ViabiBitSetTrajectoryHelper::sortPreferedControlIndexes(const double *x, double t, int strategyIndex) {
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    double **controlCoords = dynsys->getControlCoords();
    sortIndexes(preferedControlIndexes, x, controlCoords, nbCTotal, t, controlWeight, trajIndex, strategyIndex);
    return preferedControlIndexes;
}

int ViabiBitSetTrajectoryHelper::computeViableTrajectoryHeavy(double *initPosition, double *initControl, double finalTime, string fileName)
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
    const long long int *indicesDecalCell = grid->getIndicesDecalCell();
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
                cout << " La position initiale sélectionnée n'appartiant pas au noyau de viabilité\n";
                succes = 0;
            }
            else
            {
                TrajectoryStorage traj = TrajectoryStorage::createNonFlagSavingStorage(initPosition, finalTime, dynsys);
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
                int currentCu = getClosestControlTo(currentControl);

                while (time < finalTime && nbIter < NB_MAX_TRAJ_ITER)
                {
                    //  cout<< " temps= "<<newTrajPoint[dim]<<endl;

                    rho = dynsys->calculRho_local(xCoordsDouble);

                    rho = min(rho, finalTime - time);
                    // cout<< " rho= "<<rho<<endl;
                    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, currentControl, imageVect, rho);
                    testNonVide = false;
                    time += rho;
                    if (grid->isPointInGrid(imageVect))
                    {
                        /*
                         *  le sucesseur est dans la grille de calcul
                         */
                        if (dynsys->constraintsX(imageVect) < PLUS_INF)
                        {
                            cellNum = grid->localizePoint(imageVect);
                            // cout<< " num cellule "<<cellNum<<endl;

                            cptOK = 0;
                            for (int ii = 0; ii < nbPointsCube; ii++)
                            {
                                posTemp = cellNum + indicesDecalCell[ii];
                                grid->numToIntAndDoubleCoords(posTemp, testI, testV);
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
                            for (int i = 0; i < dim; i++)
                            {
                                xCoordsDouble[i] = imageVect[i];
                            }
                            traj.addPoint(currentCu, imageVect, time);
                        }
                        else
                        {
                            cout << " ======================= Recalage pendant la phase controle constant  =======================\n";
                            grid->findNearestViabPointInCell(xCoordsDouble, imageVect, xCoordsDouble, dynsys->dynConstraintsForTraj);
                            traj.addPoint(currentCu, imageVect, time);

                        }
                    }
                    else
                    {
                        cout << " ======================= CHANGEMENT DE CONTROLE =======================\n";

                        currentCu = this->findViabControl(xCoordsDouble, rho, 1, 1.0, imageVect, maxnbViabPoints, testNonVide);

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
                                    currentControl[dc] = controlCoords[currentCu][dc];
                                }
                                for (int i = 0; i < dim; i++)
                                {
                                    xCoordsDouble[i] = imageVect[i];
                                }
                                traj.addPoint(currentCu, imageVect, time);
                            }
                            else
                            {
                                cout << " ======================= Recalage =======================\n";
                                grid->findNearestViabPointInCell(xCoordsDouble, imageVect, xCoordsDouble, dynsys->dynConstraintsForTraj);

                                cout << " controle ";
                                for (int dc = 0; dc < dimC; dc++)
                                {
                                    currentControl[dc] = controlCoords[currentCu][dc];
                                    cout << " " << currentControl[dc];
                                }
                                traj.addPoint(currentCu, imageVect, time);

                            }
                        }
                        else
                        {
                            cout << "   Echec! Sortie de l'ensemble viable \n";
                            break;

                        }
                    }                    
                    nbIter++;
                }
                if (time >= finalTime)
                {
                    succes = 1.0;
                }
                traj.writeToFile(fileName);
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

    return succes;
}

int ViabiBitSetTrajectoryHelper::findViabControl(double *currentPos, double &dt, int nbStepIter, double stepCoeff, double *resPos, int &nbViabVoisins, bool &succes)
    {

    int nbPointsCube = (int) pow(2.0, dim);			//pow(2.0, dim);
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */
    unsigned long long int testI[dim];
    double testV[dim];
    double imageVect[dim];
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
            if (dynsys->isViableControl(xCoordsDouble, controlCoords[cu], imageVect, rho))
            {
                /*
                 * le successeur vérifie les contraintes
                 * On identifie la maille où il se trouve
                 */
                cellNum = grid->localizePoint(imageVect);

                /*
                 * On parcours les sommets de la maille
                 * autour du sucesseur pour voir s'il y a des
                 * points viables
                 */
                cptOK = 0;
                for (int ii = 0; ii < nbPointsCube; ii++)
                {
                    posTemp = cellNum + indicesDecalCell[ii];
                    grid->numToIntAndDoubleCoords(posTemp, testI, testV);
                    if (dynsys->dynConstraintsForTraj(xCoordsDouble, testV) < PLUS_INF)	// on v�rifie si la projection corrspond aux contraintes de dynamique
                    {
                        double dist = 0.0;
                        for (int k = 0; k < dim; k++)
                        {
                            dist = max(dist, abs(testV[k] - imageVect[k]));
                        }
                        testviabInt = (grid->isInSet(testI) && (dist <= hMax / 2.0));
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
            cu++;
	    }				//fin de parcours de tous les contrôles
        // la boucle s'arête ici u premier contrôle
        // qui donne un successeur viable
        iter++;
        rho = rho * stepCoeff;
        //   cout<< " iteration  "<<iter<<" rho= "<<rho<<endl;
	}

    succes = testNonVide;
    nbViabVoisins = maxnbViabPoints;
    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[bestCu], imageVect, rho);
    for (int i = 0; i < dim; i++)
	{
        resPos[i] = imageVect[i];
	}
    return bestCu;
    }

unsigned long long int ViabiBitSetTrajectoryHelper::findViableDiscreteSuccessor(unsigned long long int pos, double time, double dt, int &usedCu)
    {
    unsigned long long int viabSuccessor = grid->nbTotalPoints + 1;
    unsigned long long int intCoordsOnDiscreteTraj[dim];
    /*
     * coordonnées du point courant de la trajectoire
     */
    double xCoordsDouble[dim];
    grid->numToIntAndDoubleCoords(pos, intCoordsOnDiscreteTraj, xCoordsDouble);

    double imageVect[dim];

    /*
     * tableau de coordonnées  définissant
     * la grille de contrôles
     */
    
    int *preferedControlIndexes = sortPreferedControlIndexes(xCoordsDouble, time);
    double **controlCoords = dynsys->getControlCoords();

    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();

    /*
     * le pas de temps déterminé localement pour chaque point de la trajectoire
     */
    double rho;

    rho = dt;

    unsigned long long int cPrefu = 0;

    bool testNonVide = false;
    while (cPrefu < nbCTotal && !testNonVide) {
        /*
         * on ne choisit que ceux qui vérifient les éventuelles contraintes mixtes
         * au point en cours
         */                
        if (dynsys->isViableControl(xCoordsDouble, controlCoords[preferedControlIndexes[cPrefu]], imageVect, rho)) {
            viabSuccessor = grid->getNearestPointInSet(imageVect);
            testNonVide = (viabSuccessor < grid->nbTotalPoints + 1);
            usedCu = preferedControlIndexes[cPrefu];
        }
        cPrefu++;
    }
    return viabSuccessor;
    }

unsigned long long int ViabiBitSetTrajectoryHelper::findViableDiscreteSuccessor_tych(unsigned long long int pos, double time, double dt, int tychIndex, int &usedCu) {
        unsigned long long int viabSuccessor = grid->nbTotalPoints + 1;
    unsigned long long int intCoordsOnDiscreteTraj[dim];
    /*
     * coordonnées du point courant de la trajectoire
     */
    double xCoordsDouble[dim];
    grid->numToIntAndDoubleCoords(pos, intCoordsOnDiscreteTraj, xCoordsDouble);

    double imageVect[dim];

    /*
     * tableau de coordonnées  définissant
     * la grille de contrôles
     */
    
    int *preferedControlIndexes = sortPreferedControlIndexes(xCoordsDouble, time);
    double **controlCoords = dynsys->getControlCoords();
    double **tycheCoords = dynsys->getTychCoords();

    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();

    /*
     * le pas de temps déterminé localement pour chaque point de la trajectoire
     */
    double rho;

    rho = dt;

    unsigned long long int cPrefu = 0;

    bool testNonVide = false;
    while (cPrefu < nbCTotal && !testNonVide) {
        /*
         * on ne choisit que ceux qui vérifient les éventuelles contraintes mixtes
         * au point en cours
         */                
        if (dynsys->isViableControl_tych(xCoordsDouble, controlCoords[preferedControlIndexes[cPrefu]], tycheCoords[tychIndex], imageVect, rho)) {
            viabSuccessor = grid->getNearestPointInSet(imageVect);
            testNonVide = (viabSuccessor < grid->nbTotalPoints + 1);
            usedCu = preferedControlIndexes[cPrefu];
        }
        cPrefu++;
    }
    return viabSuccessor;
}


int ViabiBitSetTrajectoryHelper::findViabControl_bis(const double *currentPos, unsigned long long int optimDiscreteSuccessor, double &dt, int nbStepIter, double stepCoeff,
	double *resPos, bool &succes)
    {
    double realCoordsOnDiscreteTraj[dim];
    unsigned long long int intCoordsOnDiscreteTraj[dim];
    double rho0 = dt;
    int bestCu;
    const long long int *indicesDecalCell = grid->getIndicesDecalCell();
    grid->numToIntAndDoubleCoords(optimDiscreteSuccessor, intCoordsOnDiscreteTraj, realCoordsOnDiscreteTraj);

    double imageVect[dim], testV[dim];
    unsigned long long int testI[dim];
    int nbPointsCube = (int) pow(2.0, dim);
    double **controlCoords = dynsys->getControlCoords();
    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int posTemp;
    /*
     * coordonnées du points corent de la trajectoire
     */
    double xCoordsDouble[dim];
    ostringstream os;
    string msg;

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;

    spdlog::debug("Start of computing of the real viable control");

    for (int i = 0; i < (int) dim; i++)
	{
	xCoordsDouble[i] = currentPos[i];
	}
    os << "Time step  : " << rho0 << ". Current point = ";
    msg = os.str();
    os.str("");
    logVector(msg, xCoordsDouble, dim);

    bestCu = 0;
    int iter = 0;
    double rhoOptim = rho0;
    double rhoMin = rho0 * (1 - stepCoeff);			// we will test steps from rhoMin to  rhoMax = rho0 * (1 + stepCoeff);
    double rho = rhoMin;
    double step = rho0 * stepCoeff / nbStepIter;
    testNonVide = false;
    double minDistSquared = PLUS_INF;
    while (iter < 2 * nbStepIter + 1)
	{
        unsigned long long int cu = 0;

        while (cu < nbCTotal)
        {
            /*
             * on ne choisit que ceux qui vérifient les éventuelles contraintes mixtes
             * au point en cours
             */
            if (dynsys->isViableControl(xCoordsDouble, controlCoords[cu], imageVect, rho)) {
                unsigned long long int cellNum = grid->localizePoint(imageVect);

                /*
                 * On parcours les sommets de la maille
                 * autour du sucesseur pour voir s'il y a des
                 * points viables
                 */
                bool allInSet = false;
                for (int ii = 0; ii < nbPointsCube; ii++)
                {
                    posTemp = cellNum + indicesDecalCell[ii];
                    grid->numToIntAndDoubleCoords(posTemp, testI, testV);

                    allInSet |= grid->isInSet(testI);

                }
                if (allInSet)
                {
                    testNonVide = true;
                    double distSquared = 0.0;
                    for (int k = 0; k < dim; k++)
                    {
                        distSquared += (imageVect[k] - realCoordsOnDiscreteTraj[k]) * (imageVect[k] - realCoordsOnDiscreteTraj[k]);
                    }

                    if (distSquared < minDistSquared)
                    {

                        minDistSquared = distSquared;
                        bestCu = cu;
                        rhoOptim = rho;
                    }
                }
            }
            cu++;
        }		//fin de parcours de tous les contrôles

        iter++;
        rho += step;
	}

    if (testNonVide)
	{
        spdlog::info("Best real control found : Distance {}, optimal time step {}", sqrt(minDistSquared), rhoOptim);
	}
    succes = testNonVide;

    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[bestCu], imageVect, rhoOptim);
    for (int i = 0; i < (int) dim; i++)
	{
        resPos[i] = imageVect[i];
	}
    dt = rhoOptim;
    return bestCu;
    }

int ViabiBitSetTrajectoryHelper::findViabControl_bis_tych(const double *currentPos, unsigned long long int optimDiscreteSuccessor, double &dt, int tychIndex, int nbStepIter, double stepCoeff,
	double *resPos, bool &succes)
    {
    double realCoordsOnDiscreteTraj[dim];
    unsigned long long int intCoordsOnDiscreteTraj[dim];
    double rho0 = dt;
    int bestCu;
    const long long int *indicesDecalCell = grid->getIndicesDecalCell();
    grid->numToIntAndDoubleCoords(optimDiscreteSuccessor, intCoordsOnDiscreteTraj, realCoordsOnDiscreteTraj);

    double imageVect[dim], testV[dim];
    unsigned long long int testI[dim];
    int nbPointsCube = (int) pow(2.0, dim);
    double **controlCoords = dynsys->getControlCoords();
    double **tycheCoords = dynsys->getTychCoords();
    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int posTemp;
    /*
     * coordonnées du points corent de la trajectoire
     */
    double xCoordsDouble[dim];
    ostringstream os;
    string msg;

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;

    spdlog::debug("Start of computing of the real viable control");

    for (int i = 0; i < (int) dim; i++)
	{
	xCoordsDouble[i] = currentPos[i];
	}
    os << "Time step  : " << rho0 << ". Current point = ";
    msg = os.str();
    os.str("");
    logVector(msg, xCoordsDouble, dim);

    bestCu = 0;
    int iter = 0;
    double rhoOptim = rho0;
    double rhoMin = rho0 * (1 - stepCoeff);			// we will test steps from rhoMin to  rhoMax = rho0 * (1 + stepCoeff);
    double rho = rhoMin;
    double step = rho0 * stepCoeff / nbStepIter;
    testNonVide = false;
    double minDistSquared = PLUS_INF;
    while (iter < 2 * nbStepIter + 1)
	{
        unsigned long long int cu = 0;

        while (cu < nbCTotal)
        {
            /*
             * on ne choisit que ceux qui vérifient les éventuelles contraintes mixtes
             * au point en cours
             */
            if (dynsys->isViableControl_tych(xCoordsDouble, controlCoords[cu], tycheCoords[tychIndex], imageVect, rho)) {
                unsigned long long int cellNum = grid->localizePoint(imageVect);

                /*
                 * On parcours les sommets de la maille
                 * autour du sucesseur pour voir s'il y a des
                 * points viables
                 */
                bool allInSet = false;
                for (int ii = 0; ii < nbPointsCube; ii++)
                {
                    posTemp = cellNum + indicesDecalCell[ii];
                    grid->numToIntAndDoubleCoords(posTemp, testI, testV);

                    allInSet |= grid->isInSet(testI);

                }
                if (allInSet)
                {
                    testNonVide = true;
                    double distSquared = 0.0;
                    for (int k = 0; k < dim; k++)
                    {
                        distSquared += (imageVect[k] - realCoordsOnDiscreteTraj[k]) * (imageVect[k] - realCoordsOnDiscreteTraj[k]);
                    }

                    if (distSquared < minDistSquared)
                    {

                        minDistSquared = distSquared;
                        bestCu = cu;
                        rhoOptim = rho;
                    }
                }
            }
            cu++;
        }		//fin de parcours de tous les contrôles

        iter++;
        rho += step;
	}

    if (testNonVide)
	{
        spdlog::info("Best real control found : Distance {}, optimal time step {}", sqrt(minDistSquared), rhoOptim);
	}
    succes = testNonVide;

    (dynsys->*(dynsys->discretDynamics_tych))(xCoordsDouble, controlCoords[bestCu], tycheCoords[tychIndex], imageVect, rhoOptim);
    for (int i = 0; i < (int) dim; i++)
	{
        resPos[i] = imageVect[i];
	}
    dt = rhoOptim;
    return bestCu;
    }


int ViabiBitSetTrajectoryHelper::computeViableTrajectory(double *initPosition, double finalTime, string fileName)
    {

    bool succes = false;

    unsigned long long int currentPosIntCoords[dim];
    /*
     * coordonnées du points courant de la trajectoire
     */
    double xCoordsDouble[dim], imageVect[dim], doubleCoordsOnDiscreteTraj[dim];

    SysDyn::PointStatus initStatus = dynsys->checkKernelRelation(initPosition);
    if (initStatus == SysDyn::OUTSIDE_DOMAIN) {
        spdlog::error("Initial point is out of the viability domain. Stop.");
    }
    else if (initStatus == SysDyn::OUTSIDE_CONSTRAINTS) {
        spdlog::error("Initial point is out of the constraints set. Stop.");
    }
    else if (initStatus == SysDyn::OUTSIDE_GRID) {
        spdlog::error("Initial point is out of the grid. Stop.");
    }
    else {
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
        unsigned long long int currentDiscreteTrajPos = grid->getNearestPointInSet(initPosition);
        grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
        TrajectoryStorage traj = TrajectoryStorage::createNonFlagSavingStorage(initPosition, finalTime, dynsys);
        TrajectoryPointsStorage trajDiscrete(doubleCoordsOnDiscreteTraj, finalTime, dynsys);
        int realTimeStepsPerDiscreteStep = tpm->getTrajectoryParameters()->REAL_TIME_STEPS_PER_DISCRETE_STEP;
        
        int nbIter = 0;
        ostringstream os;
        string msg;
        os << "Time : " << time << ". Current point = ";
        msg = os.str();
        os.str("");
        logVector(msg, xCoordsDouble, dim);

        succes = true;
        while (succes && time < finalTime && nbIter < NB_MAX_TRAJ_ITER)
        {
            succes = computeViableTrajectoryIteration(currentDiscreteTrajPos,
                                                      xCoordsDouble,
                                                      currentPosIntCoords,
                                                      doubleCoordsOnDiscreteTraj,
                                                      imageVect,
                                                      finalTime, time,
                                                      traj, trajDiscrete,
                                                      realTimeStepsPerDiscreteStep);
            nbIter++;
        }
        if (time < finalTime)
        {
            succes = false;
        }
        string fileNameD(fileName);
        for (int k = 0; k < 4; k++)
        {
            fileNameD.pop_back();
        }
        fileNameD += "-Discrete.dat";
    
        if (succes)
        {
            spdlog::info("Trajectory reconstruction finished with success");
            spdlog::info("The real trajectory will be saved in the file {}, the discrete trajectory saved in the file ", fileName, fileNameD);
        }
        else
        {
            spdlog::warn("Trajectory reconstruction failed before attempting the target set");
            spdlog::warn("The real partial trajectory will be saved in the file {}, the partial discrete trajectory saved in the file ", fileName,
                         fileNameD);
        }

        traj.writeToFile(fileName);
        trajDiscrete.writeToFile(fileNameD);
    }
    
    return succes;

    }

bool ViabiBitSetTrajectoryHelper::isOnBorder(unsigned long long int *currentPos) {
    // On est sur le bord uniquement si parmi les voisins de la cellules (pas les diagonales)
    // il existe un point qui n'est pas dans le noyau de viabilité
    bool onBorder = false;
    int i = 0;
    while (!onBorder && i < dim) {
        currentPos[i]++;
        onBorder |= !(grid->isPointInGrid_fd(currentPos) && grid->isInSet(currentPos));
        currentPos[i] -= 2;
        onBorder |= !(grid->isPointInGrid_fd(currentPos) && grid->isInSet(currentPos));
        currentPos[i]++;
        ++i;
    }
    
    return onBorder;
}

bool ViabiBitSetTrajectoryHelper::findPositionAtBorder(double *xCoordsDouble,
                                       Bubble &bubble,
                                       double &futureTime, double &rho,
                                       bool &onBorder, unsigned long long int &discreteBorderPos, double *borderCoords, int realTimeStepsPerDiscreteStep) {
    
    using ull = unsigned long long int;
    
    onBorder = false;
    
    EmptyTrajectory traj;
    EmptyTrajectoryPoints trajDiscrete;

    int nbIter = 0;

    ull *currentPosIntCoords = new ull[dim];
    double *doubleCoordsOnDiscreteTraj = new double[dim];
    
    discreteBorderPos = grid->getNearestPointInSet(xCoordsDouble);

    grid->numToIntAndDoubleCoords(discreteBorderPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
    
    double *imageVect = borderCoords;
    double *simulatedFuturePos = new double[dim];
    std::copy(xCoordsDouble, xCoordsDouble+dim, simulatedFuturePos);

    rho = dynsys->calculRho_local(doubleCoordsOnDiscreteTraj);
    
    while (bubble.isInBubble(currentPosIntCoords)
           && !onBorder && nbIter < NB_MAX_TRAJ_SIMULATIONS) {
        
        if (!computeViableTrajectoryIteration(
                discreteBorderPos,
                simulatedFuturePos,
                currentPosIntCoords, doubleCoordsOnDiscreteTraj, imageVect,
                PLUS_INF, rho,
                traj, trajDiscrete,
                realTimeStepsPerDiscreteStep)) {
            return false;
        }

        futureTime += rho;
        onBorder = isOnBorder(currentPosIntCoords);
        rho = dynsys->calculRho_local(doubleCoordsOnDiscreteTraj);
        ++nbIter;
    }

    if (nbIter >= NB_MAX_TRAJ_SIMULATIONS) {
        spdlog::warn("Simulation never left bubble after {} iterations. Using the default control", nbIter);
    }
    
    delete[] doubleCoordsOnDiscreteTraj;
    delete[] simulatedFuturePos;
    delete[] currentPosIntCoords;

    return true;
}


bool ViabiBitSetTrajectoryHelper::computeViableTrajectoryIteration(unsigned long long int& currentDiscreteTrajPos,
                                                   double *xCoordsDouble,
                                                   unsigned long long int *currentPosIntCoords, double *doubleCoordsOnDiscreteTraj,
                                                   double *imageVect,
                                                   double finalTime, double &time,
                                                   Trajectory &traj, TrajectoryPoints &trajDiscrete,
                                                   int realTimeStepsPerDiscreteStep) {
    bool testNonVide;
    /*
     * tableau de coordonnées  définissant
     * la grille de contrôles
     */
    double **controlCoords = dynsys->getControlCoords();

    double rho = dynsys->calculRho_local(doubleCoordsOnDiscreteTraj);

    rho = min(rho, finalTime - time);

    int usedCu;
    unsigned long long int optimSuccessor = this->findViableDiscreteSuccessor(currentDiscreteTrajPos, time, rho, usedCu);
    if (optimSuccessor > grid->nbTotalPoints)
    {
        spdlog::error("Optimal successor not found");
        return false;
    }
    else
    {
        double realTimeStep = rho;
        int bestCu = this->findViabControl_bis(xCoordsDouble, optimSuccessor, realTimeStep, 25, 0.5, imageVect, testNonVide);
        if (testNonVide)
        {
            realTimeStep = min(realTimeStep, finalTime - time);
            time += realTimeStep;
            traj.addPoints(bestCu, realTimeStepsPerDiscreteStep, time, imageVect);
            std::copy(imageVect, imageVect+dim, xCoordsDouble);
            currentDiscreteTrajPos = grid->getNearestPointInSet(xCoordsDouble);
            grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
            
            trajDiscrete.addPoint(doubleCoordsOnDiscreteTraj, time);
        }
        else
        {
            spdlog::warn("Search for real viable control failed. Reset the real trajectory with the discret position.");

            time += rho;            
            currentDiscreteTrajPos = optimSuccessor;
            grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
            std::copy(doubleCoordsOnDiscreteTraj, doubleCoordsOnDiscreteTraj+dim, xCoordsDouble);
            
            traj.addPoint(bestCu, doubleCoordsOnDiscreteTraj, time);
            trajDiscrete.addPoint(doubleCoordsOnDiscreteTraj, time);
        }
    }
    return true;
}

int ViabiBitSetTrajectoryHelper::applyClosestToControl(unsigned long long int currentDiscreteTrajPos,
                                       int cu,
                                       double rho, double time,
                                       double *xCoordsDouble,
                                       double *imageVect) {
    using ull = unsigned long long int;
    // On cherche le contrôle courant
    double **controlCoords = dynsys->getControlCoords();
    int *preferedControlIndexes = getPreferedControlIndexes();
    const unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();

    const int indexInPref = std::find(preferedControlIndexes, preferedControlIndexes+nbCTotal, cu) - preferedControlIndexes;
    bool viableControlFound = false;
    // On sait que indexInPref est invalide (précondition de la méthode)
    // On cherche parmi contrôles par ordre de préférence croissant
    // en partant du contrôle du bord lequel fonctionne
    // Ainsi, on privilégie quand même le choix de l'utilisateur en sélectionnant
    // d'abord le contrôle le plus intéressant pour lui
    int prefCu = indexInPref - 1;
    while (!viableControlFound && prefCu >= 0) {
        const double *prefControl = controlCoords[preferedControlIndexes[prefCu]];
        viableControlFound = dynsys->isViableControl(xCoordsDouble, prefControl, imageVect, rho);
        prefCu--;
    }
    prefCu = ((viableControlFound) ? prefCu : indexInPref);
    // Si l'ordre de préférence croissant n'a pas donné de contrôle viable,
    // il reste en partant du contrôle de bord à vérifier par préférence décroissante
    while (!viableControlFound && (ull) prefCu < nbCTotal) {
        const double *prefControl = controlCoords[preferedControlIndexes[prefCu]];
        viableControlFound = dynsys->isViableControl(xCoordsDouble, prefControl, imageVect, rho);
        prefCu++;
    }

    
    return prefCu;
}

int ViabiBitSetTrajectoryHelper::computeCautiousTrajectory(double *initPosition, double finalTime, string fileName) {

    using ull = unsigned long long int;
    
    double **controlCoords = dynsys->getControlCoords();
    const unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    bool succes = false;

    SysDyn::PointStatus initStatus = dynsys->checkKernelRelation(initPosition);
    if (initStatus == SysDyn::OUTSIDE_DOMAIN) {
        spdlog::error("Initial point is out of the viability domain. Stop.");
    }
    else if (initStatus == SysDyn::OUTSIDE_CONSTRAINTS) {
        spdlog::error("Initial point is out of the constraints set. Stop.");
    }
    else if (initStatus == SysDyn::OUTSIDE_GRID) {
        spdlog::error("Initial point is out of the grid. Stop.");
    }
    else {
        int borderCu = 0;
        
        double time = 0.0;
        int nbIter = 0;

        valarray<double> newTrajPoint(dim + 1);
        valarray<double> trajControlCoords(dimC);
        
        ull currentDiscreteTrajPos = grid->getNearestPointInSet(initPosition);
        ull *currentIntCoords = new ull[dim];
        double *doubleCoordsOnDiscreteTraj = new double[dim];
        grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentIntCoords, doubleCoordsOnDiscreteTraj);

        double *xCoordsDouble = new double[dim];
        std::copy(initPosition, initPosition+dim, xCoordsDouble);
        double *imageVect = new double[dim];
        
        TrajectoryStorage traj = TrajectoryStorage::createNonFlagSavingStorage(xCoordsDouble, finalTime, dynsys);
        TrajectoryPointsStorage trajDiscrete(doubleCoordsOnDiscreteTraj, finalTime, dynsys);
        
        double *borderCoordsDouble = new double[dim];        

        succes = true;

        int realTimeStepsPerDiscreteStep = tpm->getTrajectoryParameters()->REAL_TIME_STEPS_PER_DISCRETE_STEP;
        Bubble bubble(tpm, grid);
        
        while (succes && time < finalTime && nbIter < NB_MAX_TRAJ_ITER) {
            // La bulle a touché un bord
            bubble.setCenter(currentIntCoords);
            if (bubble.isTouchingBoundaryInwards()) {

                spdlog::info("Bubble touching");
                bool onBorder = false;
                double rho = dynsys->calculRho_local(doubleCoordsOnDiscreteTraj);
                rho = min(rho, finalTime - time);

                double futureTime = time;
                double futureRho;
                
                ull discreteBorderPos;
                if (!findPositionAtBorder(xCoordsDouble, bubble,
                                          futureTime, futureRho,
                                          onBorder, discreteBorderPos, borderCoordsDouble, realTimeStepsPerDiscreteStep)) {
                    succes = false;
                    break;
                }
                if (!onBorder) {
                    succes = computeViableTrajectoryIteration(currentDiscreteTrajPos,
                                                              xCoordsDouble,
                                                              currentIntCoords, doubleCoordsOnDiscreteTraj, imageVect,
                                                              finalTime, time, traj, trajDiscrete, realTimeStepsPerDiscreteStep);
                }
                else {
                    
                    int usedCu;
                    ull optimBorderSuccessor = findViableDiscreteSuccessor(discreteBorderPos, time, futureRho, usedCu);
                    
                    // Si on a pas trouvé de contrôle viable au bord, on reprend le précédent
                    if (optimBorderSuccessor > grid->nbTotalPoints) {
                        spdlog::warn("Optimal discrete border successor not found. Using previously viable border control");
                    }
                    else {
                        double realTimeStep = futureRho;
                        bool foundBorderCu;
                        borderCu = findViabControl_bis(borderCoordsDouble, optimBorderSuccessor, realTimeStep, 25, 0.5, imageVect, foundBorderCu);
                        
                        if (!foundBorderCu) {
                            spdlog::warn("No viable real control found on border with approximated real trajectory position. Using control from discrete grid position instead");
                            borderCu = usedCu;
                        }
                        rho = realTimeStep;
                    }
                    rho = min(rho, finalTime - time);
                    // Si le contrôle au bord n'est pas viable, on prend le premier contrôle viable
                    // en partant de celui du bord vers celui de l'emplacement actuel
                    if (!dynsys->isViableControl(xCoordsDouble, controlCoords[borderCu], imageVect, rho)) {
                        if ((ull) (borderCu = applyClosestToControl(currentDiscreteTrajPos, borderCu,
                                                                     rho, time,
                                                                     xCoordsDouble, imageVect)) >= nbCTotal) {
                            succes = false;
                            break;
                        }
                    }
                    
                    time += rho;
                    traj.addPoints(borderCu, realTimeStepsPerDiscreteStep, time, imageVect);
                    std::copy(imageVect, imageVect+dim, xCoordsDouble);
                    currentDiscreteTrajPos = grid->getNearestPointInSet(xCoordsDouble);
                    grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentIntCoords, doubleCoordsOnDiscreteTraj);
                    
                    trajDiscrete.addPoint(doubleCoordsOnDiscreteTraj, time);
                }
            }
            else {
                succes = computeViableTrajectoryIteration(currentDiscreteTrajPos,
                                                          xCoordsDouble, currentIntCoords, doubleCoordsOnDiscreteTraj, imageVect,
                                                          finalTime, time, traj, trajDiscrete,
                                                          realTimeStepsPerDiscreteStep);
            }
            ++nbIter;
        }
        
        delete[] borderCoordsDouble;
        delete[] doubleCoordsOnDiscreteTraj;
        delete[] xCoordsDouble;
        delete[] currentIntCoords;
        delete[] imageVect;
        
        if (time >= finalTime) {
            succes = true;
        }

        string fileNameD(fileName);
        for (int k = 0; k < 4; k++)
        {
            fileNameD.pop_back();
        }
        fileNameD += "-Discrete.dat";
    
        if (succes)
        {
            spdlog::info("Trajectory reconstruction finished with success");
            spdlog::info("The real trajectory will be saved in the file {}, the discrete trajectory saved in the file {}",
                         fileName, fileNameD);
        }
        else
        {
            spdlog::warn("Trajectory reconstruction failed before attempting the target set");
            spdlog::warn("The real partial trajectory will be saved in the file {}, the partial discrete trajectory saved in the file {}",
                         fileName, fileNameD);
        }

        traj.writeToFile(fileName);
        trajDiscrete.writeToFile(fileNameD);
    }
    
    return succes;
}


int ViabiBitSetTrajectoryHelper::getClosestControlTo(const double *u) {
    double **controlCoords = dynsys->getControlCoords();
    unsigned long long nbCTotal = dynsys->getTotalNbPointsC();

    double minDist = PLUS_INF;
    int argMinDist = -1;
    
    for (unsigned long long i = 0; i < nbCTotal; ++i) {
        double euclideanDistanceSquared = 0;
        double *uCandidate = controlCoords[i];
        for (int d = 0; d < dimC; ++d) {
            euclideanDistanceSquared += (u[d] - uCandidate[d])*(u[d] - uCandidate[d]);
        }
        if (euclideanDistanceSquared < minDist) {
            minDist = euclideanDistanceSquared;
            argMinDist = i;
        }
    }
    return argMinDist;
}

bool ViabiBitSetTrajectoryHelper::findPositionAtBorderHeavy(double *xCoordsDouble, const double *control,
                                            Bubble &bubble,
                                            double &futureTime, double &rho,
                                            bool &onBorder, unsigned long long int &discreteBorderPos, double *borderCoords) {
    using ull = unsigned long long int;
    
    onBorder = false;

    int nbIter = 0;

    ull *currentPosIntCoords = new ull[dim];
    double *doubleCoordsOnDiscreteTraj = new double[dim];

    discreteBorderPos = grid->getNearestPointInSet(xCoordsDouble);
    grid->numToIntAndDoubleCoords(discreteBorderPos,
                                  currentPosIntCoords,
                                  doubleCoordsOnDiscreteTraj);
    double *imageVect = borderCoords;
    double *simulatedFuturePos = new double[dim];
    std::copy(xCoordsDouble, xCoordsDouble+dim, simulatedFuturePos);

    rho = dynsys->calculRho_local(doubleCoordsOnDiscreteTraj);

    while (bubble.isInBubble(currentPosIntCoords)
           && !onBorder && nbIter < NB_MAX_TRAJ_SIMULATIONS) {

        // Il se peut que l'on dépasse au delà du bord parce que le pas de temps est trop grand
        bool viable = dynsys->isViableControl(xCoordsDouble, control, imageVect, rho);
        bool differentGridPos = grid->getNearestPointInSet(imageVect) != discreteBorderPos;
        while (!viable && differentGridPos) {
            rho /= 2;
            viable = dynsys->isViableControl(xCoordsDouble, control, imageVect, rho);        
            differentGridPos = (grid->getNearestPointInSet(imageVect) != discreteBorderPos);
            // S'il n'y a aucune différence en sortie de boucle,
            // ou que les seules différences sont des sorties du noyau de viabilité
            // on peut supposer que l'on est sur un bord
            onBorder |= !differentGridPos;
        }

        std::copy(imageVect, imageVect+dim, simulatedFuturePos);
        discreteBorderPos = grid->getNearestPointInSet(simulatedFuturePos);
        grid->numToIntAndDoubleCoords(discreteBorderPos,
                                      currentPosIntCoords,
                                      doubleCoordsOnDiscreteTraj);
        onBorder |= isOnBorder(currentPosIntCoords);

        futureTime += rho;
        rho = dynsys->calculRho_local(doubleCoordsOnDiscreteTraj);
        ++nbIter;
    }

    if (nbIter >= NB_MAX_TRAJ_SIMULATIONS) {
        spdlog::warn("Simulation never left bubble after {} iterations. Using the default control", nbIter);
    }
    
    delete[] doubleCoordsOnDiscreteTraj;
    delete[] simulatedFuturePos;
    delete[] currentPosIntCoords;
    return true;
}


int ViabiBitSetTrajectoryHelper::computeCautiousTrajectoryHeavy(double *initPosition, double *initControl, double finalTime, string fileName) {

    using ull = unsigned long long int;
    
    double **controlCoords = dynsys->getControlCoords();
    const ull nbCTotal = dynsys->getTotalNbPointsC();
    bool succes = false;

    SysDyn::PointStatus initStatus = dynsys->checkKernelRelation(initPosition);
    if (initStatus == SysDyn::OUTSIDE_DOMAIN) {
        spdlog::error("Initial point is out of the viability domain. Stop.");
    }
    else if (initStatus == SysDyn::OUTSIDE_CONSTRAINTS) {
        spdlog::error("Initial point is out of the constraints set. Stop.");
    }
    else if (initStatus == SysDyn::OUTSIDE_GRID) {
        spdlog::error("Initial point is out of the grid. Stop.");
    }
    else {
        double time = 0.0;
        int nbIter = 0;
       
        double *imageVect = new double[dim];
        double *xCoordsDouble = new double[dim];
        std::copy(initPosition, initPosition+dim, xCoordsDouble);
        
        double *borderCoordsDouble = new double[dim];

        int currentCu = getClosestControlTo(initControl);
        const double *currentControl = controlCoords[currentCu];
        
        ull currentDiscreteTrajPos = grid->getNearestPointInSet(xCoordsDouble);
        ull *currentPosIntCoords = new ull[dim];
        double *doubleCoordsOnDiscreteTraj = new double[dim];
        grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);

        int realTimeStepsPerDiscreteStep = tpm->getTrajectoryParameters()->REAL_TIME_STEPS_PER_DISCRETE_STEP;
        Bubble bubble(tpm, grid);

        TrajectoryStorage traj = TrajectoryStorage::createNonFlagSavingStorage(xCoordsDouble, finalTime, dynsys);
        TrajectoryPointsStorage trajDiscrete(doubleCoordsOnDiscreteTraj, finalTime, dynsys);

        succes = true;
        
        while (succes && time < finalTime && nbIter < NB_MAX_TRAJ_ITER) {
            
            bubble.setCenter(currentPosIntCoords);
            double rho = dynsys->calculRho_local(doubleCoordsOnDiscreteTraj);
            rho = min(rho, finalTime - time);
            
            if (bubble.isTouchingBoundaryInwards()) {
                spdlog::info("Bubble touching");
                double futureTime = time;
                double futureRho;
                
                ull discreteBorderPos;
                bool onBorder = false;
                findPositionAtBorderHeavy(xCoordsDouble, currentControl, bubble,
                                          futureTime, futureRho,
                                          onBorder, discreteBorderPos, borderCoordsDouble);
                futureRho = min(futureRho, finalTime - time);
                if (!onBorder) {
                    // Si on n'est pas sur le bord, on continue à utiliser le contrôle courant s'il est valide
                    if (dynsys->isViableControl(xCoordsDouble, currentControl, imageVect, futureRho)) {
                        spdlog::info("Not reaching border after simulation. Using heavy control");
                        time += futureRho;
                        std::copy(imageVect, imageVect+dim, xCoordsDouble);
                        currentDiscreteTrajPos = grid->getNearestPointInSet(xCoordsDouble);
                        grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);

                        traj.addPoints(currentCu, realTimeStepsPerDiscreteStep, time, imageVect);
                        trajDiscrete.addPoint(doubleCoordsOnDiscreteTraj, time);
                    }
                    else {
                        spdlog::error("End position of simulation not on border but control not viable ? Unexpected programm state reached");
                        succes = false;
                        break;
                    }  
                }
                else {
                    spdlog::info("Control leads to valid border position");
                    int usedCu;                    
                    ull optimBorderSuccessor = findViableDiscreteSuccessor(discreteBorderPos, time, futureRho, usedCu);
                    // Si on a pas trouvé de contrôle viable au bord, on reprend le précédent
                    if (optimBorderSuccessor >= grid->nbTotalPoints) {
                        spdlog::warn("Optimal discrete border successor not found. Continuing with heavy control");
                    }
                    else {
                        double realTimeStep = futureRho;
                        bool foundCu;
                        int previousCu = currentCu;
                        currentCu = findViabControl_bis(borderCoordsDouble, optimBorderSuccessor, realTimeStep, 10, 0.5, imageVect, foundCu);
                        if (!foundCu) {
                            spdlog::warn("No viable real control found on border with approximated real trajectory position. Using control from discrete grid position instead");
                            currentCu = usedCu;
                        }
                        if (currentCu != previousCu) {
                            currentControl = controlCoords[currentCu];
                            spdlog::warn("Changing control, switching from u[{}] to u[{}]", previousCu, currentCu);
                        }
                        rho = realTimeStep;
                    }

                    rho = min(rho, finalTime - time);
                    
                    // Si le contrôle au bord n'est pas viable, on prend le premier contrôle viable
                    // en partant de celui du bord vers celui de l'emplacement actuel
                    if (!dynsys->isViableControl(xCoordsDouble, currentControl, imageVect, rho)) {
                        if ((ull) (currentCu = applyClosestToControl(currentDiscreteTrajPos, currentCu,
                                                                     rho, time,
                                                                     xCoordsDouble, imageVect)) >= nbCTotal) {
                            succes = false;
                            break;
                        }
                        currentControl = controlCoords[currentCu];
                    }
                    
                    time += rho;
                    std::copy(imageVect, imageVect+dim, xCoordsDouble);
                    currentDiscreteTrajPos = grid->getNearestPointInSet(xCoordsDouble);
                    grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
                    
                    traj.addPoints(currentCu, realTimeStepsPerDiscreteStep, time, imageVect);
                    trajDiscrete.addPoint(doubleCoordsOnDiscreteTraj, time);
                }
            }
            else if (dynsys->isViableControl(xCoordsDouble, currentControl, imageVect, rho)) {
                time += rho;
                std::copy(imageVect, imageVect+dim, xCoordsDouble);
                currentDiscreteTrajPos = grid->getNearestPointInSet(xCoordsDouble);
                grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);

                traj.addPoints(currentCu, realTimeStepsPerDiscreteStep, time, imageVect);
                trajDiscrete.addPoint(doubleCoordsOnDiscreteTraj, time);
            }
            else {
                spdlog::error("Bubble not touching but current control not viable ? Try using a larger bubble");
                succes = false;
                break;
            }                            
        }

        delete[] imageVect;
        delete[] borderCoordsDouble;
        delete[] xCoordsDouble;
        delete[] doubleCoordsOnDiscreteTraj;
        delete[] currentPosIntCoords;

        string fileNameD(fileName);
        for (int k = 0; k < 4; k++)
        {
            fileNameD.pop_back();
        }
        fileNameD += "-Discrete.dat";
    
        if (succes)
        {
            spdlog::info("Trajectory reconstruction finished with success");
            spdlog::info("The real trajectory will be saved in the file {}, the discrete trajectory saved in the file {}", fileName, fileNameD);
        }
        else
        {
            spdlog::warn("Trajectory reconstruction failed before attempting the target set");
            spdlog::warn("The real partial trajectory will be saved in the file {}, the partial discrete trajectory saved in the file {}", fileName,
                         fileNameD);
        }

        traj.writeToFile(fileName);
        trajDiscrete.writeToFile(fileNameD);
    }
    
    return succes;
}

int ViabiBitSetTrajectoryHelper::computeStrategyTrajectory(double *initPosition, double finalTime, string fileName) {

    using ull = unsigned long long int;

    bool succes = false;

    SysDyn::PointStatus initStatus = dynsys->checkKernelRelation(initPosition);
    if (initStatus == SysDyn::OUTSIDE_DOMAIN) {
        spdlog::error("Initial point is out of the viability domain. Stop.");
    }
    else if (initStatus == SysDyn::OUTSIDE_CONSTRAINTS) {
        spdlog::error("Initial point is out of the constraints set. Stop.");
    }
    else if (initStatus == SysDyn::OUTSIDE_GRID) {
        spdlog::error("Initial point is out of the grid. Stop.");
    }
    else {
        ControlPickerBitSet picker(this, dynsys, tpm);
        double **controlCoords = dynsys->getControlCoords();
        
        double time = 0.0;
        int nbIter = 0;

        double *xCoordsDouble = new double[dim];
        std::copy(initPosition, initPosition+dim, xCoordsDouble);
        double *imageVect = new double[dim];
        double *doubleCoordsOnDiscreteTraj = new double[dim];        
        ull *currentPosIntCoords = new ull[dim];
                
        unsigned long long int currentDiscreteTrajPos = grid->getNearestPointInSet(xCoordsDouble);
        grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
        TrajectoryStorage traj = (tpm->getTrajectoryParameters()->SAVE_PICKING_STRATEGY)
            ? TrajectoryStorage::createFlagSavingStorage(initPosition, finalTime, dynsys)
            : TrajectoryStorage::createNonFlagSavingStorage(initPosition, finalTime, dynsys);
        
        TrajectoryPointsStorage trajDiscrete(doubleCoordsOnDiscreteTraj, finalTime, dynsys);        

        int realTimeStepsPerDiscreteStep = tpm->getTrajectoryParameters()->REAL_TIME_STEPS_PER_DISCRETE_STEP;
        succes = true;
        
        while (succes && time < finalTime && nbIter < NB_MAX_TRAJ_ITER) {

            ControlPicker::StrategyIndexBitFlag flag = 0;
            double rho = dynsys->calculRho_local(doubleCoordsOnDiscreteTraj);
            rho = min(rho, finalTime - time);
            
            OptionalCu opt = picker.pickControl(traj, trajDiscrete, rho, flag)
                // On vérifie que le contrôle soit viable (pas garanti par pickControl)
                .applyLeft([&, this](pickedControl p) {
                    p.timeStep = min(p.timeStep, finalTime - time);
                    // structured binding
                    auto [pickedRho, cu] = p;
                    
                    if (!dynsys->isViableControl(xCoordsDouble, controlCoords[cu], imageVect, pickedRho)) {
                        spdlog::info("No viable real control found for strategy's pick.");
                        return OptionalCu(UNSATISFIED_STRATEGY);
                    }
                    else {                
                        return OptionalCu(p);
                    }
                });
            
            if (opt.isRight()) {
                succes = false;
                break;
            }
            pickedControl defaultC = {rho, 0};
            // Structured binding
            auto [pickedRho, cu] = opt.fromLeft(defaultC);
            rho = pickedRho;
            
            time += pickedRho;
            traj.addPoints(cu, realTimeStepsPerDiscreteStep, time, imageVect, flag);
            std::copy(imageVect, imageVect+dim, xCoordsDouble);
            currentDiscreteTrajPos = grid->getNearestPointInSet(xCoordsDouble);
            grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
            trajDiscrete.addPoint(doubleCoordsOnDiscreteTraj, time);
            
            ++nbIter;
        }

        delete[] xCoordsDouble;
        delete[] imageVect;
        delete[] doubleCoordsOnDiscreteTraj;
        delete[] currentPosIntCoords;
        
        string fileNameD(fileName);
        for (int k = 0; k < 4; k++)
        {
            fileNameD.pop_back();
        }
        fileNameD += "-Discrete.dat";
    
        if (succes)
        {
            spdlog::info("Trajectory reconstruction finished with success");
            spdlog::info("The real trajectory will be saved in the file {}, the discrete trajectory saved in the file {}", fileName, fileNameD);
        }
        else
        {
            spdlog::warn("Trajectory reconstruction failed before attempting the target set");
            spdlog::warn("The real partial trajectory will be saved in the file {}, the partial discrete trajectory saved in the file {}", fileName, fileNameD);
        }

        traj.writeToFile(fileName, picker.getStrategyNames());
        trajDiscrete.writeToFile(fileNameD);
    }
    
    return succes;
}

int ViabiBitSetTrajectoryHelper::computeTychasticStrategyTrajectory(double *initPosition, double finalTime, string fileName, TychePicker &tychePicker) {

    using ull = unsigned long long int;

    bool succes = false;

    SysDyn::PointStatus initStatus = dynsys->checkKernelRelation(initPosition);
    if (initStatus == SysDyn::OUTSIDE_DOMAIN) {
        spdlog::error("Initial point is out of the viability domain. Stop.");
    }
    else if (initStatus == SysDyn::OUTSIDE_CONSTRAINTS) {
        spdlog::error("Initial point is out of the constraints set. Stop.");
    }
    else if (initStatus == SysDyn::OUTSIDE_GRID) {
        spdlog::error("Initial point is out of the grid. Stop.");
    }
    else {
        TychasticControlPickerBitSet picker(this, dynsys, tpm);
        double **controlCoords = dynsys->getControlCoords();
        double **tycheCoords = dynsys->getTychCoords();
        
        double time = 0.0;
        int nbIter = 0;

        double *xCoordsDouble = new double[dim];
        std::copy(initPosition, initPosition+dim, xCoordsDouble);
        double *imageVect = new double[dim];
        double *doubleCoordsOnDiscreteTraj = new double[dim];        
        ull *currentPosIntCoords = new ull[dim];
                
        unsigned long long int currentDiscreteTrajPos = grid->getNearestPointInSet(xCoordsDouble);
        grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
        TychasticTrajectoryStorage traj = (tpm->getTrajectoryParameters()->SAVE_PICKING_STRATEGY)
            ? TychasticTrajectoryStorage::createFlagSavingStorage(initPosition, finalTime, dynsys)
            : TychasticTrajectoryStorage::createNonFlagSavingStorage(initPosition, finalTime, dynsys);
        
        TrajectoryPointsStorage trajDiscrete(doubleCoordsOnDiscreteTraj, finalTime, dynsys);        

        int realTimeStepsPerDiscreteStep = tpm->getTrajectoryParameters()->REAL_TIME_STEPS_PER_DISCRETE_STEP;
        succes = true;
        
        while (succes && time < finalTime && nbIter < NB_MAX_TRAJ_ITER) {

            ControlPicker::StrategyIndexBitFlag flag = 0;
            double rho = dynsys->calculRho_local(doubleCoordsOnDiscreteTraj);
            rho = min(rho, finalTime - time);
            
            unsigned long long int tyIndex = tychePicker.pickTyche(xCoordsDouble, time);
            if (tyIndex >= dynsys->getTotalNbPointsTy()) {
                succes = false;
                break;
            }

            picker.setTycheIndex(tyIndex);
            OptionalCu opt = picker.pickControl(traj, trajDiscrete, rho, flag)
                // On vérifie que le contrôle soit viable (pas garanti par pickControl)
                .applyLeft([&, this](pickedControl p) {
                    p.timeStep = min(p.timeStep, finalTime - time);
                    // structured binding
                    auto [pickedRho, cu] = p;
                    
                    if (!dynsys->isViableControl_tych(xCoordsDouble, controlCoords[cu], tycheCoords[tyIndex], imageVect, pickedRho)) {
                        spdlog::info("No viable real control found for strategy's pick.");
                        return OptionalCu(UNSATISFIED_STRATEGY);
                    }
                    else {
                        return OptionalCu(p);
                    }
                });
            
            if (opt.isRight()) {
                succes = false;
                break;
            }
            pickedControl defaultC = {rho, 0};
            // Structured binding
            auto [pickedRho, cu] = opt.fromLeft(defaultC);
            rho = pickedRho;


            dynsys->getTychasticImage(xCoordsDouble, controlCoords[cu], tycheCoords[tyIndex], imageVect, pickedRho);
            
            time += pickedRho;
            traj.addPoints(cu, tyIndex, realTimeStepsPerDiscreteStep, time, imageVect, flag);
            std::copy(imageVect, imageVect+dim, xCoordsDouble);
            currentDiscreteTrajPos = grid->getNearestPointInSet(xCoordsDouble);
            grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
            trajDiscrete.addPoint(doubleCoordsOnDiscreteTraj, time);
            
            ++nbIter;
        }

        delete[] xCoordsDouble;
        delete[] imageVect;
        delete[] doubleCoordsOnDiscreteTraj;
        delete[] currentPosIntCoords;
        
        string fileNameD(fileName);
        for (int k = 0; k < 4; k++)
        {
            fileNameD.pop_back();
        }
        fileNameD += "-Discrete.dat";
    
        if (succes)
        {
            spdlog::info("Trajectory reconstruction finished with success");
            spdlog::info("The real trajectory will be saved in the file {}, the discrete trajectory saved in the file {}", fileName, fileNameD);
        }
        else
        {
            spdlog::warn("Trajectory reconstruction failed before attempting the target set");
            spdlog::warn("The real partial trajectory will be saved in the file {}, the partial discrete trajectory saved in the file {}", fileName, fileNameD);
        }

        traj.writeToFile(fileName, picker.getStrategyNames());
        trajDiscrete.writeToFile(fileNameD);
    }
    
    return succes;
}
