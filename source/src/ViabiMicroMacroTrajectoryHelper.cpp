/*
 * ViabiMicroMacroTrajectoryHelperTrajectoryHelper.cpp
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
 *  Created on: 7 juil. 2024
 *      Author: Anna DESILLES
 */

#include <numeric> // Pour std::iota

#include "../include/ViabiMicroMacroTrajectoryHelper.h"
#include "../include/Params.h"

ViabiMicroMacroTrajectoryHelper::ViabiMicroMacroTrajectoryHelper()
{
    // TODO Auto-generated constructor stub

}

ViabiMicroMacroTrajectoryHelper::ViabiMicroMacroTrajectoryHelper(GridMicroMacro *gr, SysDyn *ds, TrajectoryParametersManager *tpm)
{

    const trajectoryParams *tp = tpm->getTrajectoryParameters();
    TypeTraj typeTraj = tp->TRAJECTORY_TYPE;
        
    grid = gr;
    dynsys = ds;
    dim = grid->dim;
    dimC = dynsys->getDimC();
    filePrefix = grid->filePrefix;
    this->typeTraj = typeTraj;
    vTab = grid->getGridPtr();

    unsigned long long int nbTotalC = dynsys->getTotalNbPointsC();    
    preferedControlIndexes = new int[nbTotalC];
    std::iota(preferedControlIndexes, preferedControlIndexes+nbTotalC, 0);
    controlWeight = tp->CONTROL_WEIGHT;
    sortIndexes = tp->SORT_INDEXES;
    
    switch (typeTraj)
	{
    case VD:
	{
        ViabiMicroMacroTrajectoryHelper::findfViableControl_DD = &ViabiMicroMacroTrajectoryHelper::findViabControlDefault_DD;
        break;
	}
    case VDI:
	{
        ViabiMicroMacroTrajectoryHelper::findfViableControl_DD = &ViabiMicroMacroTrajectoryHelper::findViabControlDiffControl_DD;
        break;
	}
    case VMM:
	{
        ViabiMicroMacroTrajectoryHelper::findfViableControl_DD = &ViabiMicroMacroTrajectoryHelper::findViabControlMinValue_DD;
        break;
	}
    default:
	{
        ViabiMicroMacroTrajectoryHelper::findfViableControl_DD = &ViabiMicroMacroTrajectoryHelper::findViabControlDefault_DD;
        break;
	}
	}
}

double ViabiMicroMacroTrajectoryHelper::computeOptimalTrajectory(double *initPosition, double T, string fileName, bool &succes)
    {
    /*
     * tableau de coordonnées  définissant
     * la grille de contrôles
     */
    double **controlCoords = dynsys->getControlCoords();

    /*
     * coordonnées du points corent de la trajectoire
     */
    double xCoordsDouble[dim], imageVect[dim], doubleCoordsOnDiscreteTraj[dim];
    unsigned long long int intCoordsOnDiscreteTraj[dim];
    /*
     * le pas de temps déterminé localement pour chaque point de la trajectoire
     */
    double rho;
    /*
     * listes  pour contenir la trajectoire ( temps-position) et les contrôles
     */
    list<valarray<double> > traj, trajC, trajDiscrete;
    list<double> realTrajBudget, valueAlongDiscreteTrajectory;
    /*
     * structures accumulables dansune liste
     */
    valarray<double> newTrajPoint(dim + 1);
    valarray<double> newTrajDiscretePoint(dim + 1);
    valarray<double> trajControlCoords(dimC);
    double currentBudget, nextBudget, valueAtCurrentDiscretePoint;

    spdlog::info("Start of optimal trajectory reconstruction");
    logVector("Initial point : ", initPosition, dim);
    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;
    double minCellVal = PLUS_INF;
    double time = 0.0;
    ostringstream os;
    string msg;

    if (grid->isPointInGrid(initPosition))
	{
	if (dynsys->constraintsX(initPosition) < PLUS_INF)
	    {

	    minCellVal = grid->getOptimalValue(initPosition);
	    testNonVide = (minCellVal < PLUS_INF);

	    if (!testNonVide)
		{
		spdlog::error("The initial point selected is not in the viable set");
		succes = 0;
		}
	    else
		{
		/*
		 * la position initiale se trouve dans le noyau de viabilité
		 * on initialise le temps à 0  et recopie la pos initiale
		 * dans le coordonnées temporaires du point en cours de la trajectoire
		 */
		unsigned long long int intDiscretPos = 0;

		time = 0.0;
		intDiscretPos = grid->getNearestPointInSet(initPosition);
		rho = dynsys->calculRho_local(initPosition);
		currentBudget = vTab[intDiscretPos];
		valueAtCurrentDiscretePoint = vTab[intDiscretPos];

		spdlog::info(" Value at initial point {}", currentBudget);
		int nbIter = 0;
		/*
		 * On itère tant que le temps n'a pas dépassé l'horizon donné
		 */

		double c, realTimeStep;
		int bestCu;
		unsigned long long int currentDiscreteTrajPos = intDiscretPos;
		unsigned long long int optimSuccessor;
		c = dynsys->target(initPosition);
		for (int l1 = 0; l1 < (int) dim; l1++)
		    {
		    xCoordsDouble[l1] = initPosition[l1];
		    }


		while ((time < T) && (c >= PLUS_INF) && (nbIter < NB_MAX_TRAJ_ITER))
		    {

		    os << "Time : " << time << ". Current point = ";
		    msg = os.str();
		    os.str("");
		    logVector(msg, xCoordsDouble, dim);

		    for (int i = 0; i < (int) dim; i++)
			{
			newTrajPoint[i] = xCoordsDouble[i];
			}
		    newTrajPoint[dim] = time;

		    traj.push_back(newTrajPoint);
		    realTrajBudget.push_back(currentBudget);
		    valueAlongDiscreteTrajectory.push_back(valueAtCurrentDiscretePoint);
		    grid->numToIntAndDoubleCoords(currentDiscreteTrajPos, intCoordsOnDiscreteTraj, doubleCoordsOnDiscreteTraj);

		    for (int i = 0; i < (int) dim; i++)
			{
			newTrajDiscretePoint[i] = doubleCoordsOnDiscreteTraj[i];
			}
		    newTrajDiscretePoint[dim] = time;

		    trajDiscrete.push_back(newTrajDiscretePoint);

		    rho = dynsys->calculRho_local(doubleCoordsOnDiscreteTraj);

		    rho = min(rho, T - time);
		    realTimeStep = rho;

		    optimSuccessor = this->findOptimalDiscreteSuccessor(currentDiscreteTrajPos, time, rho);
		    if (optimSuccessor > grid->nbTotalPoints)
			{
			spdlog::error("Optimal successor not found");
			break;
			}
		    else
			{
			realTimeStep = rho;
			bestCu = this->findOptiControl(xCoordsDouble, currentBudget, optimSuccessor, realTimeStep, 15, 0.6, imageVect, testNonVide, nextBudget);

			if (testNonVide)
			    {
			    currentDiscreteTrajPos = optimSuccessor;
			    valueAtCurrentDiscretePoint = vTab[currentDiscreteTrajPos];
			    currentBudget = nextBudget;
			    time += realTimeStep;
			    for (int dc = 0; dc < (int) dimC; dc++)
				{
				trajControlCoords[dc] = controlCoords[bestCu][dc];
				}
			    trajC.push_back(trajControlCoords);
			    for (int i = 0; i < (int) dim; i++)
				{
				xCoordsDouble[i] = imageVect[i];
				}
			    }
			else
			    {
			    spdlog::error("Impossible to get a control to keep real trajectory close to optimal discrete one");
			    break;
			    }
			}
		    c = dynsys->target(xCoordsDouble);
		    nbIter++;
		    }
		if (c < PLUS_INF)
		    {
		    succes = 1;
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


    string fileNameD(fileName);
    for(int k=0;k<4;k++)
	{
	fileNameD.pop_back();
	}
    fileNameD+="-Discrete.dat";

    if(succes)
	{
	spdlog::info("Trajectory reconstruction finished with success");
	spdlog::info("The real trajectory will be saved in the file {}, the discrete trajectory saved in the file {}", fileName, fileNameD);
	}
    else
	{
	spdlog::warn("Trajectory reconstruction failed before attempting the target set");
	spdlog::warn("The real partial trajectory will be saved in the file {}, the partial discrete trajectory saved in the file {}", fileName, fileNameD);
	}


    FILE *fi;
    fi = fopen(fileName.c_str(), "w");
    if (fi == NULL)
	{
	spdlog::error("Impossible to open the file {}", fileName);
	}
    else
	{
	list<valarray<double> >::iterator it = traj.begin();
	list<double>::iterator itRealBudget = realTrajBudget.begin();
	list<valarray<double> >::iterator itc = trajC.end();
	itc--;
	trajC.push_back(*itc);
	itc = trajC.begin();

	while (it != traj.end())
	    {

	    for (int l1 = 0; l1 < (int) dim; l1++)
		{
		fprintf(fi, "%15.8f ", (*it)[l1]);
		}
	    fprintf(fi, "%15.8f ", (*it)[dim]);
	    for (int dc = 0; dc < (int) dimC; dc++)
		{
		fprintf(fi, "%15.8f ", (*itc)[dc]);
		}
	    fprintf(fi, "%15.8f ", (*itRealBudget));
	    fprintf(fi, "\n");
	    it++;
	    itc++;
	    itRealBudget++;
	    }
	fclose(fi);
	}

    FILE *fiD;


    fiD = fopen(fileNameD.c_str(), "w");
    if (fiD == NULL)
	{
	spdlog::error("Impossible to open the file {}", fileNameD);
	}
    else
	{
	list<valarray<double> >::iterator itD = trajDiscrete.begin();
	list<double>::iterator itDiscretTrajBudget = valueAlongDiscreteTrajectory.begin();

	while (itD != trajDiscrete.end())
	    {

	    for (int l1 = 0; l1 < (int) dim; l1++)
		{
		fprintf(fiD, "%15.8f ", (*itD)[l1]);
		}
	    fprintf(fiD, "%15.8f ", (*itD)[dim]);
	    fprintf(fiD, "%15.8f ", (*itDiscretTrajBudget));
	    fprintf(fiD, "\n");
	    itD++;
	    itDiscretTrajBudget++;
	    }
	fclose(fiD);
	}

    return time;

    }

int ViabiMicroMacroTrajectoryHelper::findOptiControl(double *currentPos, double budget, unsigned long long int optimDiscreteSuccessor, double &dt, int nbStepIter,
	double stepCoeff, double *resPos, bool &succes, double &newBudget)
    {
    double realCoordsOnDiscreteTraj[dim];
    unsigned long long int intCoordsOnDiscreteTraj[dim];
    double rho0 = dt;
    int bestCu;

    grid->numToIntAndDoubleCoords(optimDiscreteSuccessor, intCoordsOnDiscreteTraj, realCoordsOnDiscreteTraj);

    double successorVal = vTab[optimDiscreteSuccessor];
    double imageVect[dim];


    double **controlCoords = dynsys->getControlCoords();
    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();

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

    spdlog::info("Start of computing of the real optimal control");

    for (int i = 0; i < (int) dim; i++)
	{
	xCoordsDouble[i] = currentPos[i];
	}
    os << "Initial budget  : " << budget << "Time step  : " << rho0 << ". Current point = ";
    msg = os.str();
    os.str("");
    logVector(msg, xCoordsDouble, dim);

    bestCu = 0;
    int iter = 0;
    double rhoOptim = rho0;
    double rhoMin = rho0 * (1 - stepCoeff);// we will test steps from rhoMin to  rhoMax = rho0 * (1 + stepCoeff);
    double rho = rhoMin;
    double step = rho0 * stepCoeff / nbStepIter;
    double optimBudget = budget;
    testNonVide = false;
    double minDist = PLUS_INF;
    while (iter < 2 * nbStepIter + 1)
	{
	unsigned long long int cu = 0;

	while (cu < nbCTotal)
	    {
	    /*
	     * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	     * au point en cours
	     */
	    if (dynsys->constraintsXU(xCoordsDouble, controlCoords[cu]) < PLUS_INF)
		{
		(dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[cu], imageVect, rho);
		if (grid->isPointInGrid(imageVect))
		    {
		    /*
		     *  le sucesseur est dans la grille de calcul
		     */
		    if (dynsys->constraintsX(imageVect) < PLUS_INF)
			{
			double tempL = dynsys->lFunc(xCoordsDouble, controlCoords[cu]);

			double tempL1 = dynsys->lFunc(imageVect, controlCoords[cu]);

			double newVal = budget - rho * 0.5 * (tempL + tempL1);

			testNonVide = true;
			double dist = 0.0;
			for (int k = 0; k < dim; k++)
			    {
			    dist += (imageVect[k] - realCoordsOnDiscreteTraj[k]) * (imageVect[k] - realCoordsOnDiscreteTraj[k]);
			    }
			dist = sqrt(dist);

			if (dist < minDist)
			    {
			    minDist = dist;
			    bestCu = cu;
			    rhoOptim = rho;
			    optimBudget = newVal;
			    }
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
	spdlog::info("Best real control found");
	spdlog::info("Discrete successor val {}, realSuccessor val {}, Distance {}, optimal time step {}",successorVal, optimBudget, minDist, rhoOptim );
	}
    succes = testNonVide;

    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[bestCu], imageVect, rhoOptim);
    for (int i = 0; i < (int) dim; i++)
	{
	resPos[i] = imageVect[i];
	}
    dt = rhoOptim;
    newBudget = optimBudget;
    return bestCu;
    }

unsigned long long int ViabiMicroMacroTrajectoryHelper::findOptimalDiscreteSuccessor(unsigned long long int pos, double time, double dt)
{
    int nbPointsCube = (int) pow(2.0, dim);
    unsigned long long int optimSuccessor = grid->nbTotalPoints + 1;
    unsigned long long int intCoordsOnDiscreteTraj[dim];

    /*
     * coordonnées du points corent de la trajectoire
     */
    double xCoordsDouble[dim];
    grid->numToIntAndDoubleCoords(pos, intCoordsOnDiscreteTraj, xCoordsDouble);

    double imageVect[dim];

    /*
     * tableau de coordonnées  définissant
     * la grille de contrôles
     */
    double **controlCoords = dynsys->getControlCoords();
    int *preferedControlIndexes = sortPreferedControlIndexes(xCoordsDouble, time);
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
     * le pas de temps déterminé localement pour chaque point de la trajectoire
     */
    double rho;
    /*
     * numéros de mailles
     */
    int cellNum = 0;

    int posTemp;

    /*
     * tests de validité de point initial
     */

    double currentVal = vTab[pos];
    rho = dt;
    /*
     * on parcours tous les contrôles
     */

    double minVal = PLUS_INF;
    unsigned long long int cPrefu = 0;

    minVal = PLUS_INF;
    while (cPrefu < nbCTotal)
	{
        /*
         * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
         * au point en cours
         */
        double *preferedControl = controlCoords[preferedControlIndexes[cPrefu]];
        
        if (dynsys->constraintsXU(xCoordsDouble, preferedControl) < PLUS_INF)
	    {
            (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, preferedControl, imageVect, rho);
            if (grid->isPointInGrid(imageVect))
            {
                /*
                 *  le sucesseur est dans la grille de calcul
                 */
                if (dynsys->constraintsX(imageVect) < PLUS_INF)
                {
                    /*
                     * le successeur vérifie les contraintes
                     * On identifie la maille où il se trouve
                     */
                    cellNum = grid->localizePoint(imageVect);
                    //cout<< " num cellule "<<cellNum<<endl;

                    /*
                     * On parcours les sommets de la maille
                     * autour du sucesseur pour voir s'il y a des
                     * points viables
                     */
                    double tempL = dynsys->lFunc(xCoordsDouble, preferedControl);

                    double tempL1 = dynsys->lFunc(imageVect, preferedControl);

                    double newVal = currentVal - rho * 0.5 * (tempL + tempL1);
                    for (int ii = 0; ii < nbPointsCube; ii++)
                    {
                        posTemp = cellNum + indicesDecalCell[ii];
                        if (newVal <= vTab[posTemp] + rho && vTab[posTemp] < minVal)
                        {
                            optimSuccessor = posTemp;
                            minVal = vTab[posTemp];
                        }
                    }
                }
            }
	    }
        cPrefu++;
	}
    if(optimSuccessor < grid->nbTotalPoints + 1)
	{
        spdlog::info("Optimal discrete sucessot found. Sucessorr value : {}, currentpoint value : {}", minVal, currentVal);
	}
    return optimSuccessor;
}

unsigned long long int ViabiMicroMacroTrajectoryHelper::findViabControlDefault_DD(double budget, unsigned long long int *currentPos,
	unsigned long long int currentControl, unsigned long long int *resPos, double &newBudget, bool &succes)
    {
    unsigned long long int **controlCoords = dynsys->getControlIntCoords();

    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    /*
     * coordonnées du points corent de la trajectoire
     */
    unsigned long long int xCoords[dim], intVect1[dim];

    unsigned long long int posTemp, imagePos;

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;

    for (int i = 0; i < (int) dim; i++)
	{
	xCoords[i] = currentPos[i];
	}

    int bestCu;

    /*
     * on parcours tous les contrôles
     */

    bestCu = 0;
    testNonVide = false;

    posTemp = grid->intCoordsToNum(xCoords);

    testNonVide = false;
    succes = false;

    unsigned long long int cu = 0;
    double tempM;
    bool eligibleControl;
    while (cu < nbCTotal && !testNonVide)
	{
	if (currentControl < nbCTotal)
	    {
	    eligibleControl = (dynsys->controlEligibilityForTraj_fd(xCoords, controlCoords[cu], controlCoords[currentControl]) < PLUS_INF);
	    }
	else
	    {
	    eligibleControl = true;
	    }

	if (eligibleControl)
	    {
	    if (dynsys->constraintsXU_fd(xCoords, controlCoords[cu]) < PLUS_INF)
		{
		dynsys->dynamics_fd(xCoords, controlCoords[cu], intVect1);

		if (grid->isPointInGrid_fd(intVect1))
		    {

		    /*!
		     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
		     */
		    if (dynsys->constraintsX_fd(intVect1) < PLUS_INF)
			{

			/*
			 * le successeur vérifie les contraintes
			 * On identifie la maille où il se trouve
			 */

			tempM = dynsys->muFunc_fd(xCoords, controlCoords[cu]);
			if (budget >= tempM)
			    {
			    double tempL = dynsys->lFunc_fd(xCoords, controlCoords[cu]);

			    double newVal = budget + tempL;

			    imagePos = grid->intCoordsToNum(intVect1);
			    testNonVide = (newVal >= vTab[imagePos]);
			    /*
			     * On parcours les sommets de la maille
			     * autour du sucesseur pour voir s'il y a des
			     * points viables
			     */
			    if (testNonVide)
				{
				newBudget = budget + tempL;
				bestCu = cu;
				for (int i = 0; i < (int) dim; i++)
				    {
				    resPos[i] = intVect1[i];
				    }
				}
			    }
			else
			    {
			    cout << " controle ne vérifie pas la contrainte de budget\n";
			    }

			}
		    }
		}
	    }
	cu++;
	}				//fin de parcours de tous les contrôles
    // la boucle s'arête ici u premier contrôle
    // qui donne un successeur viable
    succes = testNonVide;
    cout << " succes = " << succes << endl;
    return bestCu;
    }

unsigned long long int ViabiMicroMacroTrajectoryHelper::findViabControlDefault_tych_DD(double budget, unsigned long long int *currentPos,
	unsigned long long int currentControl, unsigned long long int currentTych, unsigned long long int *resPos, double &newBudget, bool &succes)
    {

    unsigned long long int **controlCoords = dynsys->getControlIntCoords();
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */

    unsigned long long int **tychIntCoords = dynsys->getTychIntCoords();

    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();

    /*
     * indices de déplacement pour parcourir les sommets d'une malle à partir
     * du sommet inf
     */
    /*
     * coordonnées du points corent de la trajectoire
     */
    unsigned long long int xCoords[dim], intVect1[dim];

    unsigned long long int posTemp, imagePos;

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;

    for (int i = 0; i < (int) dim; i++)
	{
	xCoords[i] = currentPos[i];
	}

    int bestCu;

    /*
     * on parcours tous les contrôles
     */

    bestCu = 0;
    testNonVide = false;

    posTemp = grid->intCoordsToNum(xCoords);

    testNonVide = false;
    succes = false;

    unsigned long long int cu = 0;
    double tempM;
    bool eligibleControl;
    while (cu < nbCTotal && !testNonVide)
	{
	/*
	 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	 * au point en cours
	 */

	if (currentControl < nbCTotal)
	    {
	    eligibleControl = (dynsys->controlEligibilityForTraj_fd(xCoords, controlCoords[cu], controlCoords[currentControl]) < PLUS_INF);
	    }
	else
	    {
	    eligibleControl = true;
	    }

	if (eligibleControl)
	    {
	    if (dynsys->constraintsXU_fd(xCoords, controlCoords[cu]) < PLUS_INF)
		{
		dynsys->dynamics_tych_fd(xCoords, controlCoords[cu], tychIntCoords[currentTych], intVect1);
		if (grid->isPointInGrid_fd(intVect1))
		    {
		    /*!
		     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
		     */
		    if (dynsys->constraintsX_fd(intVect1) < PLUS_INF)
			{
			/*
			 * le successeur vérifie les contraintes
			 * On identifie la maille où il se trouve
			 */

			tempM = dynsys->muFunc_fd(xCoords, controlCoords[cu]);
			if (budget >= tempM)
			    {
			    double tempL = dynsys->lFunc_tych_fd(xCoords, controlCoords[cu], tychIntCoords[currentTych]);

			    double newVal = budget + tempL;

			    imagePos = grid->intCoordsToNum(intVect1);
			    testNonVide = (newVal >= vTab[imagePos]);
			    /*
			     * On parcours les sommets de la maille
			     * autour du sucesseur pour voir s'il y a des
			     * points viables
			     */
			    if (testNonVide)
				{
				newBudget = budget + tempL;
				bestCu = cu;
				for (int i = 0; i < (int) dim; i++)
				    {
				    resPos[i] = intVect1[i];
				    }
				}
			    }
			else
			    {
			    cout << " controle ne vérifie pas la contrainte de budget\n";
			    }

			}
		    }
		}
	    }
	cu++;
	}				//fin de parcours de tous les contrôles
    // la boucle s'arête ici u premier contrôle
    // qui donne un successeur viable
    succes = testNonVide;
    cout << " succes = " << succes << endl;
    return bestCu;
    }

unsigned long long int ViabiMicroMacroTrajectoryHelper::findViabControlMinValue_DD(double budget, unsigned long long int *currentPos,
	unsigned long long int currentControl, unsigned long long int *resPos, double &newBudget, bool &succes)
    {

    unsigned long long int **controlCoords = dynsys->getControlIntCoords();

    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    /*
     * coordonnées du points corent de la trajectoire
     */
    unsigned long long int xCoords[dim], intVect1[dim];

    unsigned long long int posTemp, imagePos;

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;

    for (int i = 0; i < (int) dim; i++)
	{
	xCoords[i] = currentPos[i];
	}

    int bestCu;

    /*
     * on parcours tous les contrôles
     */

    bestCu = 0;
    testNonVide = false;
    double minVal = PLUS_INF;
    posTemp = grid->intCoordsToNum(xCoords);

    testNonVide = false;
    succes = false;

    unsigned long long int cu = 0;
    while (cu < nbCTotal)
	{
	/*
	 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	 * au point en cours
	 */
	if (dynsys->controlEligibilityForTraj_fd(xCoords, controlCoords[cu], controlCoords[currentControl]) < PLUS_INF)
	    {
	    if (dynsys->constraintsXU_fd(xCoords, controlCoords[cu]) < PLUS_INF)
		{
		dynsys->dynamics_fd(xCoords, controlCoords[cu], intVect1);
		if (grid->isPointInGrid_fd(intVect1))
		    {
		    /*!
		     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
		     */
		    if (dynsys->constraintsX_fd(intVect1) < PLUS_INF)
			{
			/*
			 * le successeur vérifie les contraintes
			 * On identifie la maille où il se trouve
			 */
			double tempM = dynsys->muFunc_fd(xCoords, controlCoords[cu]);
			if (budget >= tempM)
			    {
			    double tempL = dynsys->lFunc_fd(xCoords, controlCoords[cu]);

			    double newVal = budget + tempL;

			    imagePos = grid->intCoordsToNum(intVect1);
			    testNonVide = (newVal >= vTab[imagePos]);
			    /*
			     * On parcours les sommets de la maille
			     * autour du sucesseur pour voir s'il y a des
			     * points viables
			     */
			    if (testNonVide)
				{
				succes = testNonVide;
				if (vTab[imagePos] < minVal)
				    {
				    minVal = vTab[imagePos];
				    newBudget = budget + tempL;
				    bestCu = cu;
				    for (int i = 0; i < (int) dim; i++)
					{
					resPos[i] = intVect1[i];
					}
				    }
				}
			    }
			else
			    {
			    cout << " controle ne vérifie pas la contrainte de budget\n";
			    }
			}
		    }
		}
	    }
	cu++;
	}				//fin de parcours de tous les contrôles
    // la boucle s'arête ici u premier contrôle
    // qui donne un successeur viable

    cout << " succes = " << succes << endl;
    return bestCu;
    }

unsigned long long int ViabiMicroMacroTrajectoryHelper::findViabControlDiffControl_DD(double budget, unsigned long long int *currentPos,
	unsigned long long int currentControl, unsigned long long int *resPos, double &newBudget, bool &succes)
    {

    unsigned long long int **controlCoords = dynsys->getControlIntCoords();

    /*
     * nombre total de points de contrôle
     */
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    /*
     * coordonnées du points corent de la trajectoire
     */
    unsigned long long int xCoords[dim], intVect1[dim];

    unsigned long long int posTemp, imagePos;

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;

    for (int i = 0; i < (int) dim; i++)
	{
	xCoords[i] = currentPos[i];
	}

    int bestCu;

    /*
     * on parcours tous les contrôles
     */

    bestCu = 0;
    testNonVide = false;

    posTemp = grid->intCoordsToNum(xCoords);

    testNonVide = false;
    succes = false;

    unsigned long long int cu = 0;
    while (cu < nbCTotal)
	{
	/*
	 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	 * au point en cours
	 */
	if (dynsys->controlEligibilityForTraj_fd(xCoords, controlCoords[cu], controlCoords[currentControl]) < PLUS_INF)
	    {
	    if (dynsys->constraintsXU_fd(xCoords, controlCoords[cu]) < PLUS_INF)
		{
		dynsys->dynamics_fd(xCoords, controlCoords[cu], intVect1);
		if (grid->isPointInGrid_fd(intVect1))
		    {
		    /*!
		     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
		     */
		    if (dynsys->constraintsX_fd(intVect1) < PLUS_INF)
			{
			/*
			 * le successeur vérifie les contraintes
			 * On identifie la maille où il se trouve
			 */
			double tempM = dynsys->muFunc_fd(xCoords, controlCoords[cu]);
			if (budget >= tempM)
			    {
			    bool testEqualControls = (cu == currentControl);

			    double tempL = dynsys->lFunc_fd(xCoords, controlCoords[cu]);

			    double newVal = budget + tempL;

			    imagePos = grid->intCoordsToNum(intVect1);
			    if (!testNonVide && newVal >= vTab[imagePos])
				{
				cout << " Controle viable trouve " << cu << endl;
				testNonVide = true;
				succes = testNonVide;
				newBudget = budget + tempL;
				bestCu = cu;
				for (int i = 0; i < (int) dim; i++)
				    {
				    resPos[i] = intVect1[i];
				    }
				if (!testEqualControls)
				    {
				    cout << " on arrete c'est le bon \n";
				    break;
				    }
				}
			    else if (newVal >= vTab[imagePos])
				{
				succes = true;

				if (!testEqualControls)
				    {
				    cout << " controle viable different trouve :new budget = " << newBudget << " image val = " << vTab[imagePos]
																       << " testNonVide= " << testNonVide << " currentCu " << currentControl << endl;

				    newBudget = budget + tempL;
				    bestCu = cu;
				    for (int i = 0; i < (int) dim; i++)
					{
					resPos[i] = intVect1[i];
					}
				    break;
				    }
				}
			    }
			else
			    {
			    cout << " controle ne vérifie pas la contrainte de budget\n";
			    }
			}
		    }
		}
	    }
	cu++;
	}				//fin de parcours de tous les contrôles
    // la boucle s'arête ici u premier contrôle
    // qui donne un successeur viable

    cout << " on retourne bestCu = " << bestCu << " succes = " << succes << endl;
    return bestCu;
    }

double ViabiMicroMacroTrajectoryHelper::computeViableTrajectory_DD(unsigned long long int *initPosition, double initValue, string fileName,
	bool &succes)
    {
    /*
     * tableau de coordonnées  définissant
     * la grille de contrôles
     */
    unsigned long long int **controlCoords = dynsys->getControlIntCoords();
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();

    /*
     * coordonnées du points corent de la trajectoire
     */

    unsigned long long int xCoordsInt[dim], imageVect[dim];
    /*
     * listes  pour contenir la trajectoire ( temps-position) et les contrôles
     */
    list<valarray<unsigned long long int> > traj;
    list<double> valsOptiTraj;
    list<double> valsRealTraj;
    list<valarray<unsigned long long int> > trajC;

    valarray<unsigned long long int> newTrajPoint(dim);
    /*
     * structures accumulables dansune liste
     */

    valarray<unsigned long long int> trajControlCoords(dimC);

    unsigned long long int posTemp;

    cout << " calcul de traj a partir de coords \n";
    cout << " Postion initiale = ";

    for (int l1 = 0; l1 < (int) dim; l1++)
	{
	cout << " " << initPosition[l1];
	}
    cout << " \n";
    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;
    double budget, newBudget;
    unsigned long long int bestCu, currentCu = nbCTotal + 1;

    if (grid->isPointInGrid_fd(initPosition))
	{
	if (dynsys->constraintsX_fd(initPosition) < PLUS_INF)
	    {
            posTemp = grid->intCoordsToNum(initPosition);

	    cout << " pos temp num = " << posTemp << endl;
	    testNonVide = (vTab[posTemp] < PLUS_INF);

	    if (!testNonVide)
		{
		cout << " La position initiale sélectionnée n'appartiant pas au noyau de viabilité\n";
		succes = 0;
		}
	    else
		{
		/*
		 * la position initiale se trouve dans le noyau de viabilité
		 * on initialise le temps à 0  et recopie la pos initiale
		 * dans le coordonnées temporaires du point en cours de la trajectoire
		 */
		if (initValue < vTab[posTemp])
		    {
		    cout
		    << "  le budget initial est inf�rieur au budget minimum requis pour cet �tat initial. Construction de trajectoire � partir du budget minimum requis\n";
		    initValue = vTab[posTemp];
		    }

		budget = initValue;
		cout << " value of init point  " << vTab[posTemp] << endl;
		cout << " budget initial  " << budget << endl;

		for (int i = 0; i < (int) dim; i++)
		    {
		    xCoordsInt[i] = initPosition[i];
		    }
		int nbIter = 0;
		/*
		 * On itère tant que le temps n'a pas dépassé l'horizon donné
		 */

		double c;
		c = dynsys->target_fd(xCoordsInt);
		cout << " cible donne au point init " << c << endl;
		while ((c >= PLUS_INF) && (testNonVide) && (nbIter < NB_MAX_TRAJ_ITER))
		    {
		    cout << " point en cours ";
		    for (int i = 0; i < (int) dim; i++)
			{
			newTrajPoint[i] = xCoordsInt[i];
			cout << " " << newTrajPoint[i];
			}
		    posTemp = grid->intCoordsToNum(xCoordsInt);
		    cout << " pos temp num = " << posTemp << endl;
		    cout << " value of   point  " << vTab[posTemp] << endl;
		    cout << " budget initial  " << budget << endl;
		    traj.push_back(newTrajPoint);
		    valsOptiTraj.push_back(vTab[posTemp]);
		    valsRealTraj.push_back(budget);

		    bestCu = (this->*findfViableControl_DD)(budget, xCoordsInt, currentCu, imageVect, newBudget, testNonVide);
		    currentCu = bestCu;
		    // la boucle s'arête ici u premier contrôle
		    // qui donne un successeur viable

		    // contrôle viable trouvé
		    // on recopie ce contrôle dans la liste et
		    // le successeur devient le point  courent
		    if (testNonVide)
			{

			cout << " image interieure tourvee \n";
			cout << "controle viable trouve : best Cu = " << bestCu << " currentCu = " << currentCu << " ";
			for (int dc = 0; dc < (int) dimC; dc++)
			    {
			    trajControlCoords[dc] = controlCoords[bestCu][dc];
			    cout << " " << trajControlCoords[dc];
			    }
			cout << endl;
			trajC.push_back(trajControlCoords);
			for (int i = 0; i < (int) dim; i++)
			    {
			    xCoordsInt[i] = imageVect[i];
			    }
			budget = newBudget;
			cout << " new budget = " << newBudget << endl;
			}
		    else
			{

			cout << "   Echec! Sortie de l'ensemble viable \n";
			break;
			}
		    c = dynsys->target_fd(xCoordsInt);
		    cout << "   valeur cible : " << c << endl;
		    nbIter++;
		    }
		if (c < PLUS_INF)
		    {
		    succes = 1;
		    cout << " dernier point en cours ";
		    for (int i = 0; i < (int) dim; i++)
			{
			newTrajPoint[i] = xCoordsInt[i];
			cout << " " << newTrajPoint[i];
			}
		    posTemp = grid->intCoordsToNum(xCoordsInt);
		    traj.push_back(newTrajPoint);
		    valsOptiTraj.push_back(vTab[posTemp]);
		    valsRealTraj.push_back(budget);
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
	printf("** error: impossible to open the file %s.\n", fileName.c_str());
	}
    else
	{
	list<valarray<unsigned long long int> >::iterator it = traj.begin();
	list<valarray<unsigned long long int> >::iterator itc = trajC.end();
	list<double>::iterator itVO = valsOptiTraj.begin();
	list<double>::iterator itVR = valsRealTraj.begin();
	itc--;
	trajC.push_back(*itc);
	itc = trajC.begin();

	while (it != traj.end())
	    {

	    for (int l1 = 0; l1 < dim; l1++)
		{
		fprintf(fi, "%15llu ", (*it)[l1]);
		}
	    for (int dc = 0; dc < dimC; dc++)
		{
		fprintf(fi, "%15llu ", (*itc)[dc]);
		}
	    fprintf(fi, "%15.8f ", (*itVR));
	    fprintf(fi, "%15.8f ", (*itVO));
	    fprintf(fi, "\n");
	    it++;
	    itc++;
	    itVR++;
	    itVO++;
	    }
	fclose(fi);
	}
    }
    return newBudget;
    }

double ViabiMicroMacroTrajectoryHelper::computeViableTrajectory_tych_DD(unsigned long long int *initPosition, double initValue, string fileName,
	bool &succes)
    {
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */
    // unsigned long long int   testI[dim];
    // double testV[dim];
    /*
     * tableau de coordonnées  définissant
     * la grille de contrôles
     */
    unsigned long long int **controlCoords = dynsys->getControlIntCoords();
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int **tychIntCoords = dynsys->getTychIntCoords();
    unsigned long long int nbTy = dynsys->getTotalNbPointsTy();

    int dimTy = dynsys->getDimTy();

    /*
     * coordonnées du points corent de la trajectoire
     */

    unsigned long long int xCoordsInt[dim], imageVect[dim];
    /*
     * listes  pour contenir la trajectoire ( temps-position) et les contrôles
     */
    list<valarray<unsigned long long int> > traj;
    list<double> valsOptiTraj;
    list<double> valsRealTraj;
    list<valarray<unsigned long long int> > trajC;
    list<valarray<unsigned long long int> > trajTy;
    valarray<unsigned long long int> newTrajPoint(dim);
    /*
     * structures accumulables dansune liste
     */

    valarray<unsigned long long int> trajControlCoords(dimC);
    valarray<unsigned long long int> trajTyControlCoords(dimTy);

    unsigned long long int posTemp;

    cout << " calcul de traj a partir de coords \n";
    cout << " Postion initiale = ";

    for (int l1 = 0; l1 < (int) dim; l1++)
	{
	cout << " " << initPosition[l1];
	}
    cout << " \n";
    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;
    double budget, newBudget;
    unsigned long long int bestCu, currentCu = nbCTotal + 1;

    if (grid->isPointInGrid_fd(initPosition))
	{
	if (dynsys->constraintsX_fd(initPosition) < PLUS_INF)
	    {
            posTemp = grid->intCoordsToNum(initPosition);

	    cout << " pos temp num = " << posTemp << endl;
	    testNonVide = (vTab[posTemp] < PLUS_INF);

	    if (!testNonVide)
		{
		cout << " La position initiale sélectionnée n'appartiant pas au noyau de viabilité\n";
		succes = 0;
		}
	    else
		{
		/*
		 * la position initiale se trouve dans le noyau de viabilité
		 * on initialise le temps à 0  et recopie la pos initiale
		 * dans le coordonnées temporaires du point en cours de la trajectoire
		 */
		if (initValue < vTab[posTemp])
		    {
		    cout
		    << "  le budget initial est inf�rieur au budget minimum requis pour cet �tat initial. Construction de trajectoire � partir du budget minimum requis\n";
		    initValue = vTab[posTemp];
		    }

		budget = initValue;
		cout << " value of init point  " << vTab[posTemp] << endl;
		cout << " budget initial  " << budget << endl;

		for (int i = 0; i < (int) dim; i++)
		    {
		    xCoordsInt[i] = initPosition[i];
		    }
		int nbIter = 0;
		/*
		 * On itère tant que le temps n'a pas dépassé l'horizon donné
		 */

		double c;
		c = dynsys->target_fd(xCoordsInt);
		cout << " cible donne au point init " << c << endl;
		while ((c >= PLUS_INF) && (testNonVide) && (nbIter < NB_MAX_TRAJ_ITER))
		    {
		    cout << " point en cours ";
		    for (int i = 0; i < (int) dim; i++)
			{
			newTrajPoint[i] = xCoordsInt[i];
			cout << " " << newTrajPoint[i];
			}
		    posTemp = grid->intCoordsToNum(xCoordsInt);
		    cout << " pos temp num = " << posTemp << endl;
		    cout << " value of   point  " << vTab[posTemp] << endl;
		    cout << " budget initial  " << budget << endl;
		    traj.push_back(newTrajPoint);
		    valsOptiTraj.push_back(vTab[posTemp]);
		    valsRealTraj.push_back(budget);

		    unsigned long long int currentTych = (unsigned long long int) (rand() % nbTy);
		    bestCu = findViabControlDefault_tych_DD(budget, xCoordsInt, currentCu, currentTych, imageVect, newBudget, testNonVide);
		    currentCu = bestCu;
		    // la boucle s'arête ici u premier contrôle
		    // qui donne un successeur viable

		    //cout<<   "Recherche de controle viable  fini parcours de controles on a  test non vide "<<testNonVide<< " bes c u= "<<bestCu<<endl;

		    // contrôle viable trouvé
		    // on recopie ce contrôle dans la liste et
		    // le successeur devient le point  courent
		    if (testNonVide)
			{

			cout << " image interieure tourvee \n";
			cout << "controle viable trouve : best Cu = " << bestCu << " currentCu = " << currentCu << " ";
			for (int dc = 0; dc < (int) dimC; dc++)
			    {
			    trajControlCoords[dc] = controlCoords[bestCu][dc];
			    cout << " " << trajControlCoords[dc];
			    }
			cout << endl;
			trajC.push_back(trajControlCoords);
			cout << "controle Ty :   = " << currentTych << " ";
			for (int dc = 0; dc < (int) dimTy; dc++)
			    {
			    trajTyControlCoords[dc] = tychIntCoords[currentTych][dc];
			    cout << " " << trajTyControlCoords[dc];
			    }
			cout << endl;
			trajTy.push_back(trajTyControlCoords);
			for (int i = 0; i < (int) dim; i++)
			    {
			    xCoordsInt[i] = imageVect[i];
			    }
			budget = newBudget;
			cout << " new budget = " << newBudget << endl;
			}
		    else
			{

			cout << "   Echec! Sortie de l'ensemble viable \n";
			break;
			}
		    c = dynsys->target_fd(xCoordsInt);
		    cout << "   valeur cible : " << c << endl;
		    nbIter++;
		    }
		if (c < PLUS_INF)
		    {
		    succes = 1;
		    cout << " dernier point en cours ";
		    for (int i = 0; i < (int) dim; i++)
			{
			newTrajPoint[i] = xCoordsInt[i];
			cout << " " << newTrajPoint[i];
			}
		    posTemp = grid->intCoordsToNum(xCoordsInt);
		    traj.push_back(newTrajPoint);
		    valsOptiTraj.push_back(vTab[posTemp]);
		    valsRealTraj.push_back(budget);
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
	printf("** error: impossible to open the file %s.\n", fileName.c_str());
	}
    else
	{
	list<valarray<unsigned long long int> >::iterator it = traj.begin();
	list<valarray<unsigned long long int> >::iterator itc = trajC.end();
	list<valarray<unsigned long long int> >::iterator itTy = trajTy.end();
	list<double>::iterator itVO = valsOptiTraj.begin();
	list<double>::iterator itVR = valsRealTraj.begin();
	itc--;
	trajC.push_back(*itc);
	itc = trajC.begin();

	itTy--;
	trajTy.push_back(*itTy);
	itTy = trajTy.begin();

	while (it != traj.end())
	    {

	    for (int l1 = 0; l1 < dim; l1++)
		{
		fprintf(fi, "%15llu ", (*it)[l1]);
		}
	    for (int dc = 0; dc < dimC; dc++)
		{
		fprintf(fi, "%15llu ", (*itc)[dc]);
		}
	    for (int dc = 0; dc < dimTy; dc++)
		{
		fprintf(fi, "%15llu ", (*itTy)[dc]);
		}
	    fprintf(fi, "%15.8f ", (*itVR));
	    fprintf(fi, "%15.8f ", (*itVO));
	    fprintf(fi, "\n");
	    it++;
	    itc++;
	    itTy++;
	    itVR++;
	    itVO++;
	    }
	fclose(fi);
	}
    }
    return newBudget;
    }

double ViabiMicroMacroTrajectoryHelper::computeViableTrajectory(double *initPosition, double initValue, string fileName, bool &succes)
    {
    succes = false;
    return 0.0;
    }

int *ViabiMicroMacroTrajectoryHelper::getPreferedControlIndexes() {
    return preferedControlIndexes;
}

int *ViabiMicroMacroTrajectoryHelper::sortPreferedControlIndexes(const double *x, double time, int strategyIndex) {
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    double **controlCoords = dynsys->getControlCoords();
    sortIndexes(preferedControlIndexes, x, controlCoords, nbCTotal, time, controlWeight, trajIndex, strategyIndex);
    return preferedControlIndexes;
}

ViabiMicroMacroTrajectoryHelper::~ViabiMicroMacroTrajectoryHelper()
    {
    // TODO Auto-generated destructor stub
        delete [] preferedControlIndexes;
    }

