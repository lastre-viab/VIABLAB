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

#include "ViabiMicroMacroTrajectoryHelper.h"

ViabiMicroMacroTrajectoryHelper::ViabiMicroMacroTrajectoryHelper()
    {
    // TODO Auto-generated constructor stub

    }

ViabiMicroMacroTrajectoryHelper::ViabiMicroMacroTrajectoryHelper(
	GridMicroMacro *gr, SysDyn *ds, int type)
    {
    grid = gr;
    dynsys = ds;
    dim = grid->dim;
    dimC = dynsys->getDimC();
    filePrefix = grid->filePrefix;
    typeTraj = type;
    vTab = grid->getGridPtr();
    switch (typeTraj)
	{
    case VD:
	{
	ViabiMicroMacroTrajectoryHelper::findfViableControl_DD =
		&ViabiMicroMacroTrajectoryHelper::findViabControlDefault_DD;
	break;
	}
    case VDI:
	{
	ViabiMicroMacroTrajectoryHelper::findfViableControl_DD =
		&ViabiMicroMacroTrajectoryHelper::findViabControlDiffControl_DD;
	break;
	}
    case VMM:
	{
	ViabiMicroMacroTrajectoryHelper::findfViableControl_DD =
		&ViabiMicroMacroTrajectoryHelper::findViabControlMinValue_DD;
	break;
	}
    default:
	{
	ViabiMicroMacroTrajectoryHelper::findfViableControl_DD =
		&ViabiMicroMacroTrajectoryHelper::findViabControlDefault_DD;
	break;
	}
	}
    }

double ViabiMicroMacroTrajectoryHelper::computeOptimalTrajectory_Lmin(
	double *initPosition, string fileName, bool &succes)
    {
    int nbPointsCube = (int) pow(2.0, dim); //number of vertexes of a cell

    /*
     * coordinates of all points in the control grid
     */
    double **controlCoords = dynsys->getControlCoords();

    /*
     * coordinates of shifts of all vertexes of cell relative to the inf corner
     */
    long long int *indicesDecalCell = grid->getIndicesDecalCell();
    /*
     * cooordinates of the current point of trajectory
     */
    double xCoordsDouble[dim], imageVect[dim]; // image vector

    double rho; //time step of the trajectory

    /*
     * lists to register the points of teh trajectory and of the corresponding controls
     */
    list<valarray<double> > traj, trajC;
    double T = dynsys->getTimeHorizon();

    valarray<double> newTrajPoint(dim + 1);
    valarray<double> trajControlCoords(dimC);

    int posTemp;
    int cellNum;

    cout << "Starting the integration of a trajectory. Initial point :  \n";
    printVector(initPosition, dim);
    /*
     * First, we check if the initial point is eligible for trajectory computation
     */
    bool testNonVide = false;

    double minCellVal = PLUS_INF;

    double time = 0.0;

    if (grid->isPointInGrid(initPosition))
	{
	if (dynsys->constraintsX(initPosition) < PLUS_INF)
	    {
	    cellNum = grid->localizePoint(initPosition);
	    minCellVal = -PLUS_INF;
	    for (int ii = 0; ii < nbPointsCube; ii++)
		{
		posTemp = cellNum + indicesDecalCell[ii];
		minCellVal = max(minCellVal, vTab[posTemp]);
		}
	    double budget = minCellVal;

	    cout << "Value of initial point  " << budget << endl;

	    testNonVide = (budget < PLUS_INF);

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
		double newBudget;
		time = 0.0;
		for (int i = 0; i < (int) dim; i++)
		    {
		    xCoordsDouble[i] = initPosition[i];
		    }
		int nbIter = 0;
		/*
		 * On itère tant que le temps n'a pas dépassé l'horizon donné
		 */

		double c;
		int bestCu;
		c = dynsys->target(xCoordsDouble);
		while ((time < T) && (c >= PLUS_INF)
			&& (nbIter <= NB_MAX_TRAJ_ITER))
		    {
		    cout << " point en cours ";
		    for (int i = 0; i < (int) dim; i++)
			{
			newTrajPoint[i] = xCoordsDouble[i];
			cout << " " << newTrajPoint[i];
			}
		    newTrajPoint[dim] = time;
		    cout << " temps= " << newTrajPoint[dim] << endl;
		    traj.push_back(newTrajPoint);

		    rho = 0.75 * dynsys->calculRho_local(xCoordsDouble);

		    rho = min(rho, T - time);

		    cout << " rho= " << rho << endl;

		    bestCu = this->findOptiControl_Lmin(budget, xCoordsDouble,
			    rho, 1, 1.0, imageVect, newBudget, testNonVide);
		    time += rho;

		    // la boucle s'arête ici u premier contrôle
		    // qui donne un successeur viable

		    cout
			    << " Premiere recherche de controle viable  fini parcours de controles on a  test non vide "
			    << testNonVide << " bes c u= " << bestCu << endl;

		    // contrôle viable trouvé
		    // on recopie ce contrôle dans la liste et
		    // le successeur devient le point  courent
		    if (testNonVide)
			{

			cout << " image interieure tourvee \n";

			for (int dc = 0; dc < (int) dimC; dc++)
			    {
			    trajControlCoords[dc] = controlCoords[bestCu][dc];
			    }
			trajC.push_back(trajControlCoords);
			for (int i = 0; i < (int) dim; i++)
			    {
			    xCoordsDouble[i] = imageVect[i];
			    }
			budget = newBudget;
			cout << " new budget = " << newBudget << endl;
			}
		    else
			{
			cout
				<< "  recherche de optimal viable avec iterations  sur le pas  de temps\n";
			rho = 0.5 * dynsys->calculRho_local(xCoordsDouble);

			rho = min(rho, T - time);

			bestCu = this->findOptiControl_Lmin(budget,
				xCoordsDouble, rho, 5, 1.1, imageVect,
				newBudget, testNonVide);
			if (testNonVide)
			    {
			    for (int dc = 0; dc < dimC; dc++)
				{
				trajControlCoords[dc] =
					controlCoords[bestCu][dc];
				}
			    trajC.push_back(trajControlCoords);
			    cout << "  coords double : ";
			    for (int i = 0; i < dim; i++)
				{
				xCoordsDouble[i] = imageVect[i];
				cout << " " << xCoordsDouble[i];
				}
			    cout << endl;
			    budget = newBudget;
			    cout << " new budget = " << newBudget << endl;
			    }

			else
			    {
			    cout << "   Echec! Sortie de l'ensemble viable \n";
			    break;

			    }
			}
		    c = dynsys->target(xCoordsDouble);
		    cout << "   valeur cible : " << c << endl;
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

    return time;

    }

double ViabiMicroMacroTrajectoryHelper::computeOptimalTrajectory(
	double *initPosition, string fileName, bool &succes)
    {

    int nbPointsCube = (int) pow(2.0, dim);			//pow(2.0, dim);

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
    double T = dynsys->getTimeHorizon();
    /*
     * structures accumulables dansune liste
     */
    valarray<double> newTrajPoint(dim + 1);
    valarray<double> trajControlCoords(dimC);

    int posTemp;

    cout << " calcul de traj TMIN a partir de coords \n";
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
    double minCellVal = PLUS_INF;
    double time = 0.0;

    if (grid->isPointInGrid(initPosition))
	{
	if (dynsys->constraintsX(initPosition) < PLUS_INF)
	    {

	    minCellVal = grid->getOptimalValue(initPosition);
	    cout << " val of initial point is " << minCellVal << endl;
	    testNonVide = (minCellVal < PLUS_INF);

	    if (!testNonVide)
		{
		cout
			<< "The selected initial point is not in the capture basin\n";
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
		for (int i = 0; i < (int) dim; i++)
		    {
		    xCoordsDouble[i] = initPosition[i];
		    }
		//intDiscretPos = grid->getNearestPointInSet(xCoordsDouble);
		intDiscretPos = grid->getBestNearPointInSet(xCoordsDouble);

		int nbIter = 0;
		/*
		 * On itère tant que le temps n'a pas dépassé l'horizon donné
		 */

		double c, realTimeStep;
		int bestCu;
		unsigned long long int currentDiscreteTrajPos = intDiscretPos;
		unsigned long long int optimSuccessor;
		c = dynsys->target(xCoordsDouble);
		while ((time < T) && (c >= PLUS_INF)
			&& (nbIter <= NB_MAX_TRAJ_ITER))
		    {
		    cout << " point en cours ";
		    for (int i = 0; i < (int) dim; i++)
			{
			newTrajPoint[i] = xCoordsDouble[i];
			cout << " " << newTrajPoint[i];
			}
		    newTrajPoint[dim] = time;
		    cout << " temps= " << newTrajPoint[dim] << endl;
		    traj.push_back(newTrajPoint);

		    rho = dynsys->calculRho_local(xCoordsDouble);

		    rho = min(rho, T - time);
		    realTimeStep = rho;
		    cout << " rho= " << rho << endl;

		    optimSuccessor = this->findOptimalDiscreteSuccessor(
			    currentDiscreteTrajPos, rho);
		    if (optimSuccessor > grid->nbTotalPoints)
			{
			cout << "   Echec! Sortie de l'ensemble viable \n";
			break;
			}
		    else
			{
			realTimeStep = 0.95 * rho;
			bestCu = this->findOptiControl(xCoordsDouble,
				optimSuccessor, realTimeStep, 5, 1.02,
				imageVect, testNonVide);

			// la boucle s'arête ici u premier contrôle
			// qui donne un successeur viable

			if (testNonVide)
			    {
			    // contrôle viable trouvé
			    // on recopie ce contrôle dans la liste et
			    // le successeur devient le point  courent
			    currentDiscreteTrajPos = optimSuccessor;
			    time += realTimeStep;
			    for (int dc = 0; dc < (int) dimC; dc++)
				{
				trajControlCoords[dc] =
					controlCoords[bestCu][dc];
				}
			    trajC.push_back(trajControlCoords);
			    for (int i = 0; i < (int) dim; i++)
				{
				xCoordsDouble[i] = imageVect[i];
				}

			    }
			else
			    {
			    realTimeStep = 1.25 * rho;
			    bestCu = this->findOptiControl(xCoordsDouble,
				    optimSuccessor, realTimeStep, 10, 0.9,
				    imageVect, testNonVide);

			    if (testNonVide)
				{
				// contrôle viable trouvé
				// on recopie ce contrôle dans la liste et
				// le successeur devient le point  courent
				currentDiscreteTrajPos = optimSuccessor;
				time += realTimeStep;
				for (int dc = 0; dc < (int) dimC; dc++)
				    {
				    trajControlCoords[dc] =
					    controlCoords[bestCu][dc];
				    }
				trajC.push_back(trajControlCoords);
				for (int i = 0; i < (int) dim; i++)
				    {
				    xCoordsDouble[i] = imageVect[i];
				    }

				}
			    else
				{
				cout
					<< "   Echec! Sortie de l'ensemble viable \n";
				break;

				}
			    }
			}
		    c = dynsys->target(xCoordsDouble);
		    cout << "   valeur cible : " << c << endl;
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

		for (int l1 = 0; l1 < (int) dim; l1++)
		    {
		    fprintf(fi, "%15.8f ", (*it)[l1]);
		    //		cout<< " "<<  (*it)[l1];
		    }
		fprintf(fi, "%15.8f ", (*it)[dim]);
		//cout<< " "<<(*it)[dim];
		for (int dc = 0; dc < (int) dimC; dc++)
		    {
		    fprintf(fi, "%15.8f ", (*itc)[dc]);
		    //		cout<< " "<<(*itc)[dc]<<endl;
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
    return time;

    }

int ViabiMicroMacroTrajectoryHelper::findOptiControl_Lmin(double budget,
	double *currentPos, double &dt, int nbStepIter, double stepCoeff,
	double *resPos, double &newBudget, bool &succes)
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
    double xCoordsDouble[dim], doubleVect1[dim];
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

    for (int i = 0; i < (int) dim; i++)
	{
	xCoordsDouble[i] = currentPos[i];
	}

    int bestCu;
    rho = dt;
    dynsys->setRho(rho);
    // cout<< " find control :  rho= "<<rho<<endl;

    /*
     * on parcours tous les contrôles
     */

    bestCu = 0;
    double minVal = PLUS_INF, minValCell;
    int iter = 0;
    testNonVide = false;
    cellNum = grid->localizePoint(currentPos);
    /*
     * On parcourt les sommets de la maille
     * autour du sucesseur pour voir s'il y a des
     * points viables
     */
    double currentVal = PLUS_INF;
    for (int ii = 0; ii < nbPointsCube; ii++)
	{
	posTemp = cellNum + indicesDecalCell[ii];
	grid->numToIntAndDoubleCoords(posTemp, testI, testV);
	currentVal = min(vTab[posTemp], currentVal);
	}

    testNonVide = false;
    succes = false;
    while (iter < nbStepIter && !testNonVide)
	{
	unsigned long long int cu = 0;
	dt = rho;
	minVal = PLUS_INF;
	while (cu < nbCTotal)
	    {
	    /*
	     * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	     * au point en cours
	     */
	    if (dynsys->constraintsXU(xCoordsDouble,
		    controlCoords[cu])<PLUS_INF)
		{
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

			double tempL = dynsys->lFunc(xCoordsDouble,
				controlCoords[cu]);

			double tempL1 = dynsys->lFunc(doubleVect1,
				controlCoords[cu]);

			double newVal = budget - rho * 0.5 * (tempL + tempL1);
			//       cout<< " budget  = "<<budget<< "L= "<<tempL<< "L1 " <<tempL1<<  " rho= "<<rho<<endl;

			cellNum = grid->localizePoint(doubleVect1);
			// cout<< " num cellule "<<cellNum<<endl;
			//   iCell=0;
			/*
			 * On parcours les sommets de la maille
			 * autour du sucesseur pour voir s'il y a des
			 * points viables
			 */
			minValCell = PLUS_INF;
			for (int ii = 0; ii < nbPointsCube; ii++)
			    {
			    posTemp = cellNum + indicesDecalCell[ii];
			    grid->numToIntAndDoubleCoords(posTemp, testI,
				    testV);
			    //    cout<< " budget  = "<<budget<< " new val= "<<newVal<< " val point cellule image " <<vTab[posTemp]<< endl;
			    minValCell = min(vTab[posTemp], minValCell);
			    /* if(vTab[posTemp]<=newVal)
			     {
			     if(  newVal-vTab[posTemp]<minValCell )
			     {
			     minValCell=min( newVal-vTab[posTemp],minValCell );
			     optNewVal=newVal;
			     }
			     }
			     */

			    }

			if (minValCell < minVal)
			    {
			    minVal = minValCell;
			    newBudget = newVal;
			    bestCu = cu;
			    cout << " min val = " << minVal << " current val= "
				    << currentVal << " test= " << testNonVide
				    << endl;

			    }
			testNonVide = (minVal <= currentVal);
			//   testNonVide=(minVal<PLUS_INF);
			}
		    }
		}
	    cu++;
	    }				//fin de parcours de tous les contrôles
					// la boucle s'arête ici u premier contrôle
					// qui donne un successeur viable

	// cout<<   " fini parcours de controles on a test interieur = "<<testviabInt<<
	//     " test non vide "<<testNonVide<< " maxnbViabPoints =  "<<maxnbViabPoints<< " bes c u= "<<bestCu<<endl;
	iter++;
	rho = rho * stepCoeff;
	dynsys->setRho(rho);
	cout << "find control  iteration  " << iter << " rho= " << rho << endl;
	succes = testNonVide;
	cout << " succes = " << succes << endl;
	}

    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[bestCu],
	    doubleVect1, rho);
    for (int i = 0; i < (int) dim; i++)
	{
	resPos[i] = doubleVect1[i];
	}
    dt = rho;
    return bestCu;
    }

int ViabiMicroMacroTrajectoryHelper::findOptiControl(double *currentPos,
	unsigned long long int optimDiscreteSuccessor, double &dt,
	int nbStepIter, double stepCoeff, double *resPos, bool &succes)
    {

    // first step : find optimal sucessor for discrete point : guaranteed by Capture Basin
    double realCoordsOnDiscreteTraj[dim];
    unsigned long long int intCoordsOnDiscreteTraj[dim];
    double rho = dt;
    dynsys->setRho(rho);
    /*
     * le pas de temps déterminé localement pour chaque point de la trajectoire
     */
    int bestCu;

    grid->numToIntAndDoubleCoords(optimDiscreteSuccessor,
	    intCoordsOnDiscreteTraj, realCoordsOnDiscreteTraj);

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
     * coordonnées du points corent de la trajectoire
     */
    double xCoordsDouble[dim];
    int posTemp;

    /*
     * tests de validité de point initial
     */
    bool testNonVide = false;
    // int cptOK=0;

    cout << "  debut find control optim tmin : current pos = ";
    for (int i = 0; i < (int) dim; i++)
	{
	xCoordsDouble[i] = currentPos[i];
	cout << " " << xCoordsDouble[i] << endl;
	}

    cout << " find control :  rho= " << rho << endl;

    /*
     * on parcours tous les contrôles
     */

    bestCu = 0;
    int iter = 0;
    double rhoOptim = rho;
    testNonVide = false;
    double minDist = PLUS_INF;
    while (iter < nbStepIter)
	{
	//cout<< "find control  iteration  "<<iter<<" rho= "<<rho<<endl;
	//cout<< " au debut on a testNonVide "<< testNonVide<<endl;
	unsigned long long int cu = 0;

	while (cu < nbCTotal)
	    {
	    /*
	     * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	     * au point en cours
	     */
	    if (dynsys->constraintsXU(xCoordsDouble,
		    controlCoords[cu])<PLUS_INF)
		{
		(dynsys->*(dynsys->discretDynamics))(xCoordsDouble,
			controlCoords[cu], imageVect, rho);
		if (grid->isPointInGrid(imageVect))
		    {
		    /*
		     *  le sucesseur est dans la grille de calcul
		     */
		    if (dynsys->constraintsX(imageVect) < PLUS_INF)
			{
			if (grid->ArePointsInTheSameCell(imageVect,
				realCoordsOnDiscreteTraj))
			    {
			    testNonVide = true;
			    double dist = 0.0;
			    for (int k = 0; k < dim; k++)
				{
				dist =
					max(dist,
						abs(
							imageVect[k]
								- realCoordsOnDiscreteTraj[k]));
				}
			    if (dist < minDist)
				{
				minDist = dist;
				bestCu = cu;
				rhoOptim = rho;
				}
			    }
			}
		    }
		}
	    cu++;
	    }		//fin de parcours de tous les contrôles
			// la boucle s'arête ici u premier contrôle
			// qui donne un successeur viable

	iter++;
	rho = rho * stepCoeff;
	dynsys->setRho(rho);

	}

    if (testNonVide)
	cout << " SUCCESS : TROUVE OK \n";
    succes = testNonVide;
    dynsys->setRho(rhoOptim);
    cout << "  fin de recherche de controle optima on a bestCu nu m " << bestCu
	    << " optimal rho = " << rhoOptim << endl;
    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[bestCu],
	    imageVect, rhoOptim);
    for (int i = 0; i < (int) dim; i++)
	{
	resPos[i] = imageVect[i];
	}
    dt = rhoOptim;
    return bestCu;
    }

unsigned long long int ViabiMicroMacroTrajectoryHelper::findOptimalDiscreteSuccessor(
	unsigned long long int pos, double dt)
    {
    int nbPointsCube = (int) pow(2.0, dim);		//pow(2.0, dim);
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */
    unsigned long long int optimSuccessor = grid->nbTotalPoints + 1;
    // first step : find optimal sucessor for discrete point : guaranteed by Capture Basin
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

    cout << "  debut find control optim tmin : current pos = ";

    rho = dt;
    dynsys->setRho(rho);
    cout << " find control :  rho= " << rho << endl;
    double currentVal = PLUS_INF;

    currentVal = vTab[pos];
    cout << " val of current point is !!!!!!!! " << currentVal << endl;
    /*
     * on parcours tous les contrôles
     */

    double minVal = PLUS_INF;
    int iter = 0;

    //cout<< "find control  iteration  "<<iter<<" rho= "<<rho<<endl;
    //cout<< " au debut on a testNonVide "<< testNonVide<<endl;
    unsigned long long int cu = 0;

    minVal = PLUS_INF;
    while (cu < nbCTotal)
	{
	/*
	 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	 * au point en cours
	 */
	if (dynsys->constraintsXU(xCoordsDouble, controlCoords[cu]) < PLUS_INF)
	    {
	    (dynsys->*(dynsys->discretDynamics))(xCoordsDouble,
		    controlCoords[cu], imageVect, rho);
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

		    for (int ii = 0; ii < nbPointsCube; ii++)
			{
			posTemp = cellNum + indicesDecalCell[ii];
			if (vTab[posTemp] < minVal)
			    {
			    optimSuccessor = posTemp;
			    minVal = vTab[posTemp];
			    }
			}
		    }
		}
	    }
	cu++;
	}	//fin de parcours de tous les contrôles
		// la boucle s'arête ici u premier contrôle
		// qui donne un successeur viable
    cout << " min val = " << minVal << " currentVal = " << currentVal << endl;

    return optimSuccessor;
    }

unsigned long long int ViabiMicroMacroTrajectoryHelper::findViabControlDefault_DD(
	double budget, unsigned long long int *currentPos,
	unsigned long long int currentControl, unsigned long long int *resPos,
	double &newBudget, bool &succes)
    {
    int nbPointsCube = (int) pow(2.0, dim);	//pow(2.0, dim);
    unsigned long long int nbPointsTotal = grid->getNbTotalPoints();

    unsigned long long int **controlCoords = dynsys->getControlIntCoords();
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */

    unsigned long long int testI[dim];
    double testV[dim];

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
    unsigned long long int xCoords[dim], intVect1[dim];

    /*
     * numéros de mailles
     */
    int cellNum;

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
    double minVal = PLUS_INF, minValCell;
    int iter = 0;
    testNonVide = false;
    /*
     * On parcourt les sommets de la maille
     * autour du sucesseur pour voir s'il y a des
     * points viables
     */
    double currentVal = PLUS_INF;

    grid->intCoordsToNum(xCoords, &posTemp);
    currentVal = vTab[posTemp];

    testNonVide = false;
    succes = false;

    unsigned long long int cu = 0;
    bool testPrint = false; // (xCoords[0] == 27) && (xCoords[1] == 27 ) && (xCoords[2] == 6)&& (xCoords[3] == 7);
    minVal = PLUS_INF;
    double tempM;
    bool eligibleControl;
    while (cu < nbCTotal && !testNonVide)
	{
	/*
	 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	 * au point en cours
	 */
	if (testPrint)
	    {
	    cout << " control ";
	    printVector(controlCoords[cu], dimC);
	    cout << " control precedent";
	    printVector(controlCoords[currentControl], dimC);
	    }
	if (currentControl < nbCTotal)
	    {
	    eligibleControl =
		    (dynsys->controlEligibilityForTraj_fd(xCoords,
			    controlCoords[cu], controlCoords[currentControl])
			    < PLUS_INF);
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
		if (testPrint)
		    {
		    cout << " retrour dnamique discrete ";
		    printVector(intVect1, dim);
		    }
		//(dynsys->*discretDynamics)(doublePointCoords, controlCoords[cu], doubleVect);
		if (grid->isPointInGrid_fd(intVect1))
		    {
		    // printf( "   le point est das la grlle\n " );

		    //////printf(" le point est das la grlle\n");
		    /*!
		     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
		     */
		    if (dynsys->constraintsX_fd(intVect1) < PLUS_INF)
			{
			/* cout<< " image du point ";
			 for(int k=0;k<dim;k++)
			 {
			 cout<< " "<<intVect1[k];
			 }*/
			/*
			 * le successeur vérifie les contraintes
			 * On identifie la maille où il se trouve
			 */

			tempM = dynsys->muFunc_fd(xCoords, controlCoords[cu]);
			//cout<< " controle "; printVector(controlCoords[cu], dimC);
			//cout<< " deficit max " << tempM<<endl;
			if (budget >= tempM)
			    {
			    double tempL = dynsys->lFunc_fd(xCoords,
				    controlCoords[cu]);

			    double newVal = budget + tempL;

			    grid->intCoordsToNum(intVect1, &imagePos);
			    testNonVide = (newVal >= vTab[imagePos]);
			    if (testPrint)
				{
				cout << " budget  = " << budget << "L= "
					<< tempL << "newVal " << newVal
					<< " vTab[imagePos]= " << vTab[imagePos]
					<< endl;
				}
			    // cout<< " num cellule "<<cellNum<<endl;
			    //   iCell=0;
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
				//	cout<< " new budget = "<<newBudget<< " image val = "<<vTab[imagePos]<< " test= "<<testNonVide<<endl;

				}
			    }
			else
			    {
			    cout
				    << " controle ne vérifie pas la contrainte de budget\n";
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

unsigned long long int ViabiMicroMacroTrajectoryHelper::findViabControlDefault_tych_DD(
	double budget, unsigned long long int *currentPos,
	unsigned long long int currentControl,
	unsigned long long int currentTych, unsigned long long int *resPos,
	double &newBudget, bool &succes)
    {
    int nbPointsCube = (int) pow(2.0, dim);	//pow(2.0, dim);
    unsigned long long int nbPointsTotal = grid->getNbTotalPoints();

    unsigned long long int **controlCoords = dynsys->getControlIntCoords();
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */

    unsigned long long int **tychIntCoords = dynsys->getTychIntCoords();
    unsigned long long int nbTy = dynsys->getTotalNbPointsTy();

    int dimTy = dynsys->getDimTy();
    double testV[dim];
    unsigned long long int testI[dim];

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
    unsigned long long int xCoords[dim], intVect1[dim];

    /*
     * numéros de mailles
     */
    int cellNum;

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
    double minVal = PLUS_INF, minValCell;
    int iter = 0;
    testNonVide = false;
    /*
     * On parcourt les sommets de la maille
     * autour du sucesseur pour voir s'il y a des
     * points viables
     */
    double currentVal = PLUS_INF;

    grid->intCoordsToNum(xCoords, &posTemp);
    currentVal = vTab[posTemp];

    testNonVide = false;
    succes = false;

    unsigned long long int cu = 0;
    bool testPrint = false; // (xCoords[0] == 27) && (xCoords[1] == 27 ) && (xCoords[2] == 6)&& (xCoords[3] == 7);
    minVal = PLUS_INF;
    double tempM;
    bool eligibleControl;
    while (cu < nbCTotal && !testNonVide)
	{
	/*
	 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	 * au point en cours
	 */
	if (testPrint)
	    {
	    cout << " control ";
	    printVector(controlCoords[cu], dimC);
	    cout << " control precedent";
	    printVector(controlCoords[currentControl], dimC);
	    }
	if (currentControl < nbCTotal)
	    {
	    eligibleControl =
		    (dynsys->controlEligibilityForTraj_fd(xCoords,
			    controlCoords[cu], controlCoords[currentControl])
			    < PLUS_INF);
	    }
	else
	    {
	    eligibleControl = true;
	    }

	if (eligibleControl)
	    {
	    if (dynsys->constraintsXU_fd(xCoords, controlCoords[cu]) < PLUS_INF)
		{
		dynsys->dynamics_tych_fd(xCoords, controlCoords[cu],
			tychIntCoords[currentTych], intVect1);
		if (testPrint)
		    {
		    cout << " retrour dnamique discrete ";
		    printVector(intVect1, dim);
		    }
		//(dynsys->*discretDynamics)(doublePointCoords, controlCoords[cu], doubleVect);
		if (grid->isPointInGrid_fd(intVect1))
		    {
		    // printf( "   le point est das la grlle\n " );

		    //////printf(" le point est das la grlle\n");
		    /*!
		     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
		     */
		    if (dynsys->constraintsX_fd(intVect1) < PLUS_INF)
			{
			/* cout<< " image du point ";
			 for(int k=0;k<dim;k++)
			 {
			 cout<< " "<<intVect1[k];
			 }*/
			/*
			 * le successeur vérifie les contraintes
			 * On identifie la maille où il se trouve
			 */

			tempM = dynsys->muFunc_fd(xCoords, controlCoords[cu]);
			//cout<< " controle "; printVector(controlCoords[cu], dimC);
			//cout<< " deficit max " << tempM<<endl;
			if (budget >= tempM)
			    {
			    double tempL = dynsys->lFunc_tych_fd(xCoords,
				    controlCoords[cu],
				    tychIntCoords[currentTych]);

			    double newVal = budget + tempL;

			    grid->intCoordsToNum(intVect1, &imagePos);
			    testNonVide = (newVal >= vTab[imagePos]);
			    if (testPrint)
				{
				cout << " budget  = " << budget << "L= "
					<< tempL << "newVal " << newVal
					<< " vTab[imagePos]= " << vTab[imagePos]
					<< endl;
				}
			    // cout<< " num cellule "<<cellNum<<endl;
			    //   iCell=0;
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
				//	cout<< " new budget = "<<newBudget<< " image val = "<<vTab[imagePos]<< " test= "<<testNonVide<<endl;

				}
			    }
			else
			    {
			    cout
				    << " controle ne vérifie pas la contrainte de budget\n";
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

unsigned long long int ViabiMicroMacroTrajectoryHelper::findViabControlMinValue_DD(
	double budget, unsigned long long int *currentPos,
	unsigned long long int currentControl, unsigned long long int *resPos,
	double &newBudget, bool &succes)
    {
    int nbPointsCube = (int) pow(2.0, dim);	//pow(2.0, dim);
    unsigned long long int nbPointsTotal = grid->getNbTotalPoints();

    unsigned long long int **controlCoords = dynsys->getControlIntCoords();
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */

    unsigned long long int testI[dim];
    double testV[dim];

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
    unsigned long long int xCoords[dim], intVect1[dim];

    /*
     * numéros de mailles
     */
    int cellNum;

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
    double minVal = PLUS_INF, minValCell;
    int iter = 0;
    testNonVide = false;
    /*
     * On parcourt les sommets de la maille
     * autour du sucesseur pour voir s'il y a des
     * points viables
     */
    double currentVal = PLUS_INF;

    grid->intCoordsToNum(xCoords, &posTemp);
    currentVal = vTab[posTemp];

    testNonVide = false;
    succes = false;

    unsigned long long int cu = 0;

    minVal = PLUS_INF;
    while (cu < nbCTotal)
	{
	/*
	 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	 * au point en cours
	 */
	if (dynsys->controlEligibilityForTraj_fd(xCoords, controlCoords[cu],
		controlCoords[currentControl]) < PLUS_INF)
	    {
	    if (dynsys->constraintsXU_fd(xCoords, controlCoords[cu]) < PLUS_INF)
		{
		dynsys->dynamics_fd(xCoords, controlCoords[cu], intVect1);

		// cout<< " retrour dnamique discrete ";
		// printVector(intVect1, dim);
		//(dynsys->*discretDynamics)(doublePointCoords, controlCoords[cu], doubleVect);
		if (grid->isPointInGrid_fd(intVect1))
		    {
		    // printf( "   le point est das la grlle\n " );

		    //////printf(" le point est das la grlle\n");
		    /*!
		     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
		     */
		    if (dynsys->constraintsX_fd(intVect1) < PLUS_INF)
			{
			/* cout<< " image du point ";
			 for(int k=0;k<dim;k++)
			 {
			 cout<< " "<<intVect1[k];
			 }*/
			/*
			 * le successeur vérifie les contraintes
			 * On identifie la maille où il se trouve
			 */
			double tempM = dynsys->muFunc_fd(xCoords,
				controlCoords[cu]);
			if (budget >= tempM)
			    {
			    double tempL = dynsys->lFunc_fd(xCoords,
				    controlCoords[cu]);

			    double newVal = budget + tempL;
			    //       cout<< " budget  = "<<budget<< "L= "<<tempL<< "L1 " <<tempL1<<  " rho= "<<rho<<endl;

			    grid->intCoordsToNum(intVect1, &imagePos);
			    testNonVide = (newVal >= vTab[imagePos]);
			    // cout<< " num cellule "<<cellNum<<endl;
			    //   iCell=0;
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
				    //cout<< " new budget = "<<newBudget<< " image val = "<<vTab[imagePos]<< " test= "<<testNonVide<<endl;
				    }
				}
			    }
			else
			    {
			    cout
				    << " controle ne vérifie pas la contrainte de budget\n";
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

unsigned long long int ViabiMicroMacroTrajectoryHelper::findViabControlDiffControl_DD(
	double budget, unsigned long long int *currentPos,
	unsigned long long int currentControl, unsigned long long int *resPos,
	double &newBudget, bool &succes)
    {
    int nbPointsCube = (int) pow(2.0, dim);	//pow(2.0, dim);
    unsigned long long int nbPointsTotal = grid->getNbTotalPoints();

    unsigned long long int **controlCoords = dynsys->getControlIntCoords();
    /*
     * tableaux temporaires pour récupérer les indices du point dans la
     * grille
     */

    unsigned long long int testI[dim];
    double testV[dim];

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
    unsigned long long int xCoords[dim], intVect1[dim];

    /*
     * numéros de mailles
     */
    int cellNum;

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
    double minVal = PLUS_INF, minValCell;
    int iter = 0;
    testNonVide = false;
    /*
     * On parcourt les sommets de la maille
     * autour du sucesseur pour voir s'il y a des
     * points viables
     */
    double currentVal = PLUS_INF;

    grid->intCoordsToNum(xCoords, &posTemp);
    currentVal = vTab[posTemp];

    testNonVide = false;
    succes = false;

    unsigned long long int cu = 0;

    minVal = PLUS_INF;
    while (cu < nbCTotal)
	{
	/*
	 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
	 * au point en cours
	 */
	if (dynsys->controlEligibilityForTraj_fd(xCoords, controlCoords[cu],
		controlCoords[currentControl]) < PLUS_INF)
	    {
	    if (dynsys->constraintsXU_fd(xCoords, controlCoords[cu]) < PLUS_INF)
		{
		dynsys->dynamics_fd(xCoords, controlCoords[cu], intVect1);

		// cout<< " retrour dnamique discrete ";
		// printVector(intVect1, dim);
		//(dynsys->*discretDynamics)(doublePointCoords, controlCoords[cu], doubleVect);
		if (grid->isPointInGrid_fd(intVect1))
		    {
		    // printf( "   le point est das la grlle\n " );

		    //////printf(" le point est das la grlle\n");
		    /*!
		     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
		     */
		    if (dynsys->constraintsX_fd(intVect1) < PLUS_INF)
			{
			/* cout<< " image du point ";
			 for(int k=0;k<dim;k++)
			 {
			 cout<< " "<<intVect1[k];
			 }*/
			/*
			 * le successeur vérifie les contraintes
			 * On identifie la maille où il se trouve
			 */
			double tempM = dynsys->muFunc_fd(xCoords,
				controlCoords[cu]);
			if (budget >= tempM)
			    {
			    bool testEqualControls = (cu == currentControl);

			    double tempL = dynsys->lFunc_fd(xCoords,
				    controlCoords[cu]);

			    double newVal = budget + tempL;
			    //       cout<< " budget  = "<<budget<< "L= "<<tempL<< "L1 " <<tempL1<<  " rho= "<<rho<<endl;

			    grid->intCoordsToNum(intVect1, &imagePos);
			    if (!testNonVide && newVal >= vTab[imagePos])
				{
				cout << " Controle viable trouve " << cu
					<< endl;
				testNonVide = true;
				succes = testNonVide;
				newBudget = budget + tempL;
				bestCu = cu;
				//cout<< " new budget = "<<newBudget<< " image val = "<<vTab[imagePos]<< " testNonVide= "<<testNonVide<< " currentCu " << currentControl<<endl;
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
				    cout
					    << " controle viable different trouve :new budget = "
					    << newBudget << " image val = "
					    << vTab[imagePos]
					    << " testNonVide= " << testNonVide
					    << " currentCu " << currentControl
					    << endl;

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
			    cout
				    << " controle ne vérifie pas la contrainte de budget\n";
			    }
			}
		    }
		}
	    }
	cu++;
	}				//fin de parcours de tous les contrôles
					// la boucle s'arête ici u premier contrôle
					// qui donne un successeur viable

    cout << " on retourne bestCu = " << bestCu << " succes = " << succes
	    << endl;
    return bestCu;
    }

double ViabiMicroMacroTrajectoryHelper::computeViableTrajectory_DD(
	unsigned long long int *initPosition, double initValue, string fileName,
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
    double T = dynsys->getTimeHorizon();
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
	    grid->intCoordsToNum(initPosition, &posTemp);

	    cout << " pos temp num = " << posTemp << endl;
	    testNonVide = (vTab[posTemp] < PLUS_INF);

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
		unsigned long long int Cu;
		c = dynsys->target_fd(xCoordsInt);
		cout << " cible donne au point init " << c << endl;
		while ((c >= PLUS_INF) && (testNonVide)
			&& (nbIter <= NB_MAX_TRAJ_ITER))
		    {
		    cout << " point en cours ";
		    for (int i = 0; i < (int) dim; i++)
			{
			newTrajPoint[i] = xCoordsInt[i];
			cout << " " << newTrajPoint[i];
			}
		    grid->intCoordsToNum(xCoordsInt, &posTemp);
		    cout << " pos temp num = " << posTemp << endl;
		    cout << " value of   point  " << vTab[posTemp] << endl;
		    cout << " budget initial  " << budget << endl;
		    traj.push_back(newTrajPoint);
		    valsOptiTraj.push_back(vTab[posTemp]);
		    valsRealTraj.push_back(budget);

		    bestCu = (this->*findfViableControl_DD)(budget, xCoordsInt,
			    currentCu, imageVect, newBudget, testNonVide);
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
			cout << "controle viable trouve : best Cu = " << bestCu
				<< " currentCu = " << currentCu << " ";
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
		    grid->intCoordsToNum(xCoordsInt, &posTemp);
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
	    printf("** error: impossible to open the file %s.\n",
		    fileName.c_str());
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
		    fprintf(fi, "%15d ", (*it)[l1]);
		    //cout<< " "<<  (*it)[l1];
		    }
		for (int dc = 0; dc < dimC; dc++)
		    {
		    fprintf(fi, "%15d ", (*itc)[dc]);
		    //cout<< " "<<(*itc)[dc]<<endl;
		    }
		fprintf(fi, "%15.8f ", (*itVR));
		fprintf(fi, "%15.8f ", (*itVO));
		fprintf(fi, "\n");
		it++;
		itc++;
		itVR++;
		itVO++;
		//   traj.pop_front();
		//   trajC.pop_front();
		}
	    fclose(fi);
	    }
	}
    return newBudget;
    }

double ViabiMicroMacroTrajectoryHelper::computeViableTrajectory_tych_DD(
	unsigned long long int *initPosition, double initValue, string fileName,
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
    double T = dynsys->getTimeHorizon();
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
	    grid->intCoordsToNum(initPosition, &posTemp);

	    cout << " pos temp num = " << posTemp << endl;
	    testNonVide = (vTab[posTemp] < PLUS_INF);

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
		unsigned long long int Cu;
		c = dynsys->target_fd(xCoordsInt);
		cout << " cible donne au point init " << c << endl;
		while ((c >= PLUS_INF) && (testNonVide)
			&& (nbIter <= NB_MAX_TRAJ_ITER))
		    {
		    cout << " point en cours ";
		    for (int i = 0; i < (int) dim; i++)
			{
			newTrajPoint[i] = xCoordsInt[i];
			cout << " " << newTrajPoint[i];
			}
		    grid->intCoordsToNum(xCoordsInt, &posTemp);
		    cout << " pos temp num = " << posTemp << endl;
		    cout << " value of   point  " << vTab[posTemp] << endl;
		    cout << " budget initial  " << budget << endl;
		    traj.push_back(newTrajPoint);
		    valsOptiTraj.push_back(vTab[posTemp]);
		    valsRealTraj.push_back(budget);

		    unsigned long long int currentTych =
			    (unsigned long long int) (rand() % nbTy);
		    bestCu = findViabControlDefault_tych_DD(budget, xCoordsInt,
			    currentCu, currentTych, imageVect, newBudget,
			    testNonVide);
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
			cout << "controle viable trouve : best Cu = " << bestCu
				<< " currentCu = " << currentCu << " ";
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
			    trajTyControlCoords[dc] =
				    tychIntCoords[currentTych][dc];
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
		    grid->intCoordsToNum(xCoordsInt, &posTemp);
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
	    printf("** error: impossible to open the file %s.\n",
		    fileName.c_str());
	    }
	else
	    {
	    list<valarray<unsigned long long int> >::iterator it = traj.begin();
	    list<valarray<unsigned long long int> >::iterator itc = trajC.end();
	    list<valarray<unsigned long long int> >::iterator itTy =
		    trajTy.end();
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
		    fprintf(fi, "%15d ", (*it)[l1]);
		    //cout<< " "<<  (*it)[l1];
		    }
		for (int dc = 0; dc < dimC; dc++)
		    {
		    fprintf(fi, "%15d ", (*itc)[dc]);
		    //cout<< " "<<(*itc)[dc]<<endl;
		    }
		for (int dc = 0; dc < dimTy; dc++)
		    {
		    fprintf(fi, "%15d ", (*itTy)[dc]);
		    //cout<< " "<<(*itc)[dc]<<endl;
		    }
		fprintf(fi, "%15.8f ", (*itVR));
		fprintf(fi, "%15.8f ", (*itVO));
		fprintf(fi, "\n");
		it++;
		itc++;
		itTy++;
		itVR++;
		itVO++;
		//   traj.pop_front();
		//   trajC.pop_front();
		}
	    fclose(fi);
	    }
	}
    return newBudget;
    }

double ViabiMicroMacroTrajectoryHelper::computeViableTrajectory(
	double *initPosition, double initValue, string fileName, bool &succes)
    {
    succes = false;
    return 0.0;
    }

ViabiMicroMacroTrajectoryHelper::~ViabiMicroMacroTrajectoryHelper()
    {
    // TODO Auto-generated destructor stub
    }

