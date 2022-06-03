/*
f * ViabiHJB.cpp
 *  *
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
 *  Created on: 9 déc. 2013
 *      Author: ANYA
 */

#include "../include/ViabiMicroMacro.h"

ViabiMicroMacro::ViabiMicroMacro(ParametersManager * pm)
{
	modelParams = pm;
	/*
	 * instanciation du système dynamique
	 */
	systemParams sp= *(pm->getSystemParameters());

	algoViabiParams avp =*(pm->getAlgoParameters());

	controlParams cp =*(pm->getControlParameters());
	gridParams gp = *(pm->getGridParameters());

	dim = gp.DIM;
	grid=new GridMicroMacro( gp);
	grid->printGrid();

	dynsys= new SysDyn(sp, dim, cp, grid);
	InitViabiMicroMacro(avp);
}

void ViabiMicroMacro::InitViabiMicroMacro(algoViabiParams avp) {
	// TODO Auto-generated constructor stub
	/*!
	 * l'initialisation commence par attribution de valeurs
	 * au pointeurs de grille et de systeme
	 */
	ostringstream os;

	vTab=grid->getGridPtr();

	dim=grid->dim;
	dimC=dynsys->getDimC();


	nbOMPThreads=avp.NB_OMP_THREADS;



	cout<< " classe viabi Micro Macro   dim  d'état = "<<dim<<endl;

	/*!
	 * On initialise la structure  qui servira au stockage temporaire de l'image discrete d'un point de grille
	 * Les tableaux sont initialisés  avec leur taille maximale possible :
	 * égale au nb  de controles.
	 * C'est la valeur effective  de nbImageCells , calculé é chaque evaluation de l'image discrete
	 * d'un point qui  servira  é lire et remplir  correctement
	 * la bonne partie de ces tableaux
	 */
	pointDI.tabCellEntrees=new unsigned long long int[dynsys->getTotalNbPointsC()+1];
	pointDI.tabImageCells=new unsigned long long int[dynsys->getTotalNbPointsC()];
	pointDI.tabImageControls=new unsigned long long int[dynsys->getTotalNbPointsC()];

	/*!
	 * On initialise la structure  qui servira au stockage temporaire de l'image discrete d'un point de grille
	 * Les tableaux sont initialisés  avec leur taille maximale possible :
	 * égale au nb  de controles.
	 * C'est la valeur effective  de nbImageCells , calculé é chaque evaluation de l'image discrete
	 * d'un point qui  servira  é lire et remplir  correctement
	 * la bonne partie de ces tableaux
	 */
	pointDI_DD.tabPointsEntrees=new unsigned long long int[dynsys->getTotalNbPointsC()+1];
	pointDI_DD.tabImagePoints=new unsigned long long int[dynsys->getTotalNbPointsC()];
	pointDI_DD.tabImageControls=new unsigned long long int[dynsys->getTotalNbPointsC()];

	/*
	 * Initialisation de tableaux servant de variables globales pour certaines fonctions
	 * Cela évte de multiples allocations/destructions de mémoire
	 */
	intPointCoords=new unsigned long long int[dim];
	intVect1=new unsigned long long int[dim];
	doublePointCoords=new double[dim];
	doubleVect=new double[dim];
	doubleVect1=new double[dim];
	intControlCoords=new unsigned long long int[dimC];
	doubleControlCoords=new double[dimC];
	imageCells=new unsigned long long int[dynsys->getTotalNbPointsC()];

	filePrefix=avp.FILE_PREFIX;


	targ_or_dep=avp.TARGET_OR_DEPARTURE;
	computeTmin=avp.COMPUTE_TMIN;


	ViabiMicroMacro::addNewPointsToSet=&ViabiMicroMacro::addNewPoints;
	ViabiMicroMacro::addDataToCurrentCell=&ViabiMicroMacro::addDataToCell;
	ViabiMicroMacro::addDataToCurrentPoint=&ViabiMicroMacro::addDataToPoint;

	if((dynsys->getDynType()==DD)   )
	{
		cout<< " on utilise les m�thodes DD"<<endl;
		ViabiMicroMacro::createCurrentPointsList=&ViabiMicroMacro::createPointsList_DD;
	}
	else
	{
		ViabiMicroMacro::createCurrentPointsList=&ViabiMicroMacro::createPointsList;
	}


	if(avp.COMPUTE_TMIN)
	{
		ViabiMicroMacro::computeFirstConvexifiedImage=&ViabiMicroMacro::computeConvexifiedImage_tmin;
		ViabiMicroMacro::computeCurrentImage=&ViabiMicroMacro::computeCurrIm_tmin;
		//  ViabiMicroMacro::computeFirstConvexifiedImage_omp=&ViabiMicroMacro::computeConvexifiedImage_tmin_omp;
	}
	else
	{
		//  ViabiMicroMacro::computeFirstConvexifiedImage_omp=&ViabiMicroMacro::computeConvexifiedImage_Lmin_omp;
		if((dynsys->getDynType()==DD)   )
		{
			cout<< " on est en mode DD et Lmin "<<endl;
			ViabiMicroMacro::computeCurrentImage=&ViabiMicroMacro::computeCurrIm_DD;
			ViabiMicroMacro::computeFirstConvexifiedImage=&ViabiMicroMacro::computeConvexifiedImage_DD;
		}
		else
		{
			ViabiMicroMacro::computeCurrentImage=&ViabiMicroMacro::computeCurrIm_Lmin;
			ViabiMicroMacro::computeFirstConvexifiedImage=&ViabiMicroMacro::computeConvexifiedImage_Lmin;
		}

	}

	tempPointsList1 = new list<imagePoint>();
	tempPointsList2 = new list<imagePoint>();;
	whichPointListToUse = 1;
	currentImagePointsList.pointsList = tempPointsList1;
}

void ViabiMicroMacro::computeDiscreteImageOfPoint(unsigned long long int num)
{


	// cout<< " calcul de l'image d'un point \n";
	double ** controlCoords=dynsys->getControlCoords();
	unsigned long long int nbCTotal=dynsys->getTotalNbPointsC();
	unsigned long long int nbCellsTotal=grid->getNbTotalPoints();




	grid->numToIntAndDoubleCoords(num,intPointCoords,doublePointCoords);
	double rho=dynsys->calculRho_local(doublePointCoords);
	// cout<<  " compute DI  pos = "<<num<< " juste apres num to int double coords\n";
	// printVector(doublePointCoords, dim);
	unsigned long long int cu;

	list<intPair> cellsList;
	//double dv[dim];
	for(cu=0;cu<nbCTotal;cu++)
	{

		//

		// cout<< "  controle numero "<<cu<<endl;
		//  	printVector(controlCoords[cu],dimC);
		/*!
		 * on calcule les coordonnées réelles de l'image du point par la dynamique discrete
		 * elles sont stockes dans le tableau doubleVect
		 * si \f$ u\in U(x)\f$  (contréle admissible) on calcule l'image réelle par la dynamique
		 *  discrétisée  en temps  du point x avec le controle u
		 */
		if(dynsys->constraintsXU(doublePointCoords,controlCoords[cu])<PLUS_INF)
		{


			// printf("  contraintes sur U et X  ok\n");
			//  printf( " x= ");
			//  printVector(doublePointCoords,dim);

			(dynsys->*(dynsys->discretDynamics))(doublePointCoords, controlCoords[cu], doubleVect1, rho);

			// cout<< " retrour dnamique discrete ";
			//printVector(doubleVect1, dim);
			//(dynsys->*discretDynamics)(doublePointCoords, controlCoords[cu], doubleVect);
			if(grid->isPointInGrid(doubleVect1))
			{
				//////printf( "   le point est das la grlle\n " );

				//////printf(" le point est das la grlle\n");
				/*!
				 * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
				 */
				if(dynsys->constraintsX(doubleVect1)<PLUS_INF)
				{
					//////printf( "   contraintes Xu ok\n " );

					//	 ////printf("  contraintes sur   X  ok\n");

					/*!
					 * Si l'image  est dans l'ensemble de contraintes sur l'état \f$ K \f$
					 * on calcule le numéro de maille qui contient cett image
					 */
					imageCells[cu]=grid->localizePoint(doubleVect1);
					/*	grid->numToIntAndDoubleCoords(imageCells[cu],testI,testV);
					if(testV[0]<doublePointCoords[0])
					{
						cout<< " num de déparrt "<<num<< " num de cell image "<<imageCells[cu]<<endl;
						cout<< " int coords de départ ";
						for(int lm=0;lm<dim;lm++)
						{
							cout<< " "<<intPointCoords[lm];
						}
						cout<< endl;
						cout<< "double coords de départ ";
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
						cout<< " int coords de result projeté ";
						for(int lm=0;lm<dim;lm++)
						{
							cout<< " "<<testI[lm];
						}
						cout<< endl;
						cout<< "double coords de result projeté  ";
						for(int lm=0;lm<dim;lm++)
						{
							cout<< " "<<testV[lm];
						}
						cout<< endl;

					}*/
					// on enregistre le numero de maille
					cellsList.push_back(intPair(imageCells[cu],cu));
					//////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,imageCells[cu]);
					//   cout<< " les coordonnées  de cin bas de la maille  sont ";
					//   grid->numToIntAndDoubleCoords(imageCells[cu], intPointCoords, doubleVect1);
					//   printVector(doubleVect1, dim);
				}
				else
				{

					/*!
					 * Si l'image n'est pas dans  \f$ K \f$ on enregistre un numéro de maille factice qui signifie que
					 * cette image est rejetée
					 */
					imageCells[cu]=nbCellsTotal+1; // sinon on enregistre un nombre convenu reconnaissanble
					cellsList.push_back(intPair(imageCells[cu],cu));
					//	////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,imageCells[cu]);
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
				if(grid->unboundedDomain && grid->isPointInGridWithConstr(doubleVect1) && (dynsys->constraintsX(doubleVect1)<PLUS_INF))
				{
					// cout<< " sortie autorisee " <<  " on est dans " ;
					// printVector(doubleVect1, dim);
					imageCells[cu]=nbCellsTotal; // sinon on enregistre un nombre convenu reconnaissanble
					cellsList.push_back(intPair(imageCells[cu],cu));
				}
				else
				{
					imageCells[cu]=nbCellsTotal+1; // sinon on enregistre un nombre convenu reconnaissanble
					cellsList.push_back(intPair(imageCells[cu],cu));
				}

				//////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,imageCells[cu]);
			}


			// cout<< " ajour d'une paire é la liste longueur "<<cellsList.size()<<"\n";

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
			imageCells[cu]=nbCellsTotal+1; // sinon on enregistre un nombre convenu reconnaissanble
			cellsList.push_back(intPair(imageCells[cu],cu));
			//printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,imageCells[cu]);
		}

	}
	/* 	cout<< " fini traitement parallele\n";
	 for( cu=0;cu<nbCTotal;cu++)
	 {
	 	cellsList.push_back(intPair(imageCells[cu],cu));
	 cout<< " cu  = "<<cu<< " image cell = "<<imageCells[cu]<<endl;
	 }*/
	////////system("pause");
	/*!
	 * Toutes les images sont calculées  et stockées dans un tableau dans l'ordre
	 *  de numérotation des controles. La tache suivante consiste é trier ce tableau, éliminer les doublons
	 *   et  compléter la structure  image de point
	 */

	cellsList.sort(pairCompare);
	list<intPair>::iterator itCell,itDouble;

	itCell=cellsList.begin();
	itDouble=itCell;
	int currentCell;

	cu=0;
	unsigned long long int iControlTab=0, iCellsTab=0, iEntreesTab=0;
	// cout<< " juste avant while de parcours de la liste  first "<<(*itCell).first<< " nb cells total "<<(int)nbCellsTotal<<endl;
	while((itCell!=cellsList.end()) &( (*itCell).first<=(int)nbCellsTotal))
	{

		currentCell=(*itCell).first;
		if(this->testConstraintesForCell(currentCell))
		{
			//cout<< " current cell "<<currentCell<<endl;
			pointDI.tabCellEntrees[iEntreesTab]=iControlTab;
			pointDI.tabImageCells[iCellsTab]=currentCell;
			iEntreesTab++;
			iCellsTab++;

			while((itDouble!=cellsList.end() )& ((*itDouble).first==currentCell))
			{

				pointDI.tabImageControls[iControlTab]=(*itDouble).second;
				iControlTab++;

				itDouble++;
			}
			/*!
			 * la boucle s'arrete au premier different
			 * é la fin de la boucle soit itDouble a atteint la fin de la liste
			 *  soit  il pointe sur une cellule différente
			 */

			itCell=itDouble;
		}
		else
		{
			while((itDouble!=cellsList.end() )& ((*itDouble).first==currentCell))
			{

				itDouble++;
			}
			/*!
			 * la boucle s'arrete au premier different
			 * é la fin de la boucle soit itDouble a atteint la fin de la liste
			 *  soit  il pointe sur une cellule différente
			 */

			itCell=itDouble;
		}

	}

	pointDI.nbImageCells=iCellsTab;

	/*!
	 * La derniére valeur  du tableau des controles indique la fin  de la liste des controles associés
	 * é la derniére cellule
	 */
	pointDI.tabCellEntrees[iCellsTab]=iControlTab;
	// cout<<  " compute DI  FINFINFINF pos = "<<num<< " double coords\n";
	//  printVector(doublePointCoords, dim);
}

void ViabiMicroMacro::computeDiscreteImageOfPoint_DD(unsigned long long int num)
{


	//  cout<< " calcul de l'image d'un point DD\n";
	unsigned long long int  ** controlCoords=dynsys->getControlIntCoords();
	unsigned long long int nbCTotal=dynsys->getTotalNbPointsC();
	unsigned long long int nbPointsTotal=grid->getNbTotalPoints();

	unsigned long long int imagePos;

	grid->numToIntAndDoubleCoords(num,intPointCoords,doublePointCoords);
	// cout<<  " compute DI  pos = "<<num<< " juste apres num to int double coords\n";
	// printVector(intPointCoords, dim);
	unsigned long long int cu;

	list<uintPair> pointsList;
	//double dv[dim];
	for(cu=0;cu<nbCTotal;cu++)
	{

		//

		//  cout<< "  controle numero "<<cu<<endl;
		// 	printVector(controlCoords[cu],dimC);
		/*!
		 * on calcule les coordonnées réelles de l'image du point par la dynamique discrete
		 * elles sont stockes dans le tableau doubleVect
		 * si \f$ u\in U(x)\f$  (contréle admissible) on calcule l'image réelle par la dynamique
		 *  discrétisée  en temps  du point x avec le controle u
		 */
		if(dynsys->constraintsXU_fd(intPointCoords,controlCoords[cu])<PLUS_INF)
		{


			// printf("  contraintes sur U et X  ok\n");
			//  printf( " x= ");
			//  printVector(doublePointCoords,dim);

			dynsys->dynamics_fd(intPointCoords, controlCoords[cu], intVect1);

			// cout<< " retrour dnamique discrete ";
			// printVector(intVect1, dim);
			//(dynsys->*discretDynamics)(doublePointCoords, controlCoords[cu], doubleVect);
			if(grid->isPointInGrid_fd(intVect1))
			{
				// printf( "   le point est das la grlle\n " );
				grid->intCoordsToNum(intVect1, &imagePos);
				//////printf(" le point est das la grlle\n");
				/*!
				 * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
				 */
				if(vTab[imagePos]<PLUS_INF)
				{
					//printf( "   contraintes Xu ok\n " );

					//	 ////printf("  contraintes sur   X  ok\n");

					/*!
					 * Si l'image  est dans l'ensemble de contraintes sur l'état \f$ K \f$
					 * on calcule le numéro de maille qui contient cett image
					 */
					imageCells[cu] = imagePos;

					/*	grid->numToIntAndDoubleCoords(imageCells[cu],testI,testV);
					if(testV[0]<doublePointCoords[0])
					{
						cout<< " num de déparrt "<<num<< " num de cell image "<<imageCells[cu]<<endl;
						cout<< " int coords de départ ";
						for(int lm=0;lm<dim;lm++)
						{
							cout<< " "<<intPointCoords[lm];
						}
						cout<< endl;
						cout<< "double coords de départ ";
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
						cout<< " int coords de result projeté ";
						for(int lm=0;lm<dim;lm++)
						{
							cout<< " "<<testI[lm];
						}
						cout<< endl;
						cout<< "double coords de result projeté  ";
						for(int lm=0;lm<dim;lm++)
						{
							cout<< " "<<testV[lm];
						}
						cout<< endl;

					}*/
					// on enregistre le numero de maille
					pointsList.push_back(uintPair(imageCells[cu],cu));
					//////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,imageCells[cu]);
					//   cout<< " les coordonnées  de cin bas de la maille  sont ";
					//   grid->numToIntAndDoubleCoords(imageCells[cu], intPointCoords, doubleVect1);
					//   printVector(doubleVect1, dim);
				}
				else
				{
					/*!
					 * Si l'image n'est pas dans  \f$ K \f$ on enregistre un numéro de maille factice qui signifie que
					 * cette image est rejetée
					 */
					//cout<< " pas dans la grille "<<endl;
					imageCells[cu]=nbPointsTotal+1; // sinon on enregistre un nombre convenu reconnaissanble
					try{
						pointsList.push_back(uintPair(imageCells[cu],cu));
					}
					catch (const std::bad_alloc& e) {
						std::cout << "Allocation failed: in discretImage_DD points push_back 1" << e.what() << '\n';
					}
					//	////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,imageCells[cu]);
				}
			}
			else
			{
				//cout<< " pas dans la grille "<<endl;
				/*!
				 * Si l'image n'est pas dans  la grille on enregistre un numéro de maille factice qui signifie que
				 * cette image est rejetée
				 */
				/*!
				 * \todo prévoir la gestion des autorisations de sortie  par axe pour les variables qui peuvent
				 * se trouver dans des intervalles non bornés
				 */
				imageCells[cu]=nbPointsTotal+1; // sinon on enregistre un nombre convenu reconnaissanble
				try{
					pointsList.push_back(uintPair(imageCells[cu],cu));
				}
				catch (const std::bad_alloc& e) {
					std::cout << "Allocation failed: in discretImage_DD points push_back 2" << e.what() << '\n';
				}
				//////printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,imageCells[cu]);
			}


			// cout<< " ajour d'une paire é la liste longueur "<<cellsList.size()<<"\n";

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
			imageCells[cu]=nbPointsTotal+1; // sinon on enregistre un nombre convenu reconnaissanble
			try{
				pointsList.push_back(uintPair(imageCells[cu],cu));
			}
			catch (const std::bad_alloc& e) {
				std::cout << "Allocation failed: in discretImage_DD points push_back 3" << e.what() << '\n';
			}
			//printf( "  thread  numero %d  controle %d   maille visee est %d\n ",ID,cu,imageCells[cu]);
		}

	}
	/* 	cout<< " fini traitement parallele\n";
	 for( cu=0;cu<nbCTotal;cu++)
	 {
	 	cellsList.push_back(intPair(imageCells[cu],cu));
	 cout<< " cu  = "<<cu<< " image cell = "<<imageCells[cu]<<endl;
	 }*/
	////////system("pause");
	/*!
	 * Toutes les images sont calculées  et stockées dans un tableau dans l'ordre
	 *  de numérotation des controles. La tache suivante consiste é trier ce tableau, éliminer les doublons
	 *   et  compléter la structure  image de point
	 */

	pointsList.sort(upairCompare);
	//cout<< " points list size "<<pointsList.size()<<endl;
	list<uintPair>::iterator itCell,itDouble;

	itCell=pointsList.begin();
	itDouble=itCell;
	unsigned long long int currentPoint;

	cu=0;
	unsigned long long int iControlTab=0, iCellsTab=0, iEntreesTab=0;
	// cout<< " juste avant while de parcours de la liste  first "<<(*itCell).first<< " nb cells total "<<(int)nbCellsTotal<<endl;
	while((itCell!=pointsList.end()) &( (*itCell).first<(int)nbPointsTotal))
	{

		currentPoint=(*itCell).first;

		//cout<< " current point "<<currentPoint<<endl;
		pointDI_DD.tabPointsEntrees[iEntreesTab]=iControlTab;
		pointDI_DD.tabImagePoints[iCellsTab]=currentPoint;
		iEntreesTab++;
		iCellsTab++;

		while((itDouble!=pointsList.end() )& ((*itDouble).first==currentPoint))
		{

			pointDI_DD.tabImageControls[iControlTab]=(*itDouble).second;
			iControlTab++;

			itDouble++;
		}
		/*!
		 * la boucle s'arrete au premier different
		 * é la fin de la boucle soit itDouble a atteint la fin de la liste
		 *  soit  il pointe sur une cellule différente
		 */

		itCell=itDouble;


	}

	pointDI_DD.nbImagePoints=iCellsTab;

	/*!
	 * La derniére valeur  du tableau des controles indique la fin  de la liste des controles associés
	 * é la derniére cellule
	 */
	pointDI_DD.tabPointsEntrees[iCellsTab]=iControlTab;
	// cout<<  " compute DI  FINFINFINF pos = "<<num<< " double coords\n";
	//  printVector(doublePointCoords, dim);
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
	if((dynsys->getDynType()==DD)   )
	{
		initialiseTargetHJB_DD();
	}
	else
	{
		initialiseTargetHJB();
	}

	//	unsigned long long int dim=grid->dim;
	//
	//	unsigned long long int  * x=new unsigned long long int [dim];
	//
	//	double * xReel =new double[dim];
	//
	//	int totalPointsX=grid->getNbTotalPoints();
	//	unsigned long long int pos ;
	//
	//
	//	/*
	//	 *  on parcourt  tous les points de l'espace discret  fini
	//	 *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
	//	 */
	//	for( pos=0;pos<(unsigned long long int)totalPointsX;pos++)
	//	{
	//		grid->numToIntAndDoubleCoords(pos,x,xReel);
	//		vTab[pos]=(double)dynsys->target(xReel);
	//	}


}


void ViabiMicroMacro::computeTrajectories()
{
	algoViabiParams * avp=modelParams->getAlgoParameters();
	if(avp->TYPE_TRAJ == OP)
		computeOptimalTrajectories();
	else if(avp->TYPE_TRAJ == VD)
		computeViableTrajectories();

}


void ViabiMicroMacro::computeOptimalTrajectories()
{
	algoViabiParams * avp=modelParams->getAlgoParameters();
	int nbTrajs=avp->NB_TRAJS;
	int typeTraj = avp->TYPE_TRAJ;
	double T = modelParams->getSystemParameters()->maxTime;
	ostringstream os;
	string fileName;
	if(nbTrajs>0)
	{
		bool success[nbTrajs];
		for(int tr=0;tr<nbTrajs;tr++)
		{
			os<<"../OUTPUT/"<<filePrefix<<"-traj-"<<tr+1<<".dat";
			fileName=os.str();
			os.str("");

			computeOptimalCaptTrajectory(avp->INIT_POINTS+tr*dim, fileName, success[tr]);
		}
	}
}


void ViabiMicroMacro::computeViableTrajectories()
{
	algoViabiParams * avp=modelParams->getAlgoParameters();
	int nbTrajs=avp->NB_TRAJS;
	int typeTraj = avp->TYPE_TRAJ;
	double T = modelParams->getSystemParameters()->maxTime;
	ostringstream os;
	string fileName;
	if(nbTrajs>0)
	{
		cout<< "nbTrajectoies is "<< nbTrajs<<endl;
		bool success[nbTrajs];
		for(int tr=0;tr<nbTrajs;tr++)
		{
			os<<"../OUTPUT/"<<filePrefix<<"-traj-"<<tr+1<<".dat";
			fileName=os.str();
			os.str("");


			if(dynsys->getDynType() == DD)
			{
				cout<< "init points : "<<endl;
				for(int k = 0; k<dim; k++)
				{
					cout<< " "<< avp->INIT_POINTS_FD[tr*dim+k];
				}
				cout<<endl;
				computeViableTrajectory_DD(avp->INIT_POINTS_FD+tr*dim, avp->INIT_VALUES_FD[tr],fileName, success[tr]);
			}
			else
				computeViableTrajectory(avp->INIT_POINTS+tr*dim, avp->INIT_VALUES[tr],fileName, success[tr]);
		}
	}
}
double ViabiMicroMacro::computeViableTrajectory(double *initPosition, double initValue,  string fileName, bool &succes)
{
	succes = false;
	return 0.0;
}

double ViabiMicroMacro::computeViableTrajectory_DD(unsigned long long int *initPosition, double initValue, string fileName, bool &succes)
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
	unsigned long long int ** controlCoords=dynsys->getControlIntCoords();


	/*
	 * coordonnées du points corent de la trajectoire
	 */

	unsigned long long int xCoordsInt[dim], imageVect[dim];
	/*
	 * listes  pour contenir la trajectoire ( temps-position) et les contrôles
	 */
	list<valarray<unsigned long long int > > traj;
	list<double > valsOptiTraj;
	list<double > valsRealTraj;
	list<valarray<unsigned long long int > > trajC;

	valarray<unsigned long long int> newTrajPoint(dim);
	double T=dynsys->getTimeHorizon();
	/*
	 * structures accumulables dansune liste
	 */

	valarray<unsigned long long int> trajControlCoords(dimC);

	unsigned long long int posTemp;

	cout<< " calcul de traj a partir de coords \n";
	cout<< " Postion initiale = ";

	for(int l1=0;l1<(int)dim;l1++)
	{
		cout<< " "<<initPosition[l1];
	}
	cout<< " \n";
	/*
	 * tests de validité de point initial
	 */
	bool testNonVide=false;
	double budget, newBudget;
	unsigned long long int bestCu;

	if(grid->isPointInGrid_fd(initPosition))
	{
		if(dynsys->constraintsX_fd(initPosition)<PLUS_INF)
		{
			grid->intCoordsToNum(initPosition, &posTemp);


			testNonVide=(vTab[posTemp]<PLUS_INF);

			if(!testNonVide)
			{
				cout<<" La position initiale sélectionnée n'appartiant pas au noyau de viabilité\n";
				succes=0;
			}
			else
			{
				/*
				 * la position initiale se trouve dans le noyau de viabilité
				 * on initialise le temps à 0  et recopie la pos initiale
				 * dans le coordonnées temporaires du point en cours de la trajectoire
				 */
				if(initValue < vTab[posTemp])
				{
					cout<<"  le budget initial est inf�rieur au budget minimum requis pour cet �tat initial. Construction de trajectoire � partir du budget minimum requis\n";
					initValue = vTab[posTemp];
				}

				budget= initValue;
				cout<< " value of init point  "<<vTab[posTemp]<<endl;
				cout<< " budget initial  "<<budget<<endl;
				for(int i=0;i<(int)dim;i++)
				{
					xCoordsInt[i]=initPosition[i];
				}
				int nbIter=0;
				/*
				 * On itère tant que le temps n'a pas dépassé l'horizon donné
				 */

				double c;
				unsigned long long int Cu;
				c=dynsys->target_fd(xCoordsInt);
				cout<< " cible donne au point init "<< c<<endl;
				while( (c>=PLUS_INF) && (testNonVide) && (nbIter<=NB_MAX_TRAJ_ITER))
				{
					cout<< " point en cours ";
					for(int i=0;i<(int)dim;i++)
					{
						newTrajPoint[i]=xCoordsInt[i];
						cout<< " "<<newTrajPoint[i];
					}
					grid->intCoordsToNum(xCoordsInt, &posTemp);
					traj.push_back(newTrajPoint);
					valsOptiTraj.push_back(vTab[posTemp]);
					valsRealTraj.push_back(budget);


					bestCu=this->findViabControl_DD(budget,xCoordsInt,  imageVect,newBudget, testNonVide );
					// la boucle s'arête ici u premier contrôle
					// qui donne un successeur viable

					//cout<<   "Recherche de controle viable  fini parcours de controles on a  test non vide "<<testNonVide<< " bes c u= "<<bestCu<<endl;

					// contrôle viable trouvé
					// on recopie ce contrôle dans la liste et
					// le successeur devient le point  courent
					if(testNonVide)
					{

						//	cout<<  " image interieure tourvee \n";

						for(int dc=0;dc<(int)dimC;dc++)
						{
							trajControlCoords[dc]=controlCoords[bestCu][dc];
						}
						trajC.push_back(trajControlCoords);
						for(int i=0;i<(int)dim;i++)
						{
							xCoordsInt[i]=imageVect[i];
						}
						budget=newBudget;
						//cout<< " new budget = "<<newBudget<<endl;
					}
					else
					{

						cout<<"   Echec! Sortie de l'ensemble viable \n";
						break;
					}
					c=dynsys->target_fd(xCoordsInt);
					cout<<"   valeur cible : "<<c<<endl;
					nbIter++;
				}
				if(c< PLUS_INF)
				{
					succes=1;
					cout<< " dernier point en cours ";
					for(int i=0;i<(int)dim;i++)
					{
						newTrajPoint[i]=xCoordsInt[i];
						cout<< " "<<newTrajPoint[i];
					}
					grid->intCoordsToNum(xCoordsInt, &posTemp);
					traj.push_back(newTrajPoint);
					valsOptiTraj.push_back(vTab[posTemp]);
					valsRealTraj.push_back(budget);
				}
			}//fin de else (reconstruction de trajectoire)

		}
		else
		{
			printf(" Point initial hors de l'ensemble de contraintes. Arret\n");
			succes=0;
		}
	}
	else
	{
		printf(" Point initial hors de grille de calcul. Arret\n");
		succes=0;
	}

	// if(succes)
	{
		printf(" trajectoire trouvée. Enregistrement\n");

		FILE * fi;
		fi = fopen( fileName.c_str(),"w");
		if(fi==NULL){
			printf("** error: impossible to open the file %s.\n", fileName.c_str());
		}
		else{
			list<valarray<unsigned long long int> >::iterator it=traj.begin();
			list<valarray<unsigned long long int> >::iterator itc=trajC.end();
			list<double >::iterator itVO=valsOptiTraj.begin();
			list<double >::iterator itVR=valsRealTraj.begin();
			itc--;
			trajC.push_back(*itc);
			itc=trajC.begin();

			while(it!=traj.end())
			{

				for(int l1=0;l1<dim;l1++)
				{
					fprintf(fi,  "%15d " ,   (*it)[l1]);
					//cout<< " "<<  (*it)[l1];
				}
				for(int dc=0;dc<dimC;dc++)
				{
					fprintf( fi, "%15d " ,   (*itc)[dc]);
					//cout<< " "<<(*itc)[dc]<<endl;
				}
				fprintf( fi, "%15.8f " ,   (*itVR));
				fprintf( fi, "%15.8f " ,   (*itVO));
				fprintf( fi, "\n" );
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


	double c;


	unsigned long long int dim=grid->dim;

	//ouverture de la base


	unsigned long long int  * x=new unsigned long long int [dim];

	double * xReel =new double[dim];



	int  totalPointsC=0;  // nombre de points  de l'espace des commandes




	int totalPointsX=grid->getNbTotalPoints();

	imagePoint currentPoint;

	currentImagePointsList.maxNum=0;
	currentImagePointsList.minNum=0;
	currentImagePointsList.pointsList= tempPointsList1;;

	unsigned long long int pos ;

	list<imagePoint>::iterator itStart=currentImagePointsList.pointsList->begin(),
			itNew;
	int cpt=0;

	cout<< " total points de grid "<<totalPointsX<<endl;
	/*
	 *  on parcourt  tous les points de l'espace discret  fini
	 *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
	 */
	for( pos=0;pos<(unsigned long long int)totalPointsX;pos++)
	{

		//cout<< " pos= "<<pos<<endl;
		/*!
		 * le compteur pos  ets l'unique numéro entier du point en cours
		 * dans la numérotation alphabétique : on parcourt axe par axe
		 */

		/*!
		 * on restitue les  coordonnées netiéres  du point é partir de son numéro
		 * ainsi que ses coordonnées réelles
		 */
		grid->numToIntAndDoubleCoords(pos,x,xReel);

		c=max( (*(dynsys->target))(xReel),(*(dynsys->constraintsX))(xReel)) ;

		vTab[pos]=c;
		if(c<PLUS_INF)
		{
			totalPointsC++;
			currentPoint.minVal=c;
			currentPoint.PointNum=pos;

			addDataToPointsList(&itStart, currentPoint, &itNew);
			//cout<<" fini ajout point  \n";
			//printVector(xReel,dim);
			//////system("pause");
			itStart=itNew;
			cpt++;
		}
	}
	cout<< " fini creation de target nb de points "<<cpt<<endl;

}

void ViabiMicroMacro::initialiseTargetHJB_DD()
{
	/*
	 *  cette fonction initialise la base de données pour . Elle
	 *  ajoute dnas la base les premiers points pour lesquels la fonction valeur
	 *
	 *   est réelle. Au début de l'algorithme de bassin de capture
	 *    seuls le spoints de la  cible  ont une fonction valeur réelle
	 *
	 */


	double c;


	unsigned long long int dim=grid->dim;

	//ouverture de la base


	unsigned long long int  * x=new unsigned long long int [dim];

	double * xReel =new double[dim];



	int  totalPointsC=0;  // nombre de points  de l'espace des commandes




	int totalPointsX=grid->getNbTotalPoints();

	imagePoint currentPoint;

	currentImagePointsList.maxNum=0;
	currentImagePointsList.minNum=0;
	currentImagePointsList.pointsList= tempPointsList1;;

	unsigned long long int pos ;

	list<imagePoint>::iterator itStart=currentImagePointsList.pointsList->begin(),
			itNew;
	int cpt=0;

	cout<< " total points de grid "<<totalPointsX<<endl;
	/*
	 *  on parcourt  tous les points de l'espace discret  fini
	 *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
	 */
	for( pos=0;pos<(unsigned long long int)totalPointsX;pos++)
	{

		//cout<< " pos= "<<pos<<endl;
		/*!
		 * le compteur pos  ets l'unique numéro entier du point en cours
		 * dans la numérotation alphabétique : on parcourt axe par axe
		 */

		/*!
		 * on restitue les  coordonnées netiéres  du point é partir de son numéro
		 * ainsi que ses coordonnées réelles
		 */
		grid->numToIntAndDoubleCoords(pos,x,xReel);

		c=max( (*(dynsys->target_fd))(x),(*(dynsys->constraintsX_fd))(x)) ;

		vTab[pos]=c;
		if(c<PLUS_INF)
		{
			totalPointsC++;
			currentPoint.minVal=c;
			currentPoint.PointNum=pos;

			addDataToPointsList(&itStart, currentPoint, &itNew);
			//cout<<" fini ajout point  \n";
			// printVector(x,dim);
			//////system("pause");
			itStart=itNew;
			cpt++;
		}
	}
	cout<< " fini creation de target nb de points "<<cpt<<endl;

}

void ViabiMicroMacro::loadViableSets()
{

	algoViabiParams * avp=modelParams->getAlgoParameters();

	ostringstream os;
	string fileName;
	// on charge dans la mémoire l'ensemble calculé et enregistré
	// correspondant au dernier raffinement
	os<<"../OUTPUT/"<<filePrefix<<"-valFunc.dat";
	fileName=os.str();
	os.str("");
	grid->loadSet(fileName);


	if(avp->SAVE_PROJECTION)
	{
		os<<"../OUTPUT/"<<filePrefix<<"-proj"<<".dat";
		fileName=os.str();
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
	algoViabiParams * avp=modelParams->getAlgoParameters();
	int refine = avp->GRID_REFINMENTS_NUMBER;

	ostringstream os;
	string fileName;
	if(dynsys->getDynType() == DD)
	{
		if(avp->SAVE_SUBLEVEL){
			os<<"../OUTPUT/"<<filePrefix<<"-subLevel.dat";
			fileName=os.str();
			os.str("");
			grid->saveSubLevelset_DD(avp->LEVEL, fileName);
		}

		if(avp->SAVE_BOUNDARY){
			os<<"../OUTPUT/"<<filePrefix<<"-valFunc.dat";
			fileName=os.str();
			os.str("");
			grid->saveValOnGrid_DD(fileName);
		}

		if(avp->SAVE_PROJECTION)
		{
			os<<"../OUTPUT/"<<filePrefix<<"-proj"<<".dat";
			fileName=os.str();
			os.str("");
			/*
			 *  calcul et sauvegarde  de la projection du  noyau
			 */
			grid->saveProjetion(fileName, avp->PROJECTION);
		}
	}
	else
	{
		if(avp->SAVE_SUBLEVEL){
			os<<"../OUTPUT/"<<filePrefix<<"-subLevel.dat";
			fileName=os.str();
			os.str("");
			grid->saveSubLevelset(avp->LEVEL, fileName);
		}

		if(avp->SAVE_BOUNDARY){
			os<<"../OUTPUT/"<<filePrefix<<"-valFunc.dat";
			fileName=os.str();
			os.str("");
			grid->saveValOnGrid(fileName);
		}

		if(avp->SAVE_PROJECTION)
		{
			os<<"../OUTPUT/"<<filePrefix<<"-proj"<<".dat";
			fileName=os.str();
			os.str("");
			/*
			 *  calcul et sauvegarde  de la projection du  noyau
			 */
			grid->saveProjetion(fileName, avp->PROJECTION);
		}

		if(avp->SAVE_SLICE_BOUND)
				{
					os<<"../OUTPUT/"<<filePrefix<<"-slice"<<".dat";
					fileName=os.str();
					os.str("");
					/*
					 *  calcul et sauvegarde  de la projection du  noyau
					 */
					grid->saveCoupeBoundary(fileName);
				}
	}

}



void ViabiMicroMacro::viabKerValFunc()
{
	unsigned long long int iCoords[dim];
	unsigned long long int cptChanged=20, nbIter=0;
	unsigned long long int pos ;
	unsigned long long int nbCellsTotal=grid->getNbTotalPoints();
	unsigned long long int compteComm, cellNum, numControl, cu;
	unsigned long long int nbCTotal=dynsys->getTotalNbPointsC();

	double rCoords[dim];

	double rho;

	uintPair image;

	int totalPointsX=grid->getNbTotalPoints();

	long long  int * indicesDecalCell=grid->getIndicesDecalCell();
	int  nbPointsCube=(int) pow(2.0,dim);


	double ** controlCoords=dynsys->getControlCoords();
	double minValCell, valAtPos;

	while((cptChanged>10) & (nbIter<15000))
	{
		cptChanged=0;

		//cout<< "  curseur en place \n";
		/*
		 *  on parcourt  tous les points de l'espace discret  fini
		 *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
		 */
		for( pos=0;pos<(unsigned long long int)totalPointsX;pos++)
		{
			double tempVal, tempL, tempM;
			valAtPos = vTab[pos];
			if(valAtPos < PLUS_INF)
			{
				grid->numToIntAndDoubleCoords(pos,iCoords,rCoords);
				rho=dynsys->calculRho_local(rCoords);
				bool testNonVide = false;
				//cout<< " pos = "<< pos<< " taille image discrete du point "<< pointDI_DD.nbImagePoints<<endl;


				this->computeDiscreteImageOfPoint(pos);
				compteComm=0;
				minValCell=PLUS_INF;

				while( (int)compteComm<pointDI.nbImageCells && !testNonVide)
				{
					cellNum=pointDI.tabImageCells[compteComm];
					if(cellNum< nbCellsTotal)
					{

						for(int iCell=0;iCell<nbPointsCube-1;iCell++)
						{
							double imageVal = vTab[cellNum+indicesDecalCell[iCell]];
							//double imageVal = vTab[cellNum];

							for(int j=pointDI.tabCellEntrees[compteComm];j<pointDI.tabCellEntrees[compteComm+1] && !testNonVide;j++)
							{
								numControl=pointDI.tabImageControls[j];
								tempL=dynsys->lFunc(rCoords, controlCoords[numControl]);
								tempM=dynsys->mFunc(rCoords, controlCoords[numControl]);
								if(tempL < PLUS_INF && tempM < PLUS_INF)
								{
									tempVal=(imageVal+rho*tempL)/(1-rho*tempM);
									testNonVide = tempVal < valAtPos;
									minValCell=min(minValCell, tempVal);
									if(false)//(rCoords[1] >= 0.8)
									{
										cout<< " pos ";
										printVector(rCoords,dim);
										cout<< " cellNum = "<< cellNum<<endl;
										cout<< " imageVal = "<<imageVal<< " currentVal = " << vTab[pos]<< endl;
										cout<< "  lTemp = "<< tempL << " Mtemp = " << tempM<< " rho = "<< rho << " tempVal = " << tempVal<<endl;
									}
								}
								if(testNonVide) break;
							}
							if(testNonVide) break;
						}
					}
					else if(cellNum == nbCellsTotal)
					{
						for(int j=pointDI.tabCellEntrees[compteComm];j<pointDI.tabCellEntrees[compteComm+1];j++)
						{
							numControl=pointDI.tabImageControls[j];
							(dynsys->*(dynsys->discretDynamics))(rCoords, controlCoords[numControl], doubleVect1, 1.0);
							//cout<< " sortie autoris�e pour control num "<<numControl << "  image "<<endl;
							//printVector(doubleVect1, dim);
							tempL=dynsys->lFunc(rCoords, controlCoords[numControl]);
							tempM=dynsys->mFunc(rCoords, controlCoords[numControl]);

							double tempV = dynsys->constraintsX(doubleVect1);
							//cout<< " value  estimee dans ce point "<< tempV << endl;
							tempVal=(tempV+rho*tempL)/(1-rho*tempM);
							testNonVide = tempVal < valAtPos;
							minValCell=min(minValCell, tempVal);
							if(testNonVide) break;
						}
					}

					compteComm++;
				}

				if(vTab[pos]< minValCell)
				{
					cptChanged++;

				}
				vTab[pos]=max(vTab[pos], minValCell);
			}
		}

		nbIter++;
		cout<< " iteration  "<<nbIter<<"  parcours de base   fini  nb  de valeurs changées  est "<<cptChanged<<endl;
		//
	}
cout<< " calculFini. Sauvegarde\n";
	saveValFunctions();
}


void ViabiMicroMacro::viabKerValFuncOld()
{
	unsigned long long int iCoords[dim];
	unsigned long long int cptChanged=10, nbIter=0;
	unsigned long long int pos ;
	unsigned long long int nbCellsTotal=grid->getNbTotalPoints();
	unsigned long long int compteComm, cellNum, numControl;

	double rCoords[dim];

	double rho;

	uintPair image;

	int totalPointsX=grid->getNbTotalPoints();

	long long  int * indicesDecalCell=grid->getIndicesDecalCell();
	int  nbPointsCube=(int) pow(2.0,dim);


	double ** controlCoords=dynsys->getControlCoords();
	double minValCell;

	while((cptChanged>0) & (nbIter<15000))
	{
		cptChanged=0;

		//cout<< "  curseur en place \n";
		/*
		 *  on parcourt  tous les points de l'espace discret  fini
		 *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
		 */
		for( pos=0;pos<(unsigned long long int)totalPointsX;pos++)
		{
			double tempVal, tempL, tempM;

			grid->numToIntAndDoubleCoords(pos,iCoords,rCoords);
			rho=dynsys->calculRho_local(rCoords);
			this->computeDiscreteImageOfPoint(pos);
			compteComm=0;
			minValCell=PLUS_INF;

			while( (int)compteComm<pointDI.nbImageCells)
			{
				cellNum=pointDI.tabImageCells[compteComm];
				if(cellNum< nbCellsTotal)
				{

					//for(int iCell=0;iCell<nbPointsCube-1;iCell++)
					{
						//double imageVal = vTab[cellNum+indicesDecalCell[iCell]];
						double imageVal = vTab[cellNum];

						for(int j=pointDI.tabCellEntrees[compteComm];j<pointDI.tabCellEntrees[compteComm+1];j++)
						{
							numControl=pointDI.tabImageControls[j];
							tempL=dynsys->lFunc(rCoords, controlCoords[numControl]);
							tempM=dynsys->mFunc(rCoords, controlCoords[numControl]);
							if(tempL < PLUS_INF && tempM < PLUS_INF)
							{
								tempVal=(imageVal+rho*tempL)/(1-rho*tempM);
								minValCell=min(minValCell, tempVal);
								if(rCoords[1] >= 0.8)
								{
									cout<< " pos ";
									printVector(rCoords,dim);
									cout<< " cellNum = "<< cellNum<<endl;
									cout<< " imageVal = "<<imageVal<< " currentVal = " << vTab[pos]<< endl;
									cout<< "  lTemp = "<< tempL << " Mtemp = " << tempM<< " rho = "<< rho << " tempVal = " << tempVal<<endl;
								}
							}
						}
					}
				}
				else if(cellNum == nbCellsTotal)
				{
					for(int j=pointDI.tabCellEntrees[compteComm];j<pointDI.tabCellEntrees[compteComm+1];j++)
					{
						numControl=pointDI.tabImageControls[j];
						(dynsys->*(dynsys->discretDynamics))(rCoords, controlCoords[numControl], doubleVect1, 1.0);
						//cout<< " sortie autoris�e pour control num "<<numControl << "  image "<<endl;
						//printVector(doubleVect1, dim);
						tempL=dynsys->lFunc(rCoords, controlCoords[numControl]);
						tempM=dynsys->mFunc(rCoords, controlCoords[numControl]);

						double tempV = dynsys->constraintsX(doubleVect1);
						//cout<< " value  estimee dans ce point "<< tempV << endl;
						tempVal=(tempV+rho*tempL)/(1-rho*tempM);
						minValCell=min(minValCell, tempVal);
					}
				}

				compteComm++;
			}

			if(vTab[pos]< minValCell)
			{
				cptChanged++;

			}
			vTab[pos]=max(vTab[pos], minValCell);

		}

		nbIter++;
		cout<< " iteration  "<<nbIter<<"  parcours de base   fini  nb  de valeurs changées  est "<<cptChanged<<endl;
		//
	}

	saveValFunctions();
}

void ViabiMicroMacro::viabKerValFunc_DD()
{
	unsigned long long int iCoords[dim];
	unsigned long long int cptChanged=10, nbIter=0;
	unsigned long long int pos ;
	unsigned long long int compteComm, cellNum, numControl;
	unsigned long long int nbPointsTotal=grid->getNbTotalPoints();

	unsigned long long int ** controlCoords=dynsys->getControlIntCoords();

	double rCoords[dim];

	double rho;

	uintPair image;

	unsigned long long int nbCTotal=dynsys->getTotalNbPointsC();
	unsigned long long int imagePos;


	// cout<<  " compute DI  pos = "<<num<< " juste apres num to int double coords\n";
	// printVector(intPointCoords, dim);
	unsigned long long int cu;



	long long  int * indicesDecalCell=grid->getIndicesDecalCell();
	int  nbPointsCube=(int) pow(2.0,dim);

	double minValCell;

	while((cptChanged>0) & (nbIter<15000))
	{
		cptChanged=0;
		cout<< " debut de l'iter "<<nbIter<<endl;
		//cout<< "  curseur en place \n";
		/*
		 *  on parcourt  tous les points de l'espace discret  fini
		 *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
		 */
		for( pos=0;pos<(unsigned long long int)nbPointsTotal;pos++)
		{
			double currentVal = vTab[pos];
			grid->numToIntAndDoubleCoords(pos,intPointCoords,doublePointCoords);
			if(currentVal<PLUS_INF)
			{
				double tempVal, tempL, tempM;

				//this->computeDiscreteImageOfPoint_DD(pos);
				compteComm=0;
				minValCell=PLUS_INF;
				unsigned long long int numPoint;
				unsigned long long int cu = 0;
				bool testNonVide = false;
				//cout<< " pos = "<< pos<< " taille image discrete du point "<< pointDI_DD.nbImagePoints<<endl;
				while( !testNonVide && ( cu<nbCTotal))
				{
					if(dynsys->constraintsXU_fd(intPointCoords,controlCoords[cu])<PLUS_INF)
					{
						// printf("  contraintes sur U et X  ok\n");
						//  printf( " x= ");
						//  printVector(doublePointCoords,dim);

						dynsys->dynamics_fd(intPointCoords, controlCoords[cu], intVect1);

						// cout<< " retrour dnamique discrete ";
						// printVector(intVect1, dim);
						//(dynsys->*discretDynamics)(doublePointCoords, controlCoords[cu], doubleVect);
						if(grid->isPointInGrid_fd(intVect1))
						{
							grid->intCoordsToNum(intVect1, &imagePos);
							//////printf(" le point est das la grlle\n");
							/*!
							 * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
							 */
							if(vTab[imagePos]<PLUS_INF)
							{
								tempL=dynsys->lFunc_fd(intPointCoords, controlCoords[cu]);
								//cout<< "  lfunc = "<<tempL << " mufunc "<<tempM<<endl;
								tempM=dynsys->muFunc_fd(intPointCoords, controlCoords[cu]);

								tempVal=max((vTab[imagePos]-tempL), tempM);
								//cout<< "  lfunc = "<<tempL << " mufunc "<<tempM<< " tempVal = "<<tempVal<<endl;
								minValCell=min(minValCell, tempVal);
								testNonVide = (tempVal <= currentVal);
							}
						}
					}
					cu++;
				}
			}
			//	cout<< "vTab[pos] "<< vTab[pos] << " minValCell = "<<minValCell<<endl;
			if(currentVal< minValCell)
			{
				cptChanged++;
			}
			vTab[pos]=max(vTab[pos], minValCell);
		}


		nbIter++;
		cout<< " iteration  "<<nbIter<<"  parcours de base   fini  nb  de valeurs changées  est "<<cptChanged<<endl;
		//
	}

	saveValFunctions();
}

void ViabiMicroMacro::CaptureBasin()
{
	cout<< " Debut de calcul : bassin de capture "<< endl;

	/*!
	 * \var nbNewPoints : nb de nouveaux points ajoutés é l'étape n
	 */
	int nbNewPoints=1;


	int iter=0;

	/*!
	 * On calcule la premiére itération, en tanant compte de la convexification de la dynamique sur la cible
	 * On appelle pour cela la fonction computeConvexifiedImage().
	 */
	(this->*computeFirstConvexifiedImage)(iter);

	//	 this->showCurrentImageList();
	////////system("pause");
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
	if(!(dynsys->getDynType()==DD)   )
	{
		nbNewPoints=(this->*addNewPointsToSet)();
		cout<< " points ajoutes =  "<<nbNewPoints<<endl;
	}

	iter++;

	/*!
	 * Tant qu'il y a de nouveaux points ajoutés on répéte les opérations suivantes.
	 */
	while( (nbNewPoints>0))
	{
		cout<<"  nbNewPoints=  "<<nbNewPoints<<endl;

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
		nbNewPoints=(this->*addNewPointsToSet)();
		cout<< " points ajoutes =  apres iter  "<<iter<< " "<<nbNewPoints<<endl;
		iter++;
	}

	saveValFunctions();
}




/*
void ViabiMicroMacro::captBasinEpi_omp()
{
  cout<< " minTimeEpiLocalRho \n";
  !
 * \var nbNewPoints : nb de nouveaux points ajoutés é l'étape n

  int nbNewPoints=1;


  this->initCnIndices();



  int iter=0;
  !
 * On calcule la premiére itération, en tanant compte de la convexification de la dynamique sur la cible
 * On appelle pour cela la fonction computeConvexifiedImage().

  (this->*computeFirstConvexifiedImage_omp)(iter);
  //     this->showCurrentImageList();
  ////////system("pause");
  !
 * -A partir de la liste des mailles on crée une liste ordonnée de points représentant l'image é l'aide de la fonction
 * createPointsList().

  (this->*createCurrentPointsList)();
  !
 * - A partir de la liste des points représentant l'image  trois bases  de données sont alimentées: La base
 * principale, associée é la grille  et représentant la fonction valeur, la base de rétroaction optimale et la base de
 *  rétro-action viable. On appelle pour cela la fonction addNewPoints().

  nbNewPoints=(this->*addNewPointsToSet)();
  cout<< " points ajoutes =  "<<nbNewPoints<<endl;

  iter++;

  !
 * Tant qu'il y a de nouveaux points ajoutés on répéte les opérations suivantes.

  while( (nbNewPoints>0))
    {
    cout<<"  nbNewPoints=  "<<nbNewPoints<<endl;

    !
 * -Calcul de l'image de \f^C_{n}\setminus C_n\f$  par la fonction computeCurrentImageLocalRho(): l'image est enregistrée en émoire vive sous forme de liste
 * de références de mailles dans lesquelles arrive au moins une évolution. Chaque référence de maille contient
 * des informations sur tous les antécédants de cette maille ainsi  que la valeur minimale de
 * temps. Dans cette version oé \f$ \rho\f$ est global, la fonction valeur prend la méme valeur
 * é chaque étape : \f$ \rho \cdot n \f$.

    (this->*computeCurrentImage)(iter);
    //this->showCurrentImageList();

    !
 * -A partir de la liste des mailles on crée une liste ordonnée de points représentant l'image en appelant la fonction createPointsList().
 *  Comme pour le mailles, chaque référence de point regroupe les informatons (regroupées é partir de différentes
 * mailles dont est vertex)  sur les antécédents  de ce point. Le but de la création de cette liste est d'éliminer
 * les doublons afin de minimiser les accés é la base de données

    (this->*createCurrentPointsList)();
    //this->showCurrentImagePointsList();
    !
 * - A partir de la liste des points représentant l'image  trois bases  de données sont alimentées: La base
 * principale, associée é la grille  et représentant la fonction valeur, la base de rétroaction optimale et la base de
 *  rétro-action viable. On appelle ici la fonction addNewPoints().

    nbNewPoints=(this->*addNewPointsToSet)();
    // cout<< " points ajoutes =  "<<nbNewPoints<<endl;
    iter++;
    }
}
 */


int  ViabiMicroMacro::findOptiControl_tmin(double *currentPos,
		double &dt,
		int nbStepIter,
		double stepCoeff,
		double *resPos,
		int & nbViabVoisins,
		bool &succes )
{



	int  nbPointsCube=(int) pow(2.0,dim);//pow(2.0, dim);
	/*
	 * tableaux temporaires pour récupérer les indices du point dans la
	 * grille
	 */

	unsigned long long int   testI[dim];
	double testV[dim];
	/*
	 * tableau de coordonnées  définissant
	 * la grille de contrôles
	 */
	double ** controlCoords=dynsys->getControlCoords();
	/*
	 * nombre total de points de contrôle
	 */
	unsigned long long int nbCTotal=dynsys->getTotalNbPointsC();

	/*
	 * indices de déplacement pour parcourir les sommets d'une malle à partir
	 * du sommet inf
	 */
	long long  int * indicesDecalCell=grid->getIndicesDecalCell();
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
	int cellNum=0;

	int posTemp;

	/*
	 * tests de validité de point initial
	 */
	bool testNonVide=false;
	// int cptOK=0;

	cout<< "  debut find control optim tmin : current pos = ";
	for(int i=0;i<(int)dim;i++)
	{
		xCoordsDouble[i]=currentPos[i];
		cout<< " "<< xCoordsDouble[i]<<endl;
	}

	int maxnbViabPoints;
	int bestCu;
	rho=dt;
	dynsys->setRho(rho);
	cout<< " find control :  rho= "<<rho<<endl;
	double       currentVal=PLUS_INF;
	cellNum=grid->localizePoint(xCoordsDouble);
	for(int ii=0;ii<nbPointsCube;ii++)
	{
		posTemp= cellNum+indicesDecalCell[ii];
		grid->numToIntAndDoubleCoords( posTemp ,testI,testV);
		currentVal=min( vTab[posTemp],currentVal );
		cout<< " posTemp = "<<posTemp << " currentVal = "<<currentVal <<endl;

	}
	cout<< " val of current point is !!!!!!!! "<<currentVal<<endl;
	/*
	 * on parcours tous les contrôles
	 */
	maxnbViabPoints=0;
	bestCu=0;
	double minVal=PLUS_INF, minValCell;
	int iter=0;

	testNonVide=false;
	while(iter<nbStepIter  && !testNonVide)
	{
		//cout<< "find control  iteration  "<<iter<<" rho= "<<rho<<endl;
		//cout<< " au debut on a testNonVide "<< testNonVide<<endl;
		unsigned long long int cu=0;
		testNonVide=false;
		dt=rho;
		minVal=PLUS_INF;
		while(cu<nbCTotal  )
		{
			/*
			 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
			 * au point en cours
			 */
			if(dynsys->constraintsXU(xCoordsDouble,controlCoords[cu])<PLUS_INF)
			{
				(dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[cu], doubleVect1, rho);
				if(grid->isPointInGrid(doubleVect1))
				{
					/*
					 *  le sucesseur est dans la grille de calcul
					 */
					if(dynsys->constraintsX(doubleVect1)<PLUS_INF)
					{
//								 cout<< " image du point ";
//						            for(int k=0;k<dim;k++)
//						              {
//						              cout<< " "<<doubleVect1[k];
//						              }
						/*
						 * le successeur vérifie les contraintes
						 * On identifie la maille où il se trouve
						 */
						cellNum=grid->localizePoint(doubleVect1);
						 //cout<< " num cellule "<<cellNum<<endl;

						/*
						 * On parcours les sommets de la maille
						 * autour du sucesseur pour voir s'il y a des
						 * points viables
						 */
						minValCell=PLUS_INF;
						for(int ii=0;ii<nbPointsCube;ii++)
						{
							posTemp= cellNum+indicesDecalCell[ii];
							minValCell=min( vTab[posTemp],minValCell );
						}
 //cout<< " current u "<<cu<< " minVal "<<minVal<< " minVal du cell " << minValCell << endl;
						if(minValCell<=minVal+1e-6)
						{

							minVal=minValCell;
							bestCu=cu;
//cout<< " best cu  num = " << bestCu << " coords "<< controlCoords[cu][0] << endl;
						}
					}
				}
			}
			cu++;
		}//fin de parcours de tous les contrôles
		// la boucle s'arête ici u premier contrôle
		// qui donne un successeur viable
		cout<< " min val = "<<minVal<< " currentVal = "<< currentVal<<endl;
		testNonVide=(minVal<currentVal);
		if(testNonVide) cout<< " TROUVE OK \n";
		iter++;
		rho=rho*stepCoeff;
		dynsys->setRho(rho);

	}

	succes=testNonVide;
	nbViabVoisins=maxnbViabPoints;
	 cout<< "  fin de recherche de controle optima on a bestCu nu m "<<bestCu<<endl;
	(dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[bestCu], doubleVect1, rho);
	for(int i=0;i<(int)dim;i++)
	{
		resPos[i]=doubleVect1[i];
	}
	dt=rho;
	return bestCu;
}


double ViabiMicroMacro::computeOptimalCaptTrajectory(double *initPosition, string fileName, bool &succes)
{
	if(computeTmin)
		return  this->computeOptimalTrajectory_tmin(initPosition,  fileName, succes);
	else
		return this->computeOptimalTrajectory_Lmin(initPosition,  fileName, succes);
}



double ViabiMicroMacro::computeOptimalTrajectory_tmin(double *initPosition, string fileName, bool &succes)
{

	int  nbPointsCube=(int) pow(2.0,dim);//pow(2.0, dim);

	/*
	 * tableau de coordonnées  définissant
	 * la grille de contrôles
	 */
	double ** controlCoords=dynsys->getControlCoords();

	/*
	 * indices de déplacement pour parcourir les sommets d'une malle à partir
	 * du sommet inf
	 */
	long long  int * indicesDecalCell=grid->getIndicesDecalCell();
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
	double T=dynsys->getTimeHorizon();
	/*
	 * structures accumulables dansune liste
	 */
	valarray<double> newTrajPoint(dim+1);
	valarray<double> trajControlCoords(dimC);

	int posTemp;

	cout<< " calcul de traj TMIN a partir de coords \n";
	cout<< " Postion initiale = ";

	for(int l1=0;l1<(int)dim;l1++)
	{
		cout<< " "<<initPosition[l1];
	}
	cout<< " \n";
	/*
	 * tests de validité de point initial
	 */
	bool testNonVide=false;
	double minCellVal=PLUS_INF;
	double time=0.0;



	if(grid->isPointInGrid(initPosition))
	{
		if(dynsys->constraintsX(initPosition)<PLUS_INF)
		{
			cellNum=grid->localizePoint(initPosition);
			minCellVal=PLUS_INF;
			cout<< " calcul de la valeur dans la cellule du point initial num "<<cellNum<<endl;
			for(int ii=0;ii<nbPointsCube;ii++  )
			{
				posTemp= cellNum+indicesDecalCell[ii];
				minCellVal=min(minCellVal, vTab[posTemp]);
				cout<< " posTemp = "<< posTemp << " celVal = "<< minCellVal<< endl;
			}
			cout<< " val of initial point is "<< minCellVal<< endl;
			testNonVide=(minCellVal<PLUS_INF);

			if(!testNonVide)
			{
				cout<<" La position initiale sélectionnée n'appartiant pas au bassin de capture\n";
				succes=0;
			}
			else
			{
				/*
				 * la position initiale se trouve dans le noyau de viabilité
				 * on initialise le temps à 0  et recopie la pos initiale
				 * dans le coordonnées temporaires du point en cours de la trajectoire
				 */
				time=0.0;
				for(int i=0;i<(int)dim;i++)
				{
					xCoordsDouble[i]=initPosition[i];
				}
				int nbIter=0;
				/*
				 * On itère tant que le temps n'a pas dépassé l'horizon donné
				 */

				double c, realTimeStep;
				bool testviabInt=true;
				int maxnbViabPoints;
				int bestCu;
				c=dynsys->target(xCoordsDouble);
				while((time<T) && (c>=PLUS_INF) && (nbIter<=NB_MAX_TRAJ_ITER))
				{
					cout<< " point en cours ";
					for(int i=0;i<(int)dim;i++)
					{
						newTrajPoint[i]=xCoordsDouble[i];
						cout<< " "<<newTrajPoint[i];
					}
					newTrajPoint[dim]=time;
					cout<< " temps= "<<newTrajPoint[dim]<<endl;
					traj.push_back(newTrajPoint);

					rho=  dynsys->calculRho_local(xCoordsDouble);

					rho=min(rho, T-time);
					realTimeStep = 2.5*rho;
					cout<< " rho= "<<rho<<endl;

					bestCu=this->findOptiControl_tmin(xCoordsDouble, realTimeStep, 1,1.0, imageVect,maxnbViabPoints, testNonVide );


					// la boucle s'arête ici u premier contrôle
					// qui donne un successeur viable

					if(testNonVide)
					{
						// contrôle viable trouvé
						// on recopie ce contrôle dans la liste et
						// le successeur devient le point  courent
						time+=realTimeStep;
						for(int dc=0;dc<(int)dimC;dc++)
						{
							trajControlCoords[dc]=controlCoords[bestCu][dc];
						}
						trajC.push_back(trajControlCoords);
						for(int i=0;i<(int)dim;i++)
						{
							xCoordsDouble[i]=imageVect[i];
						}

					}
					else
					{
						realTimeStep = 2.0*rho;
						cout<< " rho= "<<rho<<endl;

						bestCu=this->findOptiControl_tmin(xCoordsDouble, realTimeStep, 5,0.7, imageVect,maxnbViabPoints, testNonVide );

						if(testNonVide)
						{
							// contrôle viable trouvé
							// on recopie ce contrôle dans la liste et
							// le successeur devient le point  courent
							time+=realTimeStep;
							for(int dc=0;dc<(int)dimC;dc++)
							{
								trajControlCoords[dc]=controlCoords[bestCu][dc];
							}
							trajC.push_back(trajControlCoords);
							for(int i=0;i<(int)dim;i++)
							{
								xCoordsDouble[i]=imageVect[i];
							}

						}
						else
						{


							cout<<"   Echec! Sortie de l'ensemble viable \n";
							break;

						}

					}
					c=dynsys->target(xCoordsDouble);
					cout<<"   valeur cible : "<<c<<endl;
					nbIter++;
				}
				if(c< PLUS_INF)
				{
					succes=1;
				}
			}//fin de else (reconstruction de trajectoire)

		}
		else
		{
			printf(" Point initial hors de l'ensemble de contraintes. Arret\n");
			succes=0;
		}
	}
	else
	{
		printf(" Point initial hors de grille de calcul. Arret\n");
		succes=0;
	}

	// if(succes)
	{

		printf(" trajectoire trouvée. Enregistrement\n");

		FILE * fi;
		fi = fopen( fileName.c_str(),"w");
		if(fi==NULL){
			printf("** error: impossible to open the file %s.\n", fileName.c_str());

		}
		else{
			list<valarray<double> >::iterator it=traj.begin();
			list<valarray<double> >::iterator itc=trajC.end();
			itc--;
			trajC.push_back(*itc);
			itc=trajC.begin();

			while(it!=traj.end())
			{

				for(int l1=0;l1<(int)dim;l1++)
				{
					fprintf(fi,  "%15.8f " ,   (*it)[l1]);
					cout<< " "<<  (*it)[l1];
				}
				fprintf( fi, "%15.8f " ,  (*it)[dim]);
				cout<< " "<<(*it)[dim];
				for(int dc=0;dc<(int)dimC;dc++)
				{
					fprintf( fi, "%15.8f " ,   (*itc)[dc]);
					cout<< " "<<(*itc)[dc]<<endl;
				}
				fprintf( fi, "\n" );
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

int  ViabiMicroMacro::findOptiControl_Lmin(double budget, double *currentPos,
		double &dt,
		int nbStepIter,
		double stepCoeff,
		double *resPos,
		double & newBudget,
		bool &succes )
{
	int  nbPointsCube=(int) pow(2.0,dim);//pow(2.0, dim);
	/*
	 * tableaux temporaires pour récupérer les indices du point dans la
	 * grille
	 */

	unsigned long long int   testI[dim];
	double testV[dim];
	/*
	 * tableau de coordonnées  définissant
	 * la grille de contrôles
	 */
	double ** controlCoords=dynsys->getControlCoords();
	/*
	 * nombre total de points de contrôle
	 */
	unsigned long long int nbCTotal=dynsys->getTotalNbPointsC();

	/*
	 * indices de déplacement pour parcourir les sommets d'une malle à partir
	 * du sommet inf
	 */
	long long  int * indicesDecalCell=grid->getIndicesDecalCell();
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
	bool testNonVide=false;

	for(int i=0;i<(int)dim;i++)
	{
		xCoordsDouble[i]=currentPos[i];
	}

	int bestCu;
	rho=dt;
	dynsys->setRho(rho);
	// cout<< " find control :  rho= "<<rho<<endl;

	/*
	 * on parcours tous les contrôles
	 */

	bestCu=0;
	double minVal=PLUS_INF, minValCell;
	int iter=0;
	testNonVide=false;
	cellNum=grid->localizePoint(currentPos);
	/*
	 * On parcourt les sommets de la maille
	 * autour du sucesseur pour voir s'il y a des
	 * points viables
	 */
	double       currentVal=PLUS_INF;
	for(int ii=0;ii<nbPointsCube;ii++)
	{
		posTemp= cellNum+indicesDecalCell[ii];
		grid->numToIntAndDoubleCoords( posTemp ,testI,testV);
		currentVal=min( vTab[posTemp],currentVal );
	}

	testNonVide=false;
	succes=false;
	while(iter<nbStepIter  && ! testNonVide)
	{
		unsigned long long  int cu=0;
		dt=rho;
		minVal=PLUS_INF;
		while(cu<nbCTotal  )
		{
			/*
			 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
			 * au point en cours
			 */
			if(dynsys->constraintsXU(xCoordsDouble,controlCoords[cu])<PLUS_INF)
			{
				(dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[cu], doubleVect1, rho);
				if(grid->isPointInGrid(doubleVect1))
				{
					/*
					 *  le sucesseur est dans la grille de calcul
					 */
					if(dynsys->constraintsX(doubleVect1)<PLUS_INF)
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

						double tempL=dynsys->lFunc(xCoordsDouble, controlCoords[cu]);

						double tempL1=dynsys->lFunc(doubleVect1, controlCoords[cu]);

						double newVal= budget-rho*0.5*(tempL+tempL1);
						//       cout<< " budget  = "<<budget<< "L= "<<tempL<< "L1 " <<tempL1<<  " rho= "<<rho<<endl;

						cellNum=grid->localizePoint(doubleVect1);
						// cout<< " num cellule "<<cellNum<<endl;
						//   iCell=0;
						/*
						 * On parcours les sommets de la maille
						 * autour du sucesseur pour voir s'il y a des
						 * points viables
						 */
						minValCell=PLUS_INF;
						for(int ii=0;ii<nbPointsCube;ii++)
						{
							posTemp= cellNum+indicesDecalCell[ii];
							grid->numToIntAndDoubleCoords( posTemp ,testI,testV);
							//    cout<< " budget  = "<<budget<< " new val= "<<newVal<< " val point cellule image " <<vTab[posTemp]<< endl;
							minValCell=min( vTab[posTemp],minValCell );
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

						if(minValCell<minVal)
						{
							minVal=minValCell;
							newBudget=newVal;
							bestCu=cu;
							cout<< " min val = "<<minVal<< " current val= "<<currentVal<< " test= "<<testNonVide<<endl;

						}
						testNonVide=(minVal<=currentVal);
						//   testNonVide=(minVal<PLUS_INF);
					}
				}
			}
			cu++;
		}//fin de parcours de tous les contrôles
		// la boucle s'arête ici u premier contrôle
		// qui donne un successeur viable

		// cout<<   " fini parcours de controles on a test interieur = "<<testviabInt<<
		//     " test non vide "<<testNonVide<< " maxnbViabPoints =  "<<maxnbViabPoints<< " bes c u= "<<bestCu<<endl;
		iter++;
		rho=rho*stepCoeff;
		dynsys->setRho(rho);
		cout<< "find control  iteration  "<<iter<<" rho= "<<rho<<endl;
		succes=testNonVide;
		cout<< " succes = "<<succes<<endl;
	}


	(dynsys->*(dynsys->discretDynamics))(xCoordsDouble, controlCoords[bestCu], doubleVect1, rho);
	for(int i=0;i<(int)dim;i++)
	{
		resPos[i]=doubleVect1[i];
	}
	dt=rho;
	return bestCu;
}

unsigned long long int  ViabiMicroMacro::findViabControl_DD(double budget, unsigned long long int *currentPos,
		unsigned long long int *resPos,
		double & newBudget,
		bool &succes )
{
	int  nbPointsCube=(int) pow(2.0,dim);//pow(2.0, dim);
	unsigned long long int nbPointsTotal=grid->getNbTotalPoints();

	unsigned long long int ** controlCoords=dynsys->getControlIntCoords();
	/*
	 * tableaux temporaires pour récupérer les indices du point dans la
	 * grille
	 */

	unsigned long long int   testI[dim];
	double testV[dim];

	/*
	 * nombre total de points de contrôle
	 */
	unsigned long long int nbCTotal=dynsys->getTotalNbPointsC();

	/*
	 * indices de déplacement pour parcourir les sommets d'une malle à partir
	 * du sommet inf
	 */
	long long  int * indicesDecalCell=grid->getIndicesDecalCell();
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
	bool testNonVide=false;

	for(int i=0;i<(int)dim;i++)
	{
		xCoords[i]=currentPos[i];
	}

	int bestCu;

	/*
	 * on parcours tous les contrôles
	 */

	bestCu=0;
	double minVal=PLUS_INF, minValCell;
	int iter=0;
	testNonVide=false;
	/*
	 * On parcourt les sommets de la maille
	 * autour du sucesseur pour voir s'il y a des
	 * points viables
	 */
	double       currentVal=PLUS_INF;

	grid->intCoordsToNum(xCoords, &posTemp);
	currentVal= vTab[posTemp];

	testNonVide=false;
	succes=false;

	unsigned long long  int cu=0;

	minVal=PLUS_INF;
	while(cu<nbCTotal  && !testNonVide)
	{
		/*
		 * on ne choisit que ceux qui vérifies les éventuelles contrainets mixtes
		 * au point en cours
		 */
		if(dynsys->constraintsXU_fd(xCoords,controlCoords[cu])<PLUS_INF)
		{
			dynsys->dynamics_fd(xCoords, controlCoords[cu], intVect1);

			// cout<< " retrour dnamique discrete ";
			// printVector(intVect1, dim);
			//(dynsys->*discretDynamics)(doublePointCoords, controlCoords[cu], doubleVect);
			if(grid->isPointInGrid_fd(intVect1))
			{
				// printf( "   le point est das la grlle\n " );

				//////printf(" le point est das la grlle\n");
				/*!
				 * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
				 */
				if(dynsys->constraintsX_fd(intVect1)<PLUS_INF)
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

					double tempL=dynsys->lFunc_fd(xCoords, controlCoords[cu]);

					double newVal= budget+tempL;
					//       cout<< " budget  = "<<budget<< "L= "<<tempL<< "L1 " <<tempL1<<  " rho= "<<rho<<endl;

					grid->intCoordsToNum(intVect1, &imagePos);
					testNonVide=(newVal >= vTab[imagePos]);
					// cout<< " num cellule "<<cellNum<<endl;
					//   iCell=0;
					/*
					 * On parcours les sommets de la maille
					 * autour du sucesseur pour voir s'il y a des
					 * points viables
					 */
					if(testNonVide)
					{
						newBudget=newVal;
						bestCu=cu;
						for(int i=0;i<(int)dim;i++)
						{
							resPos[i]=intVect1[i];
						}
						cout<< " new budget = "<<newBudget<< " image val = "<<vTab[imagePos]<< " test= "<<testNonVide<<endl;

					}
				}
			}
		}
		cu++;
	}//fin de parcours de tous les contrôles
	// la boucle s'arête ici u premier contrôle
	// qui donne un successeur viable
	succes=testNonVide;
	cout<< " succes = "<<succes<<endl;
	return bestCu;
}


double ViabiMicroMacro::computeOptimalTrajectory_Lmin(double *initPosition, string fileName, bool &succes)
{
	int  nbPointsCube=(int) pow(2.0,dim);//pow(2.0, dim);
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
	double ** controlCoords=dynsys->getControlCoords();

	/*
	 * indices de déplacement pour parcourir les sommets d'une malle à partir
	 * du sommet inf
	 */
	long long  int * indicesDecalCell=grid->getIndicesDecalCell();
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
	double T=dynsys->getTimeHorizon();
	/*
	 * structures accumulables dansune liste
	 */
	valarray<double> newTrajPoint(dim+1);
	valarray<double> trajControlCoords(dimC);

	int posTemp;

	cout<< " calcul de traj a partir de coords \n";
	cout<< " Postion initiale = ";

	for(int l1=0;l1<(int)dim;l1++)
	{
		cout<< " "<<initPosition[l1];
	}
	cout<< " \n";
	/*
	 * tests de validité de point initial
	 */
	bool testNonVide=false;

	double minCellVal=PLUS_INF;

	double time=0.0;

	if(grid->isPointInGrid(initPosition))
	{
		if(dynsys->constraintsX(initPosition)<PLUS_INF)
		{
			cellNum=grid->localizePoint(initPosition);
			minCellVal=-PLUS_INF;
			for(int ii=0;ii<nbPointsCube;ii++  )
			{
				posTemp= cellNum+indicesDecalCell[ii];
				minCellVal=max(minCellVal, vTab[posTemp]);
			}
			double budget= minCellVal, newBudget;
			cout<< " value of init point  "<<budget<<endl;
			testNonVide=(minCellVal<PLUS_INF);

			if(!testNonVide)
			{
				cout<<" La position initiale sélectionnée n'appartiant pas au noyau de viabilité\n";
				succes=0;
			}
			else
			{
				/*
				 * la position initiale se trouve dans le noyau de viabilité
				 * on initialise le temps à 0  et recopie la pos initiale
				 * dans le coordonnées temporaires du point en cours de la trajectoire
				 */
				time=0.0;
				for(int i=0;i<(int)dim;i++)
				{
					xCoordsDouble[i]=initPosition[i];
				}
				int nbIter=0;
				/*
				 * On itère tant que le temps n'a pas dépassé l'horizon donné
				 */

				double c;
				int   bestCu;
				c=dynsys->target(xCoordsDouble);
				while((time<T) && (c>=PLUS_INF) && (nbIter<=NB_MAX_TRAJ_ITER))
				{
					cout<< " point en cours ";
					for(int i=0;i<(int)dim;i++)
					{
						newTrajPoint[i]=xCoordsDouble[i];
						cout<< " "<<newTrajPoint[i];
					}
					newTrajPoint[dim]=time;
					cout<< " temps= "<<newTrajPoint[dim]<<endl;
					traj.push_back(newTrajPoint);

					rho=    0.75*dynsys->calculRho_local(xCoordsDouble);

					rho=min(rho, T-time);

					cout<< " rho= "<<rho<<endl;

					bestCu=this->findOptiControl_Lmin(budget,xCoordsDouble, rho,1,1.0, imageVect,newBudget, testNonVide );
					time+=rho;

					// la boucle s'arête ici u premier contrôle
					// qui donne un successeur viable

					cout<<   " Premiere recherche de controle viable  fini parcours de controles on a  test non vide "<<testNonVide<< " bes c u= "<<bestCu<<endl;

					// contrôle viable trouvé
					// on recopie ce contrôle dans la liste et
					// le successeur devient le point  courent
					if(testNonVide)
					{

						cout<<  " image interieure tourvee \n";

						for(int dc=0;dc<(int)dimC;dc++)
						{
							trajControlCoords[dc]=controlCoords[bestCu][dc];
						}
						trajC.push_back(trajControlCoords);
						for(int i=0;i<(int)dim;i++)
						{
							xCoordsDouble[i]=imageVect[i];
						}
						budget=newBudget;
						cout<< " new budget = "<<newBudget<<endl;
					}
					else
					{
						cout<< "  recherche de optimal viable avec iterations  sur le pas  de temps\n";
						rho=   0.5*dynsys->calculRho_local(xCoordsDouble);

						rho=min(rho, T-time);

						bestCu=this->findOptiControl_Lmin(budget, xCoordsDouble, rho,5,1.1,imageVect,newBudget, testNonVide );
						if(testNonVide)
						{
							for(int dc=0;dc<dimC;dc++)
							{
								trajControlCoords[dc]=controlCoords[bestCu][dc];
							}
							trajC.push_back(trajControlCoords);
							cout<< "  coords double : ";
							for(int i=0;i<dim;i++)
							{
								xCoordsDouble[i]=imageVect[i];
								cout<< " "<<xCoordsDouble[i];
							}
							cout<<endl;
							budget=newBudget;
							cout<< " new budget = "<<newBudget<<endl;
						}

						else
						{
							cout<<"   Echec! Sortie de l'ensemble viable \n";
							break;

						}
					}
					c=dynsys->target(xCoordsDouble);
					cout<<"   valeur cible : "<<c<<endl;
					nbIter++;
				}
				if(c< PLUS_INF)
				{
					succes=1;
				}
			}//fin de else (reconstruction de trajectoire)

		}
		else
		{
			printf(" Point initial hors de l'ensemble de contraintes. Arret\n");
			succes=0;
		}
	}
	else
	{
		printf(" Point initial hors de grille de calcul. Arret\n");
		succes=0;
	}

	// if(succes)
	{

		printf(" trajectoire trouvée. Enregistrement\n");

		FILE * fi;
		fi = fopen( fileName.c_str(),"w");
		if(fi==NULL){
			printf("** error: impossible to open the file %s.\n", fileName.c_str());

		}
		else{
			list<valarray<double> >::iterator it=traj.begin();
			list<valarray<double> >::iterator itc=trajC.end();
			itc--;
			trajC.push_back(*itc);
			itc=trajC.begin();

			while(it!=traj.end())
			{

				for(int l1=0;l1<dim;l1++)
				{
					fprintf(fi,  "%15.8f " ,   (*it)[l1]);
					cout<< " "<<  (*it)[l1];
				}
				fprintf( fi, "%15.8f " ,  (*it)[dim]);
				cout<< " "<<(*it)[dim];
				for(int dc=0;dc<dimC;dc++)
				{
					fprintf( fi, "%15.8f " ,   (*itc)[dc]);
					cout<< " "<<(*itc)[dc]<<endl;
				}
				fprintf( fi, "\n" );
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



void ViabiMicroMacro::computeCurrIm_tmin( int iter)
{


	cout<< " toto calcul de l'image de Cn local opti nb points Cn = "<<currentImagePointsList.pointsList->size()<<endl;

	double t1,t2,elapsed_time;

	timeval tim;
	gettimeofday(&tim,NULL);
	t1=tim.tv_sec+(tim.tv_usec/1000000.0);


	int posX;


	double T=dynsys->getTimeHorizon();


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
	currentImageList.maxNum=-1;
	// currentImageList.minNum=PLUS_INF;
	currentImageList.minNum=1000+grid->nbTotalCells;
	double rho;

	list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList->begin(),
			itTemp;
	imageCell tempImageCell;

	itStart=this->currentImageList.cellsList.begin();

	while(!currentImagePointsList.pointsList->empty())//(itPoint!=itLastPoint)
	{

		posX=(*itPoint).PointNum;
		// cout<< " posX="<<posX<<endl;
		grid->numToIntAndDoubleCoords(posX,intPointCoords,doublePointCoords);
		// printVector(doublePointCoords, dim);
		rho=dynsys->calculRho_local(doublePointCoords);
		//  cout<< " rho= "<<rho;
		if((*itPoint).minVal+rho<=T)
		{
			/*!
			 * On calcule l'image discrète  du point
			 */
			tempImageCell.minVal=(*itPoint).minVal+rho;

			this->computeDiscreteImageOfPoint(posX);
			/*!
			 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
			 *  de cellules
			 */

			for(int i=0;i<pointDI.nbImageCells;i++)
			{
				tempImageCell.cellNum=pointDI.tabImageCells[i];
				addDataToCurrentImage(&itStart,tempImageCell, &itNew);
				//cout<< " ares ajout  de l'image d'un point la taille de la liste est "<<currentImageList.cellsList.size()<<endl;
				itStart=itNew;
			}
		}
		itPoint++;
		currentImagePointsList.pointsList->pop_front();
	}
	gettimeofday(&tim,NULL);
	t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	elapsed_time=(double)((t2-t1));

	cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;
}

void ViabiMicroMacro::computeCurrIm_Lmin( int iter)
{

	double t1,t2,elapsed_time;

	timeval tim;
	gettimeofday(&tim,NULL);
	t1=tim.tv_sec+(tim.tv_usec/1000000.0);


	int posX;

	unsigned long long int nbCellsTotal=grid->getNbTotalPoints();

	double ** controlCoords=dynsys->getControlCoords();
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
	currentImageList.maxNum=-1;
	// currentImageList.minNum=PLUS_INF;
	currentImageList.minNum=1000+grid->nbTotalCells;
	double rho, tempL;

	list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList->begin(),
			itTemp;
	imageCell tempImageCell;

	itStart=this->currentImageList.cellsList.begin();
	double imageCoords[dim];

	while(!currentImagePointsList.pointsList->empty())//(itPoint!=itLastPoint)
	{

		posX=(*itPoint).PointNum;
		// cout<< " posX="<<posX<<endl;
		grid->numToIntAndDoubleCoords(posX,intPointCoords,doublePointCoords);
		// printVector(doublePointCoords, dim);
		rho=dynsys->calculRho_local(doublePointCoords);
		//  cout<< " rho= "<<rho;

		/*!
		 * On calcule l'image discrète  du point
		 */

		this->computeDiscreteImageOfPoint(posX);
		/*!
		 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
		 *  de cellules
		 */
		unsigned long long int numCell,numControl;


		//cout<< "  nb cells dans l'image "<<pointDI.nbImageCells<<endl;


		//cout<<  " analyse d'une image de point\n";
		for(int i=0;i<pointDI.nbImageCells;i++)
		{
			numCell=pointDI.tabImageCells[i];
			tempImageCell.cellNum=numCell;
			//cout<< " i= "<<i<<  "numCell= "<<numCell<<endl;

			if(numCell< nbCellsTotal)
			{
				tempImageCell.minVal=PLUS_INF;
				for(int j=pointDI.tabCellEntrees[i];j<pointDI.tabCellEntrees[i+1];j++)
				{
					numControl=pointDI.tabImageControls[j];
					tempL=dynsys->lFunc(doublePointCoords, controlCoords[numControl]);
					(dynsys->*(dynsys->discretDynamics))(doublePointCoords, controlCoords[numControl], imageCoords, 1.0);

					double tempL1=dynsys->lFunc(imageCoords, controlCoords[numControl]);
					// cout<< " tempL = "<<tempL<< " tempL 1 = "<<tempL1<<endl;
					tempImageCell.minVal=min(tempImageCell.minVal, (*itPoint).minVal+rho*0.5*(tempL+tempL1));
				}
				// cout<< " ajout  de cellule avec valeur "<< tempImageCell.minVal<<"\n ";
				addDataToCurrentImage(&itStart,tempImageCell, &itNew);
				itStart=itNew;
				// cout<< " ares ajout  de l'image d'un point la taille de la liste est "<<currentImageList.cellsList.size()<<endl;
			}
		}

		//        cout<< " \n";

		itPoint++;
		currentImagePointsList.pointsList->pop_front();
	}

	cout<< " parcours de base terminé\n";


	gettimeofday(&tim,NULL);
	t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	elapsed_time=(double)((t2-t1));

	cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;



}




void ViabiMicroMacro::computeCurrIm_DD( int iter)
{

	cout<< " debut compute current image DD"<<endl;
	cout<< " currentpoints list size is "<<currentImagePointsList.pointsList->size()<<endl;
	double t1,t2,elapsed_time;

	timeval tim;
	gettimeofday(&tim,NULL);
	t1=tim.tv_sec+(tim.tv_usec/1000000.0);


	unsigned long long int posX;

	unsigned long long int nbPointsTotal=grid->getNbTotalPoints();

	unsigned long long int ** controlCoords=dynsys->getControlIntCoords();


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
	tempPointsList.maxNum =0;
	tempPointsList.minNum =0;
	if(this->whichPointListToUse ==1)
	{
		tempPointsList.pointsList= tempPointsList2;
		tempPointsList.pointsList->clear();
	}
	else
	{
		tempPointsList.pointsList= tempPointsList1;
		tempPointsList.pointsList->clear();
	}



	double tempL, tempMu;

	list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList->begin(),
			itTemp;
	imagePoint tempImagePoint;

	list<imagePoint>::iterator itStart, itNew;
	itStart=tempPointsList.pointsList->begin();
	double imageCoords[dim];

	while(!currentImagePointsList.pointsList->empty())//(itPoint!=itLastPoint)
	{

		posX=(*itPoint).PointNum;
		//cout<< " posX="<<posX<<endl;
		grid->numToIntAndDoubleCoords(posX,intPointCoords,doublePointCoords);
		//  cout<< " rho= "<<rho;

		/*!
		 * On calcule l'image discrète  du point
		 */

		this->computeDiscreteImageOfPoint_DD(posX);
		/*!
		 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
		 *  de cellules
		 */
		unsigned long long int numPoint,numControl;


		//cout<< "  nb cells dans l'image "<<pointDI_DD.nbImagePoints<<endl;


		// cout<<  " analyse d'une image de point\n";
		for(int i=0;i<pointDI_DD.nbImagePoints;i++)
		{
			numPoint=pointDI_DD.tabImagePoints[i];
			tempImagePoint.PointNum = numPoint;
			// cout<< " i= "<<i<<  "numPoint= "<<numPoint<<endl;
			//cout<< " intervalle "<<pointDI_DD.tabPointsEntrees[i]<< " "<< pointDI_DD.tabPointsEntrees[i+1]<<endl;
			if(numPoint< nbPointsTotal)
			{
				tempImagePoint.minVal=PLUS_INF;
				for(int j=pointDI_DD.tabPointsEntrees[i];j<pointDI_DD.tabPointsEntrees[i+1];j++)
				{
					numControl=pointDI_DD.tabImageControls[j];
					//cout<< " numControl = "<<numControl<<endl;
					tempL=dynsys->lFunc_fd(intPointCoords, controlCoords[numControl]);
					tempMu=dynsys->muFunc_fd(intPointCoords, controlCoords[numControl]);
					// cout<< " tempL = "<<tempL<< " tempmu = "<<tempMu<<endl;
					tempImagePoint.minVal=min(tempImagePoint.minVal, max(tempMu, (*itPoint).minVal+tempL));
				}
				// cout<< " ajout  de cellule avec valeur "<< tempImageCell.minVal<<"\n ";
				addDataToGivenPointsList(&tempPointsList,&itStart,tempImagePoint, &itNew);
				itStart=itNew;
				// cout<< " ares ajout  de l'image d'un point la taille de la liste est "<<currentImageList.cellsList.size()<<endl;
			}
		}

		//        cout<< " \n";

		itPoint++;
		currentImagePointsList.pointsList->pop_front();
	}


	cout<< " parcours de base terminé\n";
	currentImagePointsList.pointsList = tempPointsList.pointsList;

	if(this->whichPointListToUse == 1)
	{
		whichPointListToUse = 2;
	}
	else
	{
		whichPointListToUse = 1;
	}

	gettimeofday(&tim,NULL);
	t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	elapsed_time=(double)((t2-t1));

	cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;


}



/*


  void ViabiMicroMacro::computeCurrIm_Lmin( int iter)
  {


    cout<< " toto calcul de l'image de Cn  LMIN local opti nb points Cn = "<<currentImagePointsList.pointsList->size()<<endl;

    double t1,t2,elapsed_time;
    double cVect[dimC];
    for(int i=0;i<dimC;i++)
      {
      cVect[i]=0.0;
      }
    timeval tim;
    gettimeofday(&tim,NULL);
    t1=tim.tv_sec+(tim.tv_usec/1000000.0);


    int posX;

    unsigned long long int nbCellsTotal=grid->getNbTotalPoints();
    unsigned long long int numControl, numCell;
    double ** controlCoords=dynsys->getControlCoords();



    list<imageCell>::iterator itStart, itNew;

    !
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


    currentImageList.cellsList.clear();
    currentImageList.maxNum=-1;
    // currentImageList.minNum=PLUS_INF;
    currentImageList.minNum=1000+grid->nbTotalCells;
    double rho, tempL;

    list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList->begin(),
        itLastPoint=currentImagePointsList.pointsList->end(),
        itTemp;
    imageCell tempImageCell;

    itStart=this->currentImageList.cellsList.begin();


    while(!currentImagePointsList.pointsList->empty())//(itPoint!=itLastPoint)
      {

      posX=(*itPoint).PointNum;
      // cout<< " posX="<<posX<<endl;
      grid->numToIntAndDoubleCoords(posX,intPointCoords,doublePointCoords);
      // printVector(doublePointCoords, dim);
      rho=dynsys->calculRho_local(doublePointCoords);
      //  cout<< " rho= "<<rho;

      !
 * On calcule l'image discrète  du point


      this->computeDiscreteImageOfPoint(posX);
      !
 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
 *  de cellules

      unsigned long long int numCell,numControl;


      //cout<< "  nb cells dans l'image "<<pointDI.nbImageCells<<endl;


      //cout<<  " analyse d'une image de point\n";
      for(int i=0;i<pointDI.nbImageCells;i++)
        {
        numCell=pointDI.tabImageCells[i];
        tempImageCell.cellNum=numCell;
        //cout<< " i= "<<i<<  "numCell= "<<numCell<<endl;

        if(numCell< nbCellsTotal)
          {
          tempImageCell.minVal=PLUS_INF;
          for(int j=pointDI.tabCellEntrees[i];j<pointDI.tabCellEntrees[i+1];j++)
            {
            numControl=pointDI.tabImageControls[j];
            tempL=dynsys->lFunc(doublePointCoords, controlCoords[numControl]);
            tempImageCell.minVal=min(tempImageCell.minVal, (*itPoint).minVal+rho*tempL);
            }
          //cout<< " ajout  de cellule \n ";
          addDataToCurrentImage(&itStart,tempImageCell, &itNew);
          itStart=itNew;
          // cout<< " ares ajout  de l'image d'un point la taille de la liste est "<<currentImageList.cellsList.size()<<endl;
          }
        }

      //        cout<< " \n";

      itPoint++;
      currentImagePointsList.pointsList->pop_front();
      }

    cout<< " parcours de base terminé\n";


    gettimeofday(&tim,NULL);
    t2=tim.tv_sec+(tim.tv_usec/1000000.0);
    elapsed_time=(double)((t2-t1));

    cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;



  }

 */

/*

void ViabiMicroMacro::initCnIndices(  )
{
  cout<< "init Cn indices   nombre de points dans Cn = "<<currentImagePointsList.pointsList->size()<<endl;


  list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList->begin(),
      itLastPoint=currentImagePointsList.pointsList->end();

  this->indicesCn.clear();
  while(itPoint!=itLastPoint)
    {
    this->indicesCn.puch_back((*itPoint).PointNum);
    }
}
 */



void ViabiMicroMacro::computeConvexifiedImage_Lmin_omp( int iter)
{
	cout<< "new:  calcul de l'image convexifiee  de Cn , rho local, opti,  nombre de points dans Cn = "<<currentImagePointsList.pointsList->size()<<endl;
	int posX;
	list<imageCell>::iterator itStart, itNew;

	currentImageList.cellsList.clear();
	currentImageList.maxNum=-1;
	// currentImageList.minNum=PLUS_INF;
	currentImageList.minNum=1000+grid->nbTotalCells;
	double rho, tempL;

	list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList->begin(),
			itLastPoint=currentImagePointsList.pointsList->end(),
			itTemp;
	imageCell tempImageCell;


	unsigned long long int nbCellsTotal=grid->getNbTotalPoints();
	double ** controlCoords=dynsys->getControlCoords();

	itStart=this->currentImageList.cellsList.begin();

	while(itPoint!=itLastPoint)
	{
		posX=(*itPoint).PointNum;
		// cout<< " posX="<<posX<<endl;
		grid->numToIntAndDoubleCoords(posX,intPointCoords,doublePointCoords);
		// printVector(doublePointCoords, dim);
		////////system("pause");
		rho=dynsys->calculRho_local(doublePointCoords);
		//   cout<< " rho= "<<rho<<endl;
		/*!
		 * On calcule l'image discrète  du point
		 */
		//   cout<< " temp image min val = "<<tempImageCell.minVal;
		this->computeDiscreteImageOfPoint(posX);
		/*!
		 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
		 *  de cellules
		 */
		unsigned long long int numCell,numControl;


		//cout<< "  nb cells dans l'image "<<pointDI.nbImageCells<<endl;


		//cout<<  " analyse d'une image de point\n";
		for(int i=0;i<pointDI.nbImageCells;i++)
		{
			numCell=pointDI.tabImageCells[i];
			tempImageCell.cellNum=numCell;
			//cout<< " i= "<<i<<  "numCell= "<<numCell<<endl;

			if(numCell< nbCellsTotal)
			{
				tempImageCell.minVal=PLUS_INF;
				for(int j=pointDI.tabCellEntrees[i];j<pointDI.tabCellEntrees[i+1];j++)
				{
					numControl=pointDI.tabImageControls[j];
					tempL=dynsys->lFunc(doublePointCoords, controlCoords[numControl]);
					tempImageCell.minVal=min(tempImageCell.minVal, (*itPoint).minVal+rho*tempL);
				}
				//  cout<< " ajout  de cellule   "<<tempImageCell.cellNum << " valeur "<<tempImageCell.minVal<<endl;;
				addDataToCurrentImage(&itStart,tempImageCell, &itNew);
				itStart=itNew;
				addConvexCombinations(itPoint, numCell, &tempImageCell, rho,&itStart);
				// cout<< " ares ajout  de l'image d'un point la taille de la liste est "<<currentImageList.cellsList.size()<<endl;
			}
		}

		//	cout<< " le point fini \n";
		//this->showCurrentImageList();
		//	//////system("pause");
		itPoint++;
	}
	cout<< " parcours de base terminé\n";
}

void ViabiMicroMacro::computeConvexifiedImage_tmin( int iter)
{


	cout<< "new:  calcul de l'image convexifiee  de Cn , rho local, opti,  nombre de points dans Cn = "<<currentImagePointsList.pointsList->size()<<endl;



	int posX;


	double T=dynsys->getTimeHorizon();


	list<imageCell>::iterator itStart, itNew;

	currentImageList.cellsList.clear();
	currentImageList.maxNum=-1;
	// currentImageList.minNum=PLUS_INF;
	currentImageList.minNum=1000+grid->nbTotalCells;
	double rho;

	list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList->begin(),
			itLastPoint=currentImagePointsList.pointsList->end(),
			itTemp;
	imageCell tempImageCell;


	itStart=this->currentImageList.cellsList.begin();

	while(itPoint!=itLastPoint)
	{
		posX=(*itPoint).PointNum;
		// cout<< " posX="<<posX<<endl;
		grid->numToIntAndDoubleCoords(posX,intPointCoords,doublePointCoords);
		// printVector(doublePointCoords, dim);
		////////system("pause");
		rho=dynsys->calculRho_local(doublePointCoords);
		//cout<< "compute convex image  rho= "<<rho<<endl;
		if((*itPoint).minVal+rho<=T)
		{
			/*!
			 * On calcule l'image discrète  du point
			 */
			tempImageCell.minVal=(*itPoint).minVal+rho;
			//   cout<< " temp image min val = "<<tempImageCell.minVal;
			this->computeDiscreteImageOfPoint(posX);
			/*!
			 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
			 *  de cellules
			 */
			unsigned long long int numCell;

			for(int i=0;i<pointDI.nbImageCells;i++)
			{
				numCell=pointDI.tabImageCells[i];
				tempImageCell.cellNum=numCell;
				//cout<< " i= "<<i<<  "numCell= "<<numCell<<endl;

				//cout<< " ajout  de cellule \n ";
				addDataToCurrentImage(&itStart,tempImageCell, &itNew);
				itStart=itNew;
				addConvexCombinations(itPoint, numCell, &tempImageCell, rho,&itStart);
				// cout<< " ares ajout  de l'image d'un point la taille de la liste est "<<currentImageList.cellsList.size()<<endl;

				tempImageCell.minVal=(*itPoint).minVal+rho;
			}
		}
		itPoint++;

	}
}


void ViabiMicroMacro::computeConvexifiedImage_DD( int iter)
{

}

void ViabiMicroMacro::computeConvexifiedImage_Lmin( int iter)
{
	int posX;
	list<imageCell>::iterator itStart, itNew;

	currentImageList.cellsList.clear();
	currentImageList.maxNum=-1;
	// currentImageList.minNum=PLUS_INF;
	currentImageList.minNum=1000+grid->nbTotalCells;
	double rho, tempL;

	list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList->begin(),
			itLastPoint=currentImagePointsList.pointsList->end(),
			itTemp;
	imageCell tempImageCell;


	unsigned long long int nbCellsTotal=grid->getNbTotalPoints();

	double ** controlCoords=dynsys->getControlCoords();

	itStart=this->currentImageList.cellsList.begin();
	double imageCoords[dim];
	while(itPoint!=itLastPoint)
	{
		posX=(*itPoint).PointNum;
		cout<< " posX="<<posX<<endl;
		grid->numToIntAndDoubleCoords(posX,intPointCoords,doublePointCoords);
		printVector(doublePointCoords, dim);
		////////system("pause");
		rho=dynsys->calculRho_local(doublePointCoords);
		cout<< " rho= "<<rho<<endl;
		cout<<  " analyse d'une image de point ";
		printVector(doublePointCoords, dim);
		/*!
		 * On calcule l'image discrète  du point
		 */
		//   cout<< " temp image min val = "<<tempImageCell.minVal;



		cout<< "juste avant DI dans convex image posX="<<posX<<endl;
		printVector(doublePointCoords, dim);
		this->computeDiscreteImageOfPoint(posX);
		cout<<  " juste apres discrete image  ";
		printVector(doublePointCoords, dim);
		/*!
		 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
		 *  de cellules
		 */
		unsigned long long int numCell,numControl;


		//cout<< "  nb cells dans l'image "<<pointDI.nbImageCells<<endl;

		//		cout<< " =======================================================\n";
		//		cout<<  " analyse d'une image de point ";
		//		printVector(doublePointCoords, dim);
		//		cout<< "  point.minVal= "<<(*itPoint).minVal<<endl;
		//		cout<< " =======================================================\n";

		for(int i=0;i<pointDI.nbImageCells;i++)
		{
			numCell=pointDI.tabImageCells[i];
			tempImageCell.cellNum=numCell;
			//cout<< " i= "<<i<<  "numCell= "<<numCell<<endl;

			if(numCell< nbCellsTotal)
			{
				tempImageCell.minVal=PLUS_INF;
				for(int j=pointDI.tabCellEntrees[i];j<pointDI.tabCellEntrees[i+1];j++)
				{
					numControl=pointDI.tabImageControls[j];
					tempL=dynsys->lFunc(doublePointCoords, controlCoords[numControl]);
					(dynsys->*(dynsys->discretDynamics))(doublePointCoords, controlCoords[numControl], imageCoords, 1.0);
					double tempL1=dynsys->lFunc(imageCoords, controlCoords[numControl]);
					//					cout<<  " rho= "<< rho<<" ";
					//					cout<< " image "; printVector(imageCoords, dim);
					//					cout<< " tempL = "<<tempL<< " tempL 1 = "<<tempL1<< " terme integrale "<<rho*0.5*(tempL+tempL1)<<endl;
					tempImageCell.minVal=min(tempImageCell.minVal, (*itPoint).minVal+rho*0.5*(tempL+tempL1));
				}
				//cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
				//cout<< " ajout  de cellule   "<<tempImageCell.cellNum << " valeur "<<tempImageCell.minVal<<endl;;
				addDataToCurrentImage(&itStart,tempImageCell, &itNew);
				//cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
				itStart=itNew;
				addConvexCombinations(itPoint, numCell, &tempImageCell, rho,&itStart);
				// cout<< " ares ajout  de l'image d'un point la taille de la liste est "<<currentImageList.cellsList.size()<<endl;
			}
		}

		//  cout<< " le point fini \n";
		//this->showCurrentImageList();
		//  //////system("pause");
		itPoint++;
	}
	cout<< " parcours de base terminé\n";

}
/*

void ViabiMicroMacro::computeConvexifiedImage_tmin_omp( int iter)
{


  cout<< "new:  calcul de l'image convexifiee  de Cn , rho local, opti,  nombre de points dans Cn = "<<currentImagePointsList.pointsList->size()<<endl;



  int posX;


  double T=dynsys->getTimeHorizon();


  list<imageCell>::iterator itStart, itNew;
  double cVect[dimC];
  for(int i=0;i<dimC;i++)
    {
    cVect[i]=0.0;
    }

  double rho;

  vector<unsigned long long int>::iterator itPoint= this->indicesCn.begin(),
      itLastPoint=this->indicesCn.end(),
      itTemp;
  imageCell tempImageCell;


  double currentVal;
  while(itPoint!=itLastPoint)
    {
    posX=(*itPoint);
    // cout<< " posX="<<posX<<endl;
    grid->numToIntAndDoubleCoords(posX,intPointCoords,doublePointCoords);
    // printVector(doublePointCoords, dim);
    ////////system("pause");
    rho=dynsys->calculRho_local(doublePointCoords);
    cout<< "compute convex image  rho= "<<rho<<endl;
    currentVal=vTab[posX];
    if( currentVal+rho<=T)
      {
      !
 * On calcule l'image discrète  du point

      //   cout<< " temp image min val = "<<tempImageCell.minVal;
      this->computeDiscreteImageOfPoint(posX);
      !
 * L'image calculée est stockée dans la structure pointDI sous forme de tableau ordonné
 *  de cellules

      unsigned long long int numCell,numControl;


      //cout<< "  nb cells dans l'image "<<pointDI.nbImageCells<<endl;


      //cout<<  " analyse d'une image de point\n";
      for(int i=0;i<pointDI.nbImageCells;i++)
        {
        numCell=pointDI.tabImageCells[i];

        for(int ii=0;ii<nbPointsCube;ii++)
          {
          posTemp= cellNum+indicesDecalCell[ii];
          vTab[posTemp]=min( vTab[posTemp], currentVal  );

          }

        }
      }
    //  cout<< " le point fini \n";
    //this->showCurrentImageList();
    //  //////system("pause");
    itPoint++;

    }

  cout<< " parcours de base terminé\n";



}


 */





void ViabiMicroMacro::addConvexCombinations(list<imagePoint>::iterator itPoint, unsigned long long int numCell, imageCell * tempImageCell, double rho ,list<imageCell>::iterator *itStart)
{

	/*!
	 * On appelle d'abord deux fois la fonction numToIntAndDoubleCoords()
	 * pour calculer les coordonnées réelles du point de départ x  et  du coin inférieur de la maille
	 *  appartenant é l'image \f$ \Phi(x)\f$ y.
	 */


	unsigned long long int posX=(*itPoint).PointNum;
	double pointVal=(*itPoint).minVal;
	double newCellVal=(*tempImageCell).minVal;

	double LVal=(newCellVal-pointVal)/rho;


	//cout<< "add convex comnbin  posX= "<<posX<< " rho ="<<rho<<endl;
	grid->numToIntAndDoubleCoords(posX,intPointCoords,doublePointCoords);
	//cout<< " point de départ ";
	//	printVector(doublePointCoords, dim);
	//	cout<< " numCell = "<<numCell<<endl;
	grid->numToIntAndDoubleCoords(numCell,intPointCoords,doubleVect);
	//cout<< " point arrivee  ";
	//printVector(doubleVect, dim);
	list<imageCell>::iterator  itNew;
	double dist=0.;
	/*!
	 * Ensuite on calcule le vecteur différence \f$ z=y-x\f$ et sa norme \f$ \|z\|_2 \f$ .
	 */
	for(int i=0;i<dim;i++)
	{
		doubleVect[i]=doubleVect[i]-doublePointCoords[i];
		dist+=doubleVect[i]*doubleVect[i];
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
	double deltat=min(0.1, grid->getMaxStep()/(2.0*sqrt(dist)));
	//  cout<< " valeu depart "<<pointVal<< " valeur arrivee "<<newCellVal<<  " deltat= "<<deltat<< " Lval= "<<LVal<<endl;
	//	cout<< "   norme de difference "<<sqrt(dist)<< " pas maxi "<<grid->getMaxStep()<<endl<< " difference vecteur ";
	//printVector(doubleVect, dim);
	double t=deltat;
	unsigned long long int newCellNum, lastVisitCellNum=grid->getNbTotalCells()+1;
	while (t<1.0)
	{

		for(int i=0;i<dim;i++)
		{
			doubleVect1[i]=doublePointCoords[i]+t*doubleVect[i];
		}

		// cout<< "  point intermediaire  ";
		// printVector(doubleVect1, dim);
		if(grid->isPointInGrid(doubleVect1) )
		{
			//cout<< " le point est das la grlle\n";
			/*!
			 * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
			 */

			if( dynsys->constraintsX(doubleVect1)<PLUS_INF)
			{
				// cout<< "  contraintes sur   X  ok\n";

				/*!
				 * Si l'image  est dans l'ensemble de contraintes sur l'état \f$ K \f$
				 * on calcule le numéro de maille qui contient cette image
				 */
				newCellNum=grid->localizePoint(doubleVect1); // on enregistre le numero de maille
				if(newCellNum!=lastVisitCellNum)
				{
					//cout<< " maille visee est "<<newCellNum<<endl;
					(*tempImageCell).cellNum=newCellNum;
					(*tempImageCell).minVal=  pointVal+t*rho*LVal;
					//itStart=this->currentImageList.cellsList.begin();
					//   cout<<  "  ajout de celule "<<newCellNum<< " avec val = "<<t*rho*LVal<<endl;
					this->addDataToCurrentImage(itStart,(*tempImageCell), &itNew);
					lastVisitCellNum=newCellNum;
					(*itStart)=itNew;
					//this->showCurrentImageList();
					////////system("pause");
				}
			}
		}
		t=t+deltat;

	}
	////////system("pause");

}
void ViabiMicroMacro::addDataToPointsList(list<imagePoint>::iterator *startIt, imagePoint newPoint,list<imagePoint>::iterator *resIt )
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
	// cout<< " l'ajout du point itStart pointe sur "<<(*(*startIt)).PointNum <<endl;

	unsigned long long int numnewPoint=newPoint.PointNum;
	// cout<< " l'ajout du point numPOint "<<numnewPoint << " maxNum = "<<currentImagePointsList.maxNum<<endl;
	list<imagePoint>::iterator itCell;

	if(currentImagePointsList.pointsList->size()==0)
	{
		currentImagePointsList.pointsList->push_back(newPoint);
		currentImagePointsList.maxNum=numnewPoint;
		currentImagePointsList.minNum=numnewPoint;
		(*resIt)=currentImagePointsList.pointsList->end();
		(*resIt)--;
		// cout<< " on ajoute le premier element nouveax min et max  de la liste sont "<<currentImagePointsList.minNum<< " "<<currentImagePointsList.maxNum<<endl;
	}
	else
	{
		if(numnewPoint>currentImagePointsList.maxNum)
		{
			currentImagePointsList.pointsList->push_back(newPoint);
			currentImagePointsList.maxNum=numnewPoint;
			(*resIt)=currentImagePointsList.pointsList->end();
			(*resIt)--;
			// cout<< " on ajoute é la fin  nouveax min et max  de la liste sont "<<currentImagePointsList.minNum<< " "<<currentImagePointsList.maxNum<<endl;
		}
		else
		{
			if(numnewPoint<currentImagePointsList.minNum)
			{
				currentImagePointsList.pointsList->push_front(newPoint);
				currentImagePointsList.minNum=numnewPoint;
				(*resIt)=currentImagePointsList.pointsList->begin();
				// cout<< " on ajoute au debut nouveax min et max  de la liste sont "<<currentImagePointsList.minNum<< " "<<currentImagePointsList.maxNum<<endl;

			}
			else
			{
				itCell=*startIt;

				if(numnewPoint<(*itCell).PointNum)
				{

					while( (numnewPoint<(*itCell).PointNum))
					{
						//cout<< " current cell num "<<(*itCell).PointNum<< " new cell num "<<numnewPoint<<endl;
						itCell--;
					}
					if(numnewPoint>(*itCell).PointNum)
					{
						itCell++;
						currentImagePointsList.pointsList->insert(itCell,newPoint);
					}
					else
					{
						(this->*addDataToCurrentPoint)(itCell, newPoint);
					}
					(*resIt)=itCell;

				}
				else
				{	while( (numnewPoint>(*itCell).PointNum))
				{
					//	cout<< " current point num "<<(*itCell).PointNum<< " new point  num "<<numnewPoint<<endl;
					itCell++;
				}
				if(numnewPoint<(*itCell).PointNum)
				{
					currentImagePointsList.pointsList->insert(itCell,newPoint);
				}
				else
				{
					(this->*addDataToCurrentPoint)(itCell, newPoint);
				}
				(*resIt)=itCell;

				}
			}
		}
	}

	// cout<< " l'ajout du point est terminé  on a ajouté juste avant "<<(*(*resIt)).PointNum<< " size " << currentImagePointsList.pointsList->size()<< endl;
	////////system("pause");
}

void ViabiMicroMacro::addDataToGivenPointsList(imagePointsList * tempImagePointsList, list<imagePoint>::iterator *startIt, imagePoint newPoint,list<imagePoint>::iterator *resIt )
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
	//cout<< " l'ajout du point itStart pointe sur "<<(*(*startIt)).PointNum <<endl;

	unsigned long long int numnewPoint=newPoint.PointNum;
	//cout<< " l'ajout du point numPOint "<<numnewPoint << " maxNum = "<<currentImagePointsList.maxNum<<endl;
	list<imagePoint>::iterator itCell;

	if(tempImagePointsList->pointsList->size()==0)
	{
		tempImagePointsList->pointsList->push_back(newPoint);
		tempImagePointsList->maxNum=numnewPoint;
		tempImagePointsList->minNum=numnewPoint;
		(*resIt)=tempImagePointsList->pointsList->end();
		(*resIt)--;
		//cout<< " on ajoute le premier element nouveax min et max  de la liste sont "<<tempImagePointsList->minNum<< " "<<tempImagePointsList->maxNum<<endl;
	}
	else
	{
		if(numnewPoint>tempImagePointsList->maxNum)
		{
			tempImagePointsList->pointsList->push_back(newPoint);
			tempImagePointsList->maxNum=numnewPoint;
			(*resIt)=tempImagePointsList->pointsList->end();
			(*resIt)--;
			//cout<< " on ajoute é la fin  nouveax min et max  de la liste sont "<<tempImagePointsList->minNum<< " "<<tempImagePointsList->maxNum<<endl;
		}
		else
		{
			if(numnewPoint<tempImagePointsList->minNum)
			{
				tempImagePointsList->pointsList->push_front(newPoint);
				tempImagePointsList->minNum=numnewPoint;
				(*resIt)=tempImagePointsList->pointsList->begin();
				//cout<< " on ajoute au debut nouveax min et max  de la liste sont "<<tempImagePointsList->minNum<< " "<<tempImagePointsList->maxNum<<endl;

			}
			else
			{
				itCell=*startIt;

				if(numnewPoint<(*itCell).PointNum)
				{

					while( (numnewPoint<(*itCell).PointNum))
					{
						//cout<< " current cell num "<<(*itCell).PointNum<< " new cell num "<<numnewPoint<<endl;
						itCell--;
					}
					if(numnewPoint>(*itCell).PointNum)
					{
						itCell++;
						tempImagePointsList->pointsList->insert(itCell,newPoint);
					}
					else
					{
						(this->*addDataToCurrentPoint)(itCell, newPoint);
					}
					(*resIt)=itCell;

				}
				else
				{	while( (numnewPoint>(*itCell).PointNum))
				{
					//	cout<< " current point num "<<(*itCell).PointNum<< " new point  num "<<numnewPoint<<endl;
					itCell++;
				}
				if(numnewPoint<(*itCell).PointNum)
				{
					tempImagePointsList->pointsList->insert(itCell,newPoint);
				}
				else
				{
					(this->*addDataToCurrentPoint)(itCell, newPoint);
				}
				(*resIt)=itCell;

				}
			}
		}
	}

	//cout<< " l'ajout du point est terminé  on a ajouté juste avant "<<(*(*resIt)).PointNum<<endl;
	////////system("pause");
}


void ViabiMicroMacro::addDataToPoint(list<imagePoint>::iterator itCell, imagePoint newPoint)
{
	(*itCell).minVal=min( (*itCell).minVal, newPoint.minVal);
}


void ViabiMicroMacro::createPointsList_DD()
{

}
void ViabiMicroMacro::createPointsList()
{

	cout<< " create point list opti new\n";
	list<triple>::iterator itR;
	double t1,t2,elapsed_time;

	timeval tim;
	gettimeofday(&tim,NULL);
	t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	unsigned long long int posX;
	double testV[dim];
	unsigned long long int testI[dim];
	currentImagePointsList.pointsList->clear();
	currentImagePointsList.maxNum=-1;
	currentImagePointsList.minNum=currentImageList.minNum;


	int  nbPointsCube=(int) pow(2.0,dim);//pow(2.0, dim);
	list<imageCell>::iterator itCell=currentImageList.cellsList.begin(),
			itLast=currentImageList.cellsList.end();
	int i;

	long long  int * indicesDecalCell=grid->getIndicesDecalCell();



	imagePoint tempPoint;
	list<imagePoint>::iterator itStart, itNew;
	itStart=this->currentImagePointsList.pointsList->begin();
	for(i=0;i<nbPointsCube-1;i++)
	{
		while(itCell!=itLast)
		{
			//numCell=(*itCell).cellNum;
			tempPoint.minVal=(*itCell).minVal;

			posX=(*itCell).cellNum+indicesDecalCell[i];
			grid->numToIntAndDoubleCoords( posX ,testI,testV);
			if(dynsys->constraintsX(testV)<PLUS_INF)
			{

				tempPoint.PointNum=(*itCell).cellNum+indicesDecalCell[i];

				this->addDataToPointsList(&itStart,tempPoint,&itNew);
				itStart=itNew;

			}
			itCell++;

		}
		itCell=currentImageList.cellsList.begin();
	}

	while(!currentImageList.cellsList.empty())//(itCell!=itLast)
	{
		//numCell=(*itCell).cellNum;
		tempPoint.minVal=(*itCell).minVal;

		//	cout<< "  maille numéro "<<numCell<<endl;

		//posX=numCell+indicesDecalCell[nbPointsCube-1];

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
		//cout<< " on ajoute le point numero "<<posX<<endl;
		tempPoint.PointNum=(*itCell).cellNum+indicesDecalCell[nbPointsCube-1];
		this->addDataToPointsList(&itStart,tempPoint,&itNew);
		itStart=itNew;


		itCell++;
		currentImageList.cellsList.pop_front();
	}



	gettimeofday(&tim,NULL);
	t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	elapsed_time=(double)((t2-t1));

	cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;

}


void ViabiMicroMacro::addDataToCurrentImage(list<imageCell>::iterator *startIt, imageCell newCell,list<imageCell>::iterator *resIt )
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

	unsigned long long int numNewCell=newCell.cellNum;

	list<imageCell>::iterator itCell, itLast=currentImageList.cellsList.end();

	if(currentImageList.cellsList.size()==0)
	{
		currentImageList.cellsList.push_back(newCell);
		currentImageList.maxNum=numNewCell;
		currentImageList.minNum=numNewCell;
		(*resIt)=currentImageList.cellsList.end();
		(*resIt)--;
	}
	else
	{


		if((int)numNewCell>currentImageList.maxNum)
		{
			currentImageList.cellsList.push_back(newCell);
			currentImageList.maxNum=numNewCell;
			currentImageList.minNum=min(currentImageList.maxNum,currentImageList.minNum);
			(*resIt)=currentImageList.cellsList.end();
			(*resIt)--;
		}
		else
		{
			if((int)numNewCell<currentImageList.minNum)
			{
				currentImageList.cellsList.push_front(newCell);
				currentImageList.minNum=numNewCell;
				currentImageList.maxNum=max(currentImageList.maxNum,currentImageList.minNum);
				(*resIt)=currentImageList.cellsList.begin();
			}
			else
			{
				itCell=*startIt;


				if((numNewCell>(*itCell).cellNum))
				{
					while((itCell!=itLast) && (numNewCell>(*itCell).cellNum))
					{
						//cout<< " current cell num "<<(*itCell).cellNum<< " new cell num "<<numNewCell<<endl;
						itCell++;
					}
					if(numNewCell<(*itCell).cellNum)
					{
						currentImageList.cellsList.insert(itCell,newCell);
					}
					else
					{
						(this->*addDataToCurrentCell)(itCell, newCell);
					}
					(*resIt)=itCell;
				}
				else
				{
					while(  (numNewCell<(*itCell).cellNum))
					{
						//cout<< " current cell num "<<(*itCell).cellNum<< " new cell num "<<numNewCell<<endl;
						itCell--;
					}
					if(numNewCell>(*itCell).cellNum)
					{
						itCell++;
						currentImageList.cellsList.insert(itCell,newCell);
					}
					else
					{
						(this->*addDataToCurrentCell)(itCell, newCell);
					}
					(*resIt)=itCell;
				}

			}
		}

	}
}


void ViabiMicroMacro::addDataToCell(list<imageCell>::iterator itCell, imageCell newCell)
{

	(*itCell).minVal=min((*itCell).minVal,newCell.minVal);

}



/*!
 * Le destructeur nettoie la mémoire réservée pour certaines variables globales é la classe et ferme
 * les bases de données de rétro-actions
 */
ViabiMicroMacro::~ViabiMicroMacro() {

	cout<< " desturction de classe viabi hjb \n";
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


	cout<< " divers tableaux OK\n";


}

void ViabiMicroMacro::printViabiInfo()
{
	grid->printGrid();
}



void ViabiMicroMacro::showCurrentImageList()
{
	imageCell c;
	list<imageCell>::iterator itCell=currentImageList.cellsList.begin(), itLast=currentImageList.cellsList.end();
	list<triple>::iterator it;
	cout<< " image  de CN calculée  est la suivante \n";
	cout<< "*******************************************************\n";
	while((itCell!=itLast))
	{
		c=(*itCell);

		cout<< " maille  num "<<c.cellNum;
		cout<< " value= "<< c.minVal<<endl;


		itCell++;
	}
	cout<< "*******************************************************\n";
}
void ViabiMicroMacro::showCurrentImagePointsList()
{
	imagePoint c;
	list<imagePoint>::iterator itCell=currentImagePointsList.pointsList->begin(),
			itLast=currentImagePointsList.pointsList->end();
	cout<< " image  de CN calculée  est la suivante : liste de POINTS \n";
	cout<< "*******************************************************\n";
	while((itCell!=itLast))
	{
		c=(*itCell);

		cout<< " point  num "<<c.PointNum;
		cout<< " value= "<< c.minVal;
		itCell++;
	}
	cout<< "*******************************************************\n";
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

	cout<< " ajout de nouveaux points\n";

	double t1,t2,elapsed_time;

	timeval tim;
	gettimeofday(&tim,NULL);
	t1=tim.tv_sec+(tim.tv_usec/1000000.0);

	int nbNewPoints=0;

	list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList->begin(),
			itLastPoint=currentImagePointsList.pointsList->end(),
			itTemp;

	imagePoint tempPoint;
	while(itPoint!=itLastPoint)
	{

		tempPoint=(*itPoint);
		if((*itPoint).minVal<vTab[(*itPoint).PointNum])
		{
			nbNewPoints++;
			vTab[(*itPoint).PointNum]=min(vTab[(*itPoint).PointNum], (*itPoint).minVal);
			itPoint++;
		}
		else
		{
			itTemp=itPoint;
			itPoint++;
			currentImagePointsList.pointsList->erase(itTemp);
		}
	}

	gettimeofday(&tim,NULL);
	t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	elapsed_time=(double)((t2-t1));

	cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;
	return nbNewPoints;
}


bool ViabiMicroMacro::testConstraintesForCell(unsigned long long int numCell)
{
	bool res=true;
	long long  int * indicesDecalCell=grid->getIndicesDecalCell();
	double doublePointCoords[dim];


	int  nbPointsCube=(int) pow(2.0,dim);//pow(2.0, dim);
	unsigned long long int posX;
	int i=0;

	//cout<< "  maille numéro "<<numCell<<endl;
	while( res & (i<nbPointsCube))
	{
		posX=numCell+indicesDecalCell[i];
		grid->numToIntAndDoubleCoords(posX, intPointCoords, doublePointCoords);
		res=(dynsys->constraintsX(doublePointCoords)<PLUS_INF);

		i++;

	}



	return res;


}

void ViabiMicroMacro::initialiseConstraints()
{
	cout<< "  initialisation de l'ensemble K0 pour le calcul du noyau\n";

	if(!(dynsys->getDynType()==DD)   )
	{
		initialiseConstraints_CC();
	}
	else
	{
		initialiseConstraints_DD();
	}
}

void ViabiMicroMacro::initialiseConstraints_CC()
{
	cout<< "  initialisation de l'ensemble K0 pour le calcul du noyau\n";

	unsigned long long int iCoords[dim];
	double xCoords[dim];
	unsigned long long int totalPointsX=grid->getNbTotalPoints();
	for(unsigned long long int pos=0;pos<totalPointsX;pos++)
	{
		grid->numToIntAndDoubleCoords(pos, iCoords,xCoords);
		// cout<< "  coords ";
		//printVector(xCoords, dim);
		vTab[pos]=dynsys->constraintsX(xCoords);
	}
}


void ViabiMicroMacro::initialiseConstraints_DD()
{
	cout<< "  initialisation de l'ensemble K0 pour le calcul du noyau\n";

	unsigned long long int iCoords[dim];
	double xCoords[dim];
	unsigned long long int totalPointsX=grid->getNbTotalPoints();
	for(unsigned long long int pos=0;pos<totalPointsX;pos++)
	{
		grid->numToIntAndDoubleCoords(pos, iCoords,xCoords);
		// cout<< "  coords ";
		//printVector(xCoords, dim);
		vTab[pos]=dynsys->constraintsX_fd(iCoords);
	}
}



void ViabiMicroMacro::ViabilityKernel( bool sortieOK,int nbArret)
{
	if(!(dynsys->getDynType()==DD)   )
	{
		viabKerValFunc();
	}
	else
	{
		viabKerValFunc_DD();
	}
}

void ViabiMicroMacro::GarantedViabilityKernel( bool sortieOK,int nbArret)
{

}

