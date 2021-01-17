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


ViabiMicroMacro::ViabiMicroMacro(GridMicroMacro * gr, SysDyn* sd, algoViabiParams avp) {
  // TODO Auto-generated constructor stub
  /*!
   * l'initialisation commence par attribution de valeurs
   * au pointeurs de grille et de systeme
   */
  grid=gr;
  dynsys = sd;

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
  pointDI.tabCellEntrees=new int[sd->getTotalNbPointsC()+1];
  pointDI.tabImageCells=new int[sd->getTotalNbPointsC()];
  pointDI.tabImageControls=new int[sd->getTotalNbPointsC()];

  /*
   * Initialisation de tableaux servant de variables globales pour certaines fonctions
   * Cela évte de multiples allocations/destructions de mémoire
   */
  intPointCoords=new unsigned long long int[dim];
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
  ViabiMicroMacro::createCurrentPointsList=&ViabiMicroMacro::createPointsList;


  if(avp.COMPUTE_TMIN)
    {
    ViabiMicroMacro::computeFirstConvexifiedImage=&ViabiMicroMacro::computeConvexifiedImage_tmin;
    ViabiMicroMacro::computeCurrentImage=&ViabiMicroMacro::computeCurrIm_tmin;
    //  ViabiMicroMacro::computeFirstConvexifiedImage_omp=&ViabiMicroMacro::computeConvexifiedImage_tmin_omp;
    }
  else
    {
    //  ViabiMicroMacro::computeFirstConvexifiedImage_omp=&ViabiMicroMacro::computeConvexifiedImage_Lmin_omp;
    ViabiMicroMacro::computeCurrentImage=&ViabiMicroMacro::computeCurrIm_Lmin;
    ViabiMicroMacro::computeFirstConvexifiedImage=&ViabiMicroMacro::computeConvexifiedImage_Lmin;
    }

}

void ViabiMicroMacro::computeDiscreteImageOfPoint(unsigned long long int num)
{


  // cout<< " calcul de l'image d'un point \n";
  double ** controlCoords=dynsys->getControlCoords();
  unsigned long long int nbCTotal=dynsys->getTotalNbPointsC();
  unsigned long long int nbCellsTotal=grid->getNbTotalPoints();




  grid->numToIntAndDoubleCoords(num,intPointCoords,doublePointCoords);
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

      (dynsys->*(dynsys->discretDynamics))(doublePointCoords, controlCoords[cu], doubleVect1, 1.0);

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
        imageCells[cu]=nbCellsTotal+1; // sinon on enregistre un nombre convenu reconnaissanble
        cellsList.push_back(intPair(imageCells[cu],cu));
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
  while((itCell!=cellsList.end()) &( (*itCell).first<(int)nbCellsTotal))
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

void ViabiMicroMacro::initialiseTarget() const
{
  /*!
   *  cette fonction initialise la base de données pour la dynamique. Elle
   *  ajoute dnas la base les premiers points pour lesquels la fonction valeur
   *  est réelle. Au début de l'algorithme de bassin de capture
   *    seuls les points de la  cible  ont une fonction valeur réelle
   *
   */
  unsigned long long int dim=grid->dim;

  unsigned long long int  * x=new unsigned long long int [dim];

  double * xReel =new double[dim];

  int totalPointsX=grid->getNbTotalPoints();
  unsigned long long int pos ;


  /*
   *  on parcourt  tous les points de l'espace discret  fini
   *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
   */
  for( pos=0;pos<(unsigned long long int)totalPointsX;pos++)
    {
    grid->numToIntAndDoubleCoords(pos,x,xReel);
    vTab[pos]=(double)dynsys->target(xReel);
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


  double c;


  unsigned long long int dim=grid->dim;

  //ouverture de la base


  unsigned long long int  * x=new unsigned long long int [dim];

  double * xReel =new double[dim];



  int  totalPointsC=0;  // nombre de points  de l'espace des commandes




  int totalPointsX=grid->getNbTotalPoints();

  imagePoint currentPoint;

  currentImagePointsList.maxNum=-1;
  currentImagePointsList.minNum=0;
  currentImagePointsList.pointsList=list<imagePoint>();

  unsigned long long int pos ;

  list<imagePoint>::iterator itStart=currentImagePointsList.pointsList.begin(),
      itNew;
  int cpt=0;


  /*
   *  on parcourt  tous les points de l'espace discret  fini
   *   et  on choisit les points oé la fonction cible  renvoie une faleur finie
   */
  for( pos=0;pos<(unsigned long long int)totalPointsX;pos++)
    {

    // cout<< " pos= "<<pos<<endl;
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
      cout<<" fini ajout point  \n";
      printVector(xReel,dim);
      //////system("pause");
      addDataToPointsList(&itStart, currentPoint, &itNew);
      itStart=itNew;
      cpt++;
      }
    }
  cout<< " fini creation de target nb de points "<<cpt<<endl;

}

void ViabiMicroMacro::viabKerValFunc()
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

          for(int iCell=0;iCell<nbPointsCube-1;iCell++)
            {
            for(int j=pointDI.tabCellEntrees[compteComm];j<pointDI.tabCellEntrees[compteComm+1];j++)
              {
              numControl=pointDI.tabImageControls[j];
              tempL=dynsys->lFunc(rCoords, controlCoords[numControl]);
              tempM=dynsys->mFunc(rCoords, controlCoords[numControl]);
              tempVal=(vTab[cellNum+indicesDecalCell[iCell]]+rho*tempL)/(1-rho*tempM);
              minValCell=min(minValCell, tempVal);
              }
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
}

void ViabiMicroMacro::captBasinEpi()
{

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
  nbNewPoints=(this->*addNewPointsToSet)();
  cout<< " points ajoutes =  "<<nbNewPoints<<endl;

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
    // cout<< " points ajoutes =  "<<nbNewPoints<<endl;
    iter++;
    }
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

  for(int i=0;i<(int)dim;i++)
    {
    xCoordsDouble[i]=currentPos[i];
    }

  int maxnbViabPoints;
  int bestCu;
  rho=dt;
  dynsys->setRho(rho);
  cout<< " find control :  rho= "<<rho<<endl;
  double       currentVal=PLUS_INF;
  for(int ii=0;ii<nbPointsCube;ii++)
    {
    posTemp= cellNum+indicesDecalCell[ii];
    grid->numToIntAndDoubleCoords( posTemp ,testI,testV);
    currentVal=min( vTab[posTemp],currentVal );

    }
  cout<< " val of current point is "<<currentVal<<endl;
  /*
   * on parcours tous les contrôles
   */
  maxnbViabPoints=0;
  bestCu=0;
  double minVal=PLUS_INF, minValCell;
  int iter=0;

  testNonVide=false;
  while(iter<nbStepIter )
    {
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
            /* cout<< " image du point ";
            for(int k=0;k<dim;k++)
              {
              cout<< " "<<doubleVect1[k];
              }*/
            /*
             * le successeur vérifie les contraintes
             * On identifie la maille où il se trouve
             */
            cellNum=grid->localizePoint(doubleVect1);
            // cout<< " num cellule "<<cellNum<<endl;

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
              double dist=0.0;
              for(int k=0;k<(int)dim;k++)
                {
                dist=max(dist, abs(testV[k]-doubleVect1[k]));
                }

              minValCell=min( vTab[posTemp],minValCell );

              }

            testNonVide=(minVal<=currentVal);
            if(minValCell<minVal)
              {
              minVal=minValCell;
              bestCu=cu;
              cout<< " min val = "<<minVal<<endl;
              }
            }
          }
        }
      cu++;
      }//fin de parcours de tous les contrôles
    // la boucle s'arête ici u premier contrôle
    // qui donne un successeur viable

    iter++;
    rho=rho*stepCoeff;
    dynsys->setRho(rho);
    cout<< "find control  iteration  "<<iter<<" rho= "<<rho<<endl;
    }

  succes=testNonVide;
  nbViabVoisins=maxnbViabPoints;
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
      minCellVal=PLUS_INF;
      for(int ii=0;ii<nbPointsCube;ii++  )
        {
        posTemp= cellNum+indicesDecalCell[ii];
        minCellVal=min(minCellVal, vTab[posTemp]);
        }

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

          rho= 1.5* dynsys->calculRho_local(xCoordsDouble);

          rho=min(rho, T-time);

          cout<< " rho= "<<rho<<endl;

          bestCu=this->findOptiControl_tmin(xCoordsDouble, rho,1,1.0, imageVect,maxnbViabPoints, testNonVide );
          time+=rho;

          // la boucle s'arête ici u premier contrôle
          // qui donne un successeur viable

          cout<<   " Premiere recherche de controle viable  fini parcours de controles on a test interieur = "<<testviabInt<<
              " test non vide "<<testNonVide<< " maxnbViabPoints =  "<<maxnbViabPoints<< " bes c u= "<<bestCu<<endl;
          if(testNonVide)
            {
            // contrôle viable trouvé
            // on recopie ce contrôle dans la liste et
            // le successeur devient le point  courent
            if(testviabInt)
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
              }
            else
              {
              cout<< "  recherche de optimal viable avec iterations  sur le pas  de temps\n";
              rho=  1.5* dynsys->calculRho_local(xCoordsDouble);

              rho=min(rho, T-time);

              bestCu=this->findOptiControl_tmin(xCoordsDouble, rho,10,0.9,imageVect,maxnbViabPoints, testNonVide );

              for(int dc=0;dc<(int)dimC;dc++)
                {
                trajControlCoords[dc]=controlCoords[bestCu][dc];
                }
              trajC.push_back(trajControlCoords);
              cout<< "  coords double : ";
              for(int i=0;i<(int)dim;i++)
                {
                xCoordsDouble[i]=imageVect[i];
                cout<< " "<<xCoordsDouble[i];
                }
              cout<<endl;
              }
            }
          else
            {
            cout<<"   Echec! Sortie de l'ensemble viable \n";
            break;

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


  cout<< " toto calcul de l'image de Cn local opti nb points Cn = "<<currentImagePointsList.pointsList.size()<<endl;

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

  list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList.begin(),
      itTemp;
  imageCell tempImageCell;

  itStart=this->currentImageList.cellsList.begin();

  while(!currentImagePointsList.pointsList.empty())//(itPoint!=itLastPoint)
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
    currentImagePointsList.pointsList.pop_front();
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

  list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList.begin(),
      itTemp;
  imageCell tempImageCell;

  itStart=this->currentImageList.cellsList.begin();
  double imageCoords[dim];

  while(!currentImagePointsList.pointsList.empty())//(itPoint!=itLastPoint)
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
    currentImagePointsList.pointsList.pop_front();
    }

  cout<< " parcours de base terminé\n";


  gettimeofday(&tim,NULL);
  t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  elapsed_time=(double)((t2-t1));

  cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;



}

/*


  void ViabiMicroMacro::computeCurrIm_Lmin( int iter)
  {


    cout<< " toto calcul de l'image de Cn  LMIN local opti nb points Cn = "<<currentImagePointsList.pointsList.size()<<endl;

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

    list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList.begin(),
        itLastPoint=currentImagePointsList.pointsList.end(),
        itTemp;
    imageCell tempImageCell;

    itStart=this->currentImageList.cellsList.begin();


    while(!currentImagePointsList.pointsList.empty())//(itPoint!=itLastPoint)
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
      currentImagePointsList.pointsList.pop_front();
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
  cout<< "init Cn indices   nombre de points dans Cn = "<<currentImagePointsList.pointsList.size()<<endl;


  list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList.begin(),
      itLastPoint=currentImagePointsList.pointsList.end();

  this->indicesCn.clear();
  while(itPoint!=itLastPoint)
    {
    this->indicesCn.puch_back((*itPoint).PointNum);
    }
}
 */



void ViabiMicroMacro::computeConvexifiedImage_Lmin_omp( int iter)
{
  cout<< "new:  calcul de l'image convexifiee  de Cn , rho local, opti,  nombre de points dans Cn = "<<currentImagePointsList.pointsList.size()<<endl;
  int posX;
  list<imageCell>::iterator itStart, itNew;

  currentImageList.cellsList.clear();
  currentImageList.maxNum=-1;
  // currentImageList.minNum=PLUS_INF;
  currentImageList.minNum=1000+grid->nbTotalCells;
  double rho, tempL;

  list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList.begin(),
      itLastPoint=currentImagePointsList.pointsList.end(),
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


  cout<< "new:  calcul de l'image convexifiee  de Cn , rho local, opti,  nombre de points dans Cn = "<<currentImagePointsList.pointsList.size()<<endl;



  int posX;


  double T=dynsys->getTimeHorizon();


  list<imageCell>::iterator itStart, itNew;

  currentImageList.cellsList.clear();
  currentImageList.maxNum=-1;
  // currentImageList.minNum=PLUS_INF;
  currentImageList.minNum=1000+grid->nbTotalCells;
  double rho;

  list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList.begin(),
      itLastPoint=currentImagePointsList.pointsList.end(),
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
    cout<< "compute convex image  rho= "<<rho<<endl;
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


void ViabiMicroMacro::computeConvexifiedImage_Lmin( int iter)
{
   int posX;
  list<imageCell>::iterator itStart, itNew;

  currentImageList.cellsList.clear();
  currentImageList.maxNum=-1;
  // currentImageList.minNum=PLUS_INF;
  currentImageList.minNum=1000+grid->nbTotalCells;
  double rho, tempL;

  list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList.begin(),
      itLastPoint=currentImagePointsList.pointsList.end(),
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

    cout<< " =======================================================\n";
    cout<<  " analyse d'une image de point ";
    printVector(doublePointCoords, dim);
    cout<< "  point.minVal= "<<(*itPoint).minVal<<endl;
    cout<< " =======================================================\n";

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
          cout<<  " rho= "<< rho<<" ";
          cout<< " image "; printVector(imageCoords, dim);
          cout<< " tempL = "<<tempL<< " tempL 1 = "<<tempL1<< " terme integrale "<<rho*0.5*(tempL+tempL1)<<endl;
          tempImageCell.minVal=min(tempImageCell.minVal, (*itPoint).minVal+rho*0.5*(tempL+tempL1));
          }
        cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        cout<< " ajout  de cellule   "<<tempImageCell.cellNum << " valeur "<<tempImageCell.minVal<<endl;;
        addDataToCurrentImage(&itStart,tempImageCell, &itNew);
        cout<< "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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


  cout<< "new:  calcul de l'image convexifiee  de Cn , rho local, opti,  nombre de points dans Cn = "<<currentImagePointsList.pointsList.size()<<endl;



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
  //cout<< " l'ajout du point itStart pointe sur "<<(*(*startIt)).PointNum <<endl;

  unsigned long long int numnewPoint=newPoint.PointNum;

  list<imagePoint>::iterator itCell;
  if(numnewPoint>currentImagePointsList.maxNum)
    {
    currentImagePointsList.pointsList.push_back(newPoint);
    currentImagePointsList.maxNum=numnewPoint;
    (*resIt)=currentImagePointsList.pointsList.end();
    (*resIt)--;
    //cout<< " on ajoute é la fin  nouveax min et max  de la liste sont "<<currentImagePointsList.minNum<< " "<<currentImagePointsList.maxNum<<endl;
    }
  else
    {
    if(numnewPoint<currentImagePointsList.minNum)
      {
      currentImagePointsList.pointsList.push_front(newPoint);
      currentImagePointsList.minNum=numnewPoint;
      (*resIt)=currentImagePointsList.pointsList.begin();
      //	cout<< " on ajoute au debut nouveax min et max  de la liste sont "<<currentImagePointsList.minNum<< " "<<currentImagePointsList.maxNum<<endl;

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
          currentImagePointsList.pointsList.insert(itCell,newPoint);
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
          currentImagePointsList.pointsList.insert(itCell,newPoint);
          }
        else
          {
          (this->*addDataToCurrentPoint)(itCell, newPoint);
          }
        (*resIt)=itCell;

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
  currentImagePointsList.pointsList.clear();
  currentImagePointsList.maxNum=-1;
  currentImagePointsList.minNum=currentImageList.minNum;


  int  nbPointsCube=(int) pow(2.0,dim);//pow(2.0, dim);
  list<imageCell>::iterator itCell=currentImageList.cellsList.begin(),
      itLast=currentImageList.cellsList.end();
  int i;

  long long  int * indicesDecalCell=grid->getIndicesDecalCell();



  imagePoint tempPoint;
  list<imagePoint>::iterator itStart, itNew;
  itStart=this->currentImagePointsList.pointsList.begin();
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

void ViabiMicroMacro::printViabiInfo() const
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
  list<imagePoint>::iterator itCell=currentImagePointsList.pointsList.begin(),
      itLast=currentImagePointsList.pointsList.end();
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

  list<imagePoint>::iterator itPoint=currentImagePointsList.pointsList.begin(),
      itLastPoint=currentImagePointsList.pointsList.end(),
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
      currentImagePointsList.pointsList.erase(itTemp);
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

void ViabiMicroMacro::initializeConstrains()
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

