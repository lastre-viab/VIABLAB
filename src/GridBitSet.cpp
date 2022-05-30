/*
 * GridBitSet.cpp
 *
 *  Created on: 10 sept. 2013
 *      Author: ANYA
 */

#include "../include/GridBitSet.h"
#include <iostream>
using namespace std;
Grid_BitSet::Grid_BitSet(gridParams gp):Grid() {
	// TODO Auto-generated constructor stub


	string dataFileName;
	filePrefix=gp.FILE_PREFIX;
	ostringstream os;


	os <<"../OUTPUT/"<<filePrefix<<"-grid_data.dat";
	dataFileName = os.str();
	os.str("");

	ofstream gridDataFile(dataFileName.c_str());

	dim=gp.DIM;
	periodic=gp.PERIODIC;
	sortieOKinf=gp.SORTIE_OK_INF;
	sortieOKsup=gp.SORTIE_OK_SUP;

	int d=0;
	unboundedDomain=false;

	while (d<dim && !unboundedDomain)
	{
		unboundedDomain |= (sortieOKinf[d]==1) | (sortieOKsup[d]==1);
		d++;
	}

	nbOMPThreads=gp.OMP_THREADS;

	//cout<< "  grid bit set nb threads= "<<nbOMPThreads<<endl;
	twoPowers=new unsigned long long int[31];
	twoPowers[0]=1;
	for(int i=1;i<31;i++)
	{
		twoPowers[i]=twoPowers[i-1]*2;
	}

	dirs.clear();
	values.clear();

	for(int k=0;k<dim;k++)
	{
		if(gp.SLICE_DIRS[k])
		{
			dirs.push_back(k);
			values.push_back(gp.SLICE_VALUES[k]);
		}
	}


	// les trames sont toujours parall�les � l'axe 1 (x1)
	// donc ce ssont les trames qui  sont stock�es sous la forme d'un tableau de cadinaux
	// les autres axes sont  repr�sent�s juste par des num�ros d leurs points
	//de grilles dans les tabaleaux de tabelaux


	limInf=new double[dim];

	limSup=new double[dim];


	tempIntCoords=new unsigned long long int[dim];

	nbPoints=new unsigned long long int[dim];
	nbPointsSubGrid=new unsigned long long int[dim-1];
	dirTramage=gp.GRID_MAIN_DIR;
	step =new double[dim];
	nbTotalPoints=1;
	nbPointsTotalSubGrid=1;

	  /*
	   * gridType=0 <=> les points  sont les noeuds
	   * gridType=1 <=> les points sont les centres des mailles
	   */
	  gridType=gp.GRID_TYPE;

	maxStep=0;
	for(int i=0;i<dim;i++)
	{
		nbPoints[i]=gp.NBPOINTS[i];
		nbTotalPoints*=nbPoints[i];

		//	px=imx/FacteurMaille;

		limInf[i]= gp.LIMINF[i];
		limSup[i]= gp.LIMSUP[i];

		step[i]=(limSup[i]-limInf[i])/(nbPoints[i]- 1);
		maxStep=max(maxStep, step[i]);
	}

	// cout<<  " nb total de points de la grile "<<  nbTotalPoints<< "   nb points de grille complet = ";
	// printVector(nbPoints, dim);

	// cout<<  "step grille     = ";
	// printVector(step, dim);


	for(int i=0;i<dirTramage;i++)
	{
		nbPointsSubGrid[i]=gp.NBPOINTS[i];
		nbPointsTotalSubGrid*=nbPoints[i];
	}

	for(int i=dirTramage+1;i<dim;i++)
	{

		nbPointsSubGrid[i-1]=gp.NBPOINTS[i];
		nbPointsTotalSubGrid*=nbPoints[i];
	}
	/*  cout<<"dir de tramage "<<dirTramage<<"\n";
  cout<<  " nb total de points de la sous-grile "<<  nbPointsTotalSubGrid<< "   nb points de la  sous-grille  = ";
  printVector(nbPointsSubGrid, dim-1);*/

	longTrame=gp.NBPOINTS[dirTramage];
	li_trame=gp.LIMINF[dirTramage];
	ls_trame=gp.LIMSUP[dirTramage];


	nbTotalCells=1;
	nbCells=new  unsigned long long int[dim];

	for(int i=0; i<dim;i++)
	{
		nbTotalCells*=(nbPoints[i]-1);
		nbCells[i]=(nbPoints[i]-1);
	}


	nbPointsCube=(int) pow(2.0,dim);//pow(2.0, dim);

	indicesDecalCell=new   long long int[nbPointsCube];

	lesDecalagesCell=new unsigned long long int* [nbPointsCube];

	for(int k=0;k<nbPointsCube;k++)
	{
		lesDecalagesCell[k]=new unsigned long long int[dim];
	}

	pow3=1;
	for(int p=0;p<dim;p++)
	{
		pow3*=3;
	}
	indicesDecal=new   long long int[pow3];

	lesDecalagesAxes=new unsigned long long int *[dim];
	indicesDecalAxes=new unsigned long long int[dim];
	for(int k=0;k<dim;k++)
	{
		lesDecalagesAxes[k]=new unsigned long long int[dim];
	}

	computeGridShifts();

	nbCellsSub=new  unsigned long long int[dim-1];
	nbPointsCubeSub=(int) pow(2.0,dim-1);

	indicesDecalCellSub=new  long long int[nbPointsCubeSub];

	lesDecalagesCellSub=new unsigned long long int* [nbPointsCubeSub];
	for(int k=0;k<nbPointsCubeSub;k++)
	{
		lesDecalagesCellSub[k]=new unsigned long long int[dim-1];
	}

	pow3Sub=1;
	for(int p=0;p<dim-1;p++)
	{
		pow3Sub*=3;
	}
	indicesDecalSub=new   long long int[pow3Sub];
	lesDecalagesAxesSub=new unsigned long long int *[dim-1];
	indicesDecalAxesSub=new unsigned long long int[dim-1];

	for(int k=0;k<dim-1;k++)
	{
		lesDecalagesAxesSub[k]=new unsigned long long int[dim-1];
	}
	computeSubGridShifts();
	gridTab= new boost::dynamic_bitset<>*[nbPointsTotalSubGrid];

	for(unsigned long long int i=0;i<nbPointsTotalSubGrid;i++)
	{
		gridTab[i]=(new boost::dynamic_bitset<>(longTrame));
	}


	gridTabNew= new boost::dynamic_bitset<>*[nbPointsTotalSubGrid];

	for(unsigned long long int i=0;i<nbPointsTotalSubGrid;i++)
	{
		gridTabNew[i]=(new boost::dynamic_bitset<>(longTrame));
	}

	gridDataFile<< dim<<endl;
	for(int k=0;k<dim;k++)
	{
		gridDataFile<< limInf[k]<<endl;
	}
	for(int k=0;k<dim;k++)
	{
		gridDataFile<< limSup[k]<<endl;
	}
	for(int k=0;k<dim;k++)
	{
		gridDataFile<< nbPoints[k]<<endl;
	}

	for(int k=0;k<dim;k++)
	{
		gridDataFile<< gp.SLICE_DIRS[k]<<endl;
	}

	for(int k=0;k<dim;k++)
	{
		gridDataFile<< gp.SLICE_VALUES[k]<<endl;
	}

	gridDataFile.close();

}


void Grid_BitSet::copyGrid(  boost::dynamic_bitset<>    **grIn,   boost::dynamic_bitset<>    **grOut)
{
#pragma omp parallel for num_threads(nbOMPThreads)  shared( grIn, grOut) default(none)
	for( unsigned long long int posX=0;posX< nbPointsTotalSubGrid;posX++)
	{
		*grOut[posX]=*grIn[posX];
	}
}

void Grid_BitSet::findNearestViabPointInCell(double *startPoint, double * currentPoint, double * newPoint, double (*dynConstraints)(double *, double *))
{
	unsigned long long int  posTemp,cellNum=localizePoint(currentPoint);
	double minDist=PLUS_INF, dist;
	double testV[dim];
	double oldPoint[dim];
	for(int i=0;i<dim;i++)
	{
		oldPoint[i]=startPoint[i];
	}
	unsigned long long int testI[dim];
	for(int ii=0;ii<nbPointsCube;ii++  )
	{
		posTemp= cellNum+indicesDecalCell[ii];


		numToIntAndDoubleCoords( posTemp ,testI,testV);
		if( (*dynConstraints)(oldPoint, testV)<PLUS_INF)
		{
			if(isInSet(testI))
			{
				dist=0.0;
				for(int i=0;i<dim;i++)
				{
					dist+=(testV[i]-currentPoint[i])*(testV[i]-currentPoint[i]);
				}
				if(dist<minDist)
				{
					for(int i=0;i<dim;i++)
					{
						newPoint[i]=testV[i];
					}
				}
			}
		}
	}
}

void Grid_BitSet::refine()
{

	/*
	 * cette m�thode r�alise une allocation dynamique  de m�moire pour le stockage de l'ensemble
	 * sur une trame raffin�e.
	 * Le pas  sur toutes les directions est divis� par deux
	 * des points sont ajout�s
	 */
	unsigned long long int nbgrille=nbPointsTotalSubGrid;//le nombre de trames dans l'ancien maillage
	unsigned long long int pos;
	unsigned long long int indiceTemp[dim-1];
	unsigned long long int parite[dim-1];

#pragma omp parallel for num_threads(nbOMPThreads) private(pos) shared( nbgrille) default(none)

	for(pos=0;pos<nbgrille;pos++)
	{
		int lTrameNew=2*nbPoints[dirTramage]-1;
		boost::dynamic_bitset<> trameTemp(lTrameNew);
		boost::dynamic_bitset<> trameDecalee1(lTrameNew);
		boost::dynamic_bitset<> trameDecalee2(lTrameNew);
		boost::dynamic_bitset<> trameDecalee3(lTrameNew);
		boost::dynamic_bitset<> trameDecalee4(lTrameNew);
		boost::dynamic_bitset<> trameDecalee5(lTrameNew);
		boost::dynamic_bitset<> trameDecalee6(lTrameNew);
		trameTemp.reset();
		//cout<<"les   anciens  ="<<tramage[pos] <<"\n";
		for(unsigned int l=0;l<longTrame;l++)
		{
			trameTemp[2*l]=(*gridTab[pos])[l];
		}
		 // cout<<"les ralonges des anciens sans decalage="<<trameTemp<<"\n";
		trameDecalee1=trameTemp>>(1);
		trameDecalee2=trameTemp<<(1);
		trameDecalee3=trameTemp>>(2);
		trameDecalee4=trameTemp<<(2);
		trameDecalee5=trameTemp>>(3);
		trameDecalee6=trameTemp<<(3);

		trameTemp|=(trameDecalee1);
		trameTemp|=(trameDecalee2);
		trameTemp|=(trameDecalee3);
		trameTemp|=(trameDecalee4);
		trameTemp|=(trameDecalee5);
		trameTemp|=(trameDecalee6);
		gridTab[pos]->resize(lTrameNew);
		// cout<<"les ralonges des anciens avec decalage="<<trameTemp<<"\n";
		*gridTab[pos]=trameTemp;
	}

	unsigned int dir=this->dirTramage;

	//op�rateur d'incr�mentation lexicographique

	int lTrameNew=2*nbPoints[dir]-1;
	boost::dynamic_bitset<> trameTemp(lTrameNew);
	longTrame=lTrameNew;
	trameTemp.reset();
	// cout<< " suppression de grid tab new ....";
	for(unsigned long long int i=0;i<nbPointsTotalSubGrid;i++)
	{
		delete  gridTabNew[i];
	}
	delete [] gridTabNew;
//	 cout<< " ok \n";

	unsigned long long int *nbPointsSubGridNew=new unsigned long long int[dim-1];
	unsigned long long int totalPointsSubgridNew=1;
	for(int kl=0;kl<dim-1;kl++)
	{
		nbPointsSubGridNew[kl]=2*nbPointsSubGrid[kl]-1;
		totalPointsSubgridNew*=nbPointsSubGridNew[kl];
	}
	// cout<< " allocation grid tab new.......";
	gridTabNew= new boost::dynamic_bitset<>*[totalPointsSubgridNew];
#pragma omp parallel for num_threads(nbOMPThreads)   shared( totalPointsSubgridNew) default(none)
	for(unsigned long long int i=0;i<totalPointsSubgridNew;i++)
	{
		gridTabNew[i]=(new boost::dynamic_bitset<>(longTrame));
		gridTabNew[i]->reset();
	}
	// cout<< " fini\n";

	// cout<< " nbGrille "<<nbgrille<<endl;

//	cout<< " nbPoints "; printVector(nbPoints, dim);
#pragma omp parallel for num_threads(nbOMPThreads)  private(pos) shared( nbgrille,nbPointsSubGridNew) default(none)

	for(pos=0;pos<nbgrille;pos++)
	{
		unsigned long long int indice[dim-1];
		unsigned long long int posNew;

		this->numToSubGridCoords(pos, indice);

		for(int kl=0;kl<dim-1;kl++)
		{
			indice[kl]=2*indice[kl];
		}
		intCoordsToNum_gen(dim-1, nbPointsSubGridNew, indice, &posNew);

		*gridTabNew[posNew]=*gridTab[pos];

	}
	 //cout<< " suppression de grid tab ancien....";
	for(unsigned long long int i=0;i<nbPointsTotalSubGrid;i++)
	{
		//cout<< " i = "<<i<< " trame= "<<*gridTab[i]<<endl;
		delete  gridTab[i];
	}
	delete [] gridTab;

	//cout<< " fini\n";
	// cout<<"ici rafinage size de grid tab "<<nbPointsTotalSubGrid<<" \n";
	for(int i=0;i<dim-1;i++)
	{
		nbPointsSubGrid[i]=nbPointsSubGrid[i]*2-1;
	}

	for(int i=0;i<dim;i++)
	{
		nbPoints[i]=nbPoints[i]*2-1;
	}

	maxStep=0;
	nbTotalPoints=1;
	for(int i=0;i<dim;i++)
	{
		nbTotalPoints*=nbPoints[i];
		step[i]=(limSup[i]-limInf[i])/(nbPoints[i]- 1);
		maxStep=max(maxStep, step[i]);
	}

	  //cout<<  " Rafine:  nb total de points de la grile "<<  nbTotalPoints<< "   nb points de grille complet = ";
	 // printVector(nbPoints, dim);

	 //cout<<  "step grille     = ";
	  //printVector(step, dim);

	nbPointsTotalSubGrid=1;
	for(int i=0;i<dirTramage;i++)
	{
		nbPointsTotalSubGrid*=nbPoints[i];
	}

	for(int i=dirTramage+1;i<dim;i++)
	{
		nbPointsTotalSubGrid*=nbPoints[i];
	}
	/* cout<<"dir de tramage "<<dirTramage<<"\n";
  cout<<  " nb total de points de la sous-grile "<<  nbPointsTotalSubGrid<< "   nb points de la  sous-grille  = ";
  printVector(nbPointsSubGrid, dim-1);
	 */

	// cout<< " allocation de gridTab \n";

	gridTab= new boost::dynamic_bitset<>*[nbPointsTotalSubGrid];

	for(unsigned long long int i=0;i<nbPointsTotalSubGrid;i++)
	{

		gridTab[i]=(new boost::dynamic_bitset<>(longTrame));
		*gridTab[i]= *gridTabNew[i];
	}

	computeGridShifts();
	computeSubGridShifts();
	//  cout<< "  compute grid shifts ok\n";
	pos=0;

	//traitement des impaires (nouveaux) � partir des paires (anciens)
	while( pos< nbPointsTotalSubGrid)
	{
		this->numToSubGridCoords(pos, indiceTemp);
		unsigned long long   int maxParite=0;
		for(int k=0;k<dim-1 ;k++)
		{
			parite[k]=indiceTemp[k]%2;
			maxParite=max(maxParite,parite[k]);
			if(maxParite==1)
				break;
		}
		if(maxParite==1)
		{
			refineTrameMasque(pos,0);
		}
		pos++;
	}
	pos=0;

	//traitement des pairs � partir des impaires
	while( pos< nbPointsTotalSubGrid)
	{
		this->numToSubGridCoords(pos, indiceTemp);
		unsigned long long   int maxParite=0;
		for(int k=0;k<dim-1 ;k++)
		{
			parite[k]=indiceTemp[k]%2;
			maxParite=max(maxParite,parite[k]);
			if(maxParite==1)
				break;
		}
		if(maxParite==0)
		{
			refineTrameMasque(pos,1);
		}
		pos++;
	}
	//traitement des impaires � partir des paires

	pos=0;

	while( pos< nbPointsTotalSubGrid)
	{
		this->numToSubGridCoords(pos, indiceTemp);
		unsigned long long   int maxParite=0;
		for(int k=0;k<dim-1 ;k++)
		{
			parite[k]=indiceTemp[k]%2;
			maxParite=max(maxParite,parite[k]);
			if(maxParite==1)
				break;

		}

		if(maxParite==1)
		{
			refineTrameMasque(pos,0);

		}

		pos++;
	}

}


unsigned long long int  Grid_BitSet::getNbPointsTotalSubGrid()
{
	return nbPointsTotalSubGrid;
}
unsigned long long int * Grid_BitSet::getNbPointsSubGrid()
{
	return nbPointsSubGrid;
}
Grid_BitSet::~Grid_BitSet() {

	for(unsigned long long int i=0;i<nbPointsTotalSubGrid;i++)
	{
		delete  gridTab[i];
	}
	delete[] gridTab;

	if(nbOMPThreads>1)
	{
		for(unsigned long long int i=0;i<nbPointsTotalSubGrid;i++)
		{
			delete  gridTabNew[i];
		}
		delete [] gridTabNew;
	}


}

void Grid_BitSet::printGrid(void)
{
	// cout<< " La grille ici est BitSet!\n";
}
void Grid_BitSet::addPointToSet(unsigned long long int pos, double val)
{
	this->numToIntCoords(pos,tempIntCoords);
	unsigned long long int posX;
	this->intCoordsToNum_dm1(tempIntCoords, &posX);
	this->setPoint(posX, tempIntCoords[dirTramage],true);

}
void Grid_BitSet::savePointsList(string fileName)
{
	// cout<< " couou\n";
}
void Grid_BitSet::addPointToSet(unsigned long long int * coords, double val)
{
	setPoint( coords ,true);
}

void Grid_BitSet::setPoint(unsigned long long int * coords, bool val)
{
	unsigned long long int posX;

	this->intCoordsToNum_dm1(coords, &posX);

	(*gridTab[posX])[coords[dirTramage]]=true;

}

void Grid_BitSet::setPoint(unsigned long long int posX, unsigned long long int posTrame, bool val)
{



	(*gridTab[posX])[posTrame]=val;
	//system("PAUSE");

}

unsigned long long int  Grid_BitSet::getDirTram()
{
	return dirTramage;
}

boost::dynamic_bitset<>  Grid_BitSet::analyseTrameMasqueBis(  unsigned long long int posX, bool print)
						{
	print=false;
	int i=0;

	unsigned long long int indice[dim-1];
	numToIntCoords_gen(posX,  dim-1,  nbPointsSubGrid,indice);
	if(print )
	{
		cout<< " posX= "<<posX<< " size "<<nbPointsTotalSubGrid<<endl;
	}
	if(print)
	{
		cout<<"\n *******************************************************************\n";

		cout<<" \n analyse dirTramage "<<dirTramage<<" indice size = "<<dim-1<<"  indice=";
		for(int jj=0;jj<dim-1;jj++)
		{
			cout<<" "<<indice[jj];
		}
		cout<<"+++++++\n";
	}
	boost::dynamic_bitset<> laTrame(longTrame);
	unsigned long long int iFront;
	unsigned long long int iNext;
	boost::dynamic_bitset<> masqueDecale(longTrame);

	boost::dynamic_bitset<> masque1(longTrame);
	boost::dynamic_bitset<> masque2(longTrame);
	boost::dynamic_bitset<> masque;
	//cout<<"\n";
	bool testBord=false;
	i=0;
	while((i<dim-1) && !testBord)
	{
		testBord=testBord|((indice[i]==nbPointsSubGrid[i]-1 )| (indice[i]==0));
		i++;
	}


	if(testBord)
	{
		if(print)  cout<<"sur la fronti�re direct \n";
		laTrame=*gridTab[posX];
		return laTrame;
	}
	else
	{

		if(print) cout<<  "pas sur la frontiere  pos X = " <<posX <<endl;
		laTrame=*gridTab[posX];

		if(print )
		{
			cout<<  " la trame                   "<< laTrame<<endl;

		}
		// on test si la trame n'est pas vide
		if(laTrame.none())
		{   if(print) cout<<"trame vide "<<laTrame;
		return laTrame;
		}
		else
		{
			// cout<< " trame pas vide "<<laTrame<<endl;
			int k=0;
			while((posX+indicesDecalSub[k]<0) & (k<pow3Sub))
			{
				if(print)  cout<< " k ="<<k<< " indices dcal="<<indicesDecalSub[k]<<endl;
				k++;
			}


			if(print)
			{
				int kk=k;
				while(kk<pow3Sub && posX+indicesDecalSub[kk]<nbPointsTotalSubGrid)
				{
					cout<<" voisin  k= "<<kk<<" trame voisine "<<((*gridTab[posX+indicesDecalSub[kk]]))<<"\n";
					kk++;
				}

			}


			if((k<pow3Sub) & (posX+indicesDecalSub[k]<nbPointsTotalSubGrid))
			{
				masque=(*gridTab[posX+indicesDecalSub[k]]);

				while(k<pow3Sub)
				{
					if((posX+indicesDecalSub[k]>=0) & (posX+indicesDecalSub[k]<nbPointsTotalSubGrid))
					{
						masque1=(((*gridTab[posX+indicesDecalSub[k]])));
						//masque1&=(((*gridTab[posX+indicesDecalSub[k])));
						//on d�cale la trame de 1 bite  et on fait le ET pour chaque
						//bit �a donne : b[i] ET b[i+1]
						//    cout<<" voisin  k= "<<k<<" size de trame voisine "<<((*gridTab[k)).size()<<"\n ";//<<((*gridTab[k))<<"\n";
						masqueDecale=masque1>>(1);

						masqueDecale[longTrame-1]=1;
						masque&=(masqueDecale);
						if(print)
						{
							cout<< " k= "<<k<<endl;
							cout<<" masque1 "<<masque1<<"\n";
							cout<<" masqueD "<<masqueDecale<<"\n";
							cout<<" masque  "<<masque<<"\n";
							cout<<"\n ============================================\n";

						}
						//on d�cale la trame de 1 bite  et on fait le ET pour chaque
						//bit �a donne : b[i] ET b[i-1]
						masqueDecale=masque1<<(1);

						masqueDecale[0]=1;

						masque&=(masqueDecale);
						if(print)
						{
							cout<<" masque1 "<<masque1<<"\n";
							cout<<" masqueD "<<masqueDecale<<"\n";
							cout<<" masque  "<<masque<<"\n";
							cout<<"\n ============================================\n";

						}
					}
					k++;
				}
			}

			masque^=(laTrame);
			if(print)
			{
				cout<<" masque  chap tram "<<masque<<"\n  ++++++++++++++++++++++\n";
			}
			masque&=(laTrame);
			if(print)
			{
				cout<<" masque   ET tram  "<<masque<<"\n  ++++++++++++++++++++++\n";
			}

			iFront=laTrame.find_first();

			masque[iFront]=1;

			while(iFront<longTrame)
			{
				iNext=laTrame.find_next(iFront);
				while(iNext-iFront==1 && iNext<longTrame)
				{
					iFront=iNext;
					iNext=laTrame.find_next(iFront);
				}
				if(iNext<longTrame)
				{
					masque[iFront]=1;
					masque[iNext]=1;
					iFront=iNext;
				}
				else
				{
					masque[iFront]=1;
					iFront=iNext;
				}
			}


			if(print) cout<< " masque final "<<masque<<endl;
			return(masque);}
	}
						}

boost::dynamic_bitset<>  Grid_BitSet::analyseTrameMasque( unsigned long long int posX)
						{

	int i=0;
	double x= limInf[0]+step[0];
	unsigned long long int indice[dim-1];
	numToIntCoords_gen(posX,  dim-1,  nbPointsSubGrid,indice);
	/*  cout<<"analyse dirTramage "<<dirTramage<<"indice size= "<<dim-1<<"indice=";
        for(int jj=0;jj<dim-1;jj++)
        {
                cout<<" "<<indice[jj];
        }
        cout<<"\n";*/
	boost::dynamic_bitset<> laTrame(longTrame);
	unsigned long long int iFront;
	unsigned long long int iNext;
	boost::dynamic_bitset<> masqueDecale(longTrame);
	boost::dynamic_bitset<> masque;
	//deuxi�me  test : si la trame n'est pas sur la fornti�re de l'espace


	//cout<<"\n";
	bool testBord=false;
	i=0;
	while((i<dim-1) && !testBord)
	{
		testBord=testBord|((indice[i]==nbPointsSubGrid[i]-1 )| (indice[i]==0));
		i++;
	}


	if(testBord)
	{
		//cout<<"sur la fronti�re direct \n";
		return((*gridTab[posX]));
	}
	else
	{

		// cout<< posX <<endl;
		laTrame=(*gridTab[posX]);
		// on test si la trame n'est pas vide
		if(laTrame.none())
		{ // cout<<"trame vide "<<laTrame;
			return laTrame;
		}
		else
		{
			int k=0;
			while((posX+indicesDecalSub[k]<0) & (k<pow3Sub))
			{
				//cout<< " k ="<<k<< " indices dcal="<<indicesDecalSub[k]<<endl;
				k++;
			}
			if((k<pow3Sub) & (posX+indicesDecalSub[k]<nbPointsTotalSubGrid))
			{
				masque=(*gridTab[posX+indicesDecalSub[k]]);

				//on d�cale la trame de 1 bite  et on fait le ET pour chaque
				//bit �a donne : b[i] ET b[i+1]
				masqueDecale=((*gridTab[posX+indicesDecalSub[k]]))>>(1);

				masqueDecale[longTrame-1]=1;
				masque&=(masqueDecale);

				//on d�cale la trame de 1 bite  et on fait le ET pour chaque
				//bit �a donne : b[i] ET b[i-1]
				masqueDecale=((*gridTab[posX+indicesDecalSub[k]]))<<(1);

				masqueDecale[0]=1;
				masque&=(masqueDecale);

				while(k<pow3Sub)
				{
					if(posX+indicesDecalSub[k]>=0)
					{
						masque&=(((*gridTab[posX+indicesDecalSub[k]])));

						//on d�cale la trame de 1 bite  et on fait le ET pour chaque
						//bit �a donne : b[i] ET b[i+1]
						masqueDecale=((*gridTab[posX+indicesDecalSub[k]]))>>(1);
						masqueDecale[longTrame-1]=1;
						masque&=(masqueDecale);
						//on d�cale la trame de 1 bite  et on fait le ET pour chaque
						//bit �a donne : b[i] ET b[i-1]
						masqueDecale=((*gridTab[posX+indicesDecalSub[k]]))<<(1);
						masqueDecale[0]=1;
						masque&=(masqueDecale);

					}
					k++;
				}
			}

			masque^=(laTrame);
			masque&=(laTrame);

			iFront=laTrame.find_first();

			masque[iFront]=1;

			while(iFront<longTrame)
			{
				iNext=laTrame.find_next(iFront);
				while(iNext-iFront==1 && iNext<longTrame)
				{
					iFront=iNext;
					iNext=laTrame.find_next(iFront);
				}
				if(iNext<longTrame)
				{
					masque[iFront]=1;
					masque[iNext]=1;
					iFront=iNext;
				}
				else
				{
					masque[iFront]=1;
					iFront=iNext;
				}
			}
			return(masque);}
	}
						}



/*
 * parite == 0 indices paires
 * paritee==1 indices impaires
 */
void  Grid_BitSet::refineTrameMasque( unsigned long long int posX, int pariteVoisin)
{
	//recherche de points de frotni�re dans la trame

	unsigned long long int parite[dim-1];
	unsigned long long int indice[dim-1];
	unsigned long long int maxParite=0;
	boost::dynamic_bitset<> masque;
	 // cout<< posX <<endl;
	masque=(*gridTab[posX]);
	//cout<<masque<<endl;
	// on test si la trame n'est pas vide

	int k=0;
	while((posX+indicesDecalSub[k]<0) & (k<pow3Sub))
	{
		// cout<< " k ="<<k<< " indices dcal="<<indicesDecalSub[k]<<endl;
		k++;
	}
	if((k<pow3Sub) & (posX+indicesDecalSub[k]<nbPointsTotalSubGrid))
	{
		numToIntCoords_gen(posX+indicesDecalSub[k],  dim-1,  nbPointsSubGrid,indice);
		maxParite=0;
		for(int k=0;k<dim-1 ;k++)
		{
			parite[k]=indice[k]%2;
			maxParite=max(maxParite,parite[k]);
			if(maxParite==1)
				break;
		}
		if((int)maxParite==pariteVoisin)
		{
			masque|=*gridTab[posX+indicesDecalSub[k]];
		}
		while(k<pow3Sub)
		{
			if(posX+indicesDecalSub[k]>=0 && posX+indicesDecalSub[k]<nbPointsTotalSubGrid)
			{
				numToIntCoords_gen(posX+indicesDecalSub[k],  dim-1,  nbPointsSubGrid,indice);
				maxParite=0;
				for(int k=0;k<dim-1 ;k++)
				{
					parite[k]=indice[k]%2;
					maxParite=max(maxParite,parite[k]);
					if(maxParite==1)
						break;
				}
				if((int)maxParite==pariteVoisin)
				{
					masque|=(*(gridTab[posX+indicesDecalSub[k]]));
				}
			}
			k++;
		}
	}

	*gridTab[posX]=masque;


}
/*
 * Transformation  de num�rotation alpha-num�rique :
 *  � partir du num�ro des coordonn�es enti�res du points dans la grille
 *  sont num�ro
 */

void Grid_BitSet:: loadSet(string fileName)
{
	string line;
	ifstream * userDataFile= new ifstream();
	double xCoords[dim], val;
	istringstream sstr;


	userDataFile->open(fileName);
	if (userDataFile->good())
	{
		for( unsigned long long  int posX=0;posX<nbPointsTotalSubGrid;posX++)
		{
			for(unsigned long long  int k=0;k<longTrame;k++)
			{
				getline((*userDataFile), line);
				line.append(" ");
				sstr.str(line);

				for(int i=0;i<dim;i++)
				{
					sstr>>xCoords[i];
				}
				sstr>>val;
				(*gridTab[posX])[k]=val;
			}
		}
	}
}

void Grid_BitSet::saveValOnGrid(string fileName)
{
	//cout<<"ecriture  de l'ensemble dans un fichier \n";
	// instructions


	unsigned long long int   indice[dim-1];
	double xCoordsDouble[dim];
	ofstream fichierB(fileName.c_str());

	if(fichierB)  // si l'ouverture a r�ussi
	{

		//        cout<<"dim etat  " <<dim<< " taille de trame est " <<nbPointsTotalSubGrid<<"\n";

		for( unsigned long long  int posX=0;posX<nbPointsTotalSubGrid;posX++)
		{
			//cout<< " posx = "<<posX<<endl;
			numToIntCoords_gen(posX,  dim-1,  nbPointsSubGrid,indice);


			for( int j=0; j<dirTramage;j++)
			{
				xCoordsDouble[j]=limInf[j]+indice[j]*step[j];
			}

			for(int j=dirTramage+1; j<(int)dim;j++)
			{
				xCoordsDouble[j]=limInf[j]+indice[j-1]*step[j];
			}

			for(unsigned long long int k=0;k<longTrame;k++)
			{
				xCoordsDouble[dirTramage]=limInf[dirTramage]+k*step[dirTramage];

				if((*gridTab[posX])[k])
				{
					for(int l1=0;l1<dim;l1++)
					{
						fichierB<<   xCoordsDouble[l1]<<" ";
					}
					fichierB<<"1.0 \n";
				}
				else
				{
					for(int l1=0;l1<dim;l1++)
					{
						fichierB<<   xCoordsDouble[l1]<<" ";
					}
					fichierB<<"0.0 \n";
				}
			}
		}// fin de for de parcours de la trame
		fichierB.close();
		// je referme le fichier
	}
	else  // sinon
		cerr << "Erreur � l'ouverture !" << endl;
}
void Grid_BitSet::intCoordsToNum_dm1( unsigned long long int * coords, unsigned long long int * res)
{
	if(dirTramage>0)
	{
		(*res)=coords[0];
		int i=0;
		while(i<(int)dirTramage-1)
		{
			(*res)=(*res)*nbPointsSubGrid[i+1]+coords[i+1];
			i++;
		}
		i=dirTramage;
		while(i<dim-1)
		{
			//cout <<"  i="<<i<<endl;
			(*res)=(*res)*nbPointsSubGrid[i]+coords[i+1];
			i++;
		}
		//cout<< "  res = "<<*res<<endl;
	}
	else
	{
		(*res)=coords[1];
		//cout<< " dir tramage est  "<<dirTramage<<endl;

		//cout<< "  res = "<<*res<<endl;
		int i=1;

		while(i<dim-1)
		{
			//cout <<"  i="<<i<<endl;
			(*res)=(*res)*nbPointsSubGrid[i]+coords[i+1];
			i++;
		}
		//cout<< "  res = "<<*res<<endl;
	}
}
void Grid_BitSet::numToIntAndDoubleCoords_dm1(unsigned long long int num,unsigned long long int *resI, double * resD)
{
	unsigned long long  int temp=num;

	for(  int d=(int)dim-1;d>=dirTramage+1;d--)
	{
		resI[d]=temp%nbPoints[d];  // coordonn�es enti�res du point
		resD[d]=limInf[d]+step[d]*resI[d]+0.5*gridType*step[d];

		temp=(temp-resI[d])/nbPoints[d];
	}
	for(  int d=(int)dirTramage-1;d>=0;d--)
	{
		resI[d]=temp%nbPoints[d];  // coordonn�es enti�res du point
		resD[d]=limInf[d]+step[d]*resI[d]+0.5*gridType*step[d];

		temp=(temp-resI[d])/nbPoints[d];
	}

}
inline void Grid_BitSet::numToIntCoords_dm1(unsigned long long int num,unsigned long long int *res)
{



	unsigned long long  int temp=num;
	// cout<< "num= "<<num<<" to int coords nbPoints "; printVector(nbPoints, dim);

	for( int d=dim-1;d>dirTramage;d--)
	{
		res[d]=temp%nbPoints[d];                       // coordonn�es enti�res du point
		temp=(temp-res[d])/nbPoints[d];
		//   cout<< " temp= "<<temp<< " res ="<<res[d]<<endl;

	}
	for( int d=dirTramage-1;d>=0;d--)
	{
		res[d]=temp%nbPoints[d];                       // coordonn�es enti�res du point
		temp=(temp-res[d])/nbPoints[d];
		//   cout<< " temp= "<<temp<< " res ="<<res[d]<<endl;

	}

}

inline void Grid_BitSet::numToSubGridCoords(unsigned long long int num,unsigned long long int *res)
{



	unsigned long long  int temp=num;
	// cout<< "num= "<<num<<" to int coords nbPoints "; printVector(nbPoints, dim);

	for( int d=dim-1;d>dirTramage;d--)
	{
		res[d-1]=temp%nbPoints[d];                       // coordonn�es enti�res du point
		temp=(temp-res[d-1])/nbPoints[d];
		//    cout<< " temp= "<<temp<< " res ="<<res[d-1]<<endl;

	}
	for( int d=dirTramage-1;d>=0;d--)
	{
		res[d]=temp%nbPoints[d];                       // coordonn�es enti�res du point
		temp=(temp-res[d])/nbPoints[d];
		//    cout<< " temp= "<<temp<< " res ="<<res[d]<<endl;

	}

}

unsigned long long int Grid_BitSet::getLongTrame()
{
	return longTrame;

}

boost::dynamic_bitset<>  ** Grid_BitSet:: getGridTab()
						{
	return gridTab;
						}


boost::dynamic_bitset<>  ** Grid_BitSet:: getGridTabNew()
						{
	return  gridTabNew;
						}

void Grid_BitSet::saveProjetion(string fileName, unsigned long long int * projection)
{
	//cout<<"projection de la fonction valeur  dans un fichier \n";
	// instructions


	unsigned long long int * reducedNbPoints=new unsigned long long int[dim-1];
	double * reducedLimInf=new double[dim-1];
	double * reducedStep=new double[dim-1];
	unsigned long long int * tempVect=new unsigned long long int[dim-1];

	int i=0;
	int j=0;
	int projAxe=0;
	while(j<dim)
	{
		//cout<< " i="<<i<< " j= "<<j<<endl;
		if(!projection[j])
		{
			reducedNbPoints[i]=nbPoints[j];
			reducedLimInf[i]=limInf[j];
			reducedStep[i]=step[j];
			i++;
			j++;
		}
		else
		{
			projAxe=j;
			j++;
		}
	}
	//cout<< " l'axe de projection est "<<projAxe<<endl;

	unsigned long long int reducedTotalPoints=1;
	for(i=0;i<dim-1;i++)
	{
		reducedTotalPoints*=reducedNbPoints[i];
	}

	unsigned long long int * x = new unsigned long long int[dim] ;


	ofstream fichierB(fileName.c_str());

	if(fichierB)  // si l'ouverture a r�ussi
	{
		int   compteB=0;
		unsigned long long int k=0;
		bool test=false;
		for(unsigned long long int posX=0;posX<reducedTotalPoints;posX++)
		{

			numToIntCoords_gen(posX, dim-1, reducedNbPoints,tempVect);
			for(j=0;j<projAxe;j++)
			{
				x[j]=tempVect[j];
			}
			for(j=projAxe+1;j<dim;j++)
			{
				x[j]=tempVect[j-1];
			}

			k=0;
			test=false;
			while((k<nbPoints[projAxe]) & !test)
			{
				x[projAxe]=k;
				if(this->isInSet(x))
				{
					test=true;
				}
				k++;
			}
			if(test)
			{
				compteB++;

				for(int l1=0;l1<dim-1;l1++)
				{
					fichierB<< reducedLimInf[l1]+reducedStep[l1]*tempVect[l1]<<" ";
				}
				fichierB<<"\n";
			}
		}
		cout<< " nbPoints  d'ensemble  ecrits="<<compteB<<"\n";
		fichierB.close();
	}
	else  // sinon
		cerr << "Erreur de4 l'ouverture !" << endl;


	//cout<<"fichier fini\n";


	delete[] reducedNbPoints;
	delete[] reducedLimInf;
	delete[] reducedStep;
	delete[] x;
	delete [] tempVect;

}

void Grid_BitSet::saveCoupe(string nomFichier)
{


	unsigned long long int pos;
	double  coordReelles[dim];
	unsigned long long int coordsInt[dim];
	ofstream fichier(nomFichier.c_str());
	if(fichier)  // si l'ouverture a r�ussi
	{
		int compte=0;
		bool test;
		pos=0;
		while(  pos<nbTotalPoints)
		{
			numToIntAndDoubleCoords(pos, coordsInt,coordReelles );
			if(this->isInSet(coordsInt))
			{
				compte++;
				//si le point appartient � l'ensemble
				//on  va enregistrer ses corrdonn�es r�elles

				test=true;
				for(int k=0;k<(int)dirs.size();k++)
				{
					test=test&(abs(coordReelles[dirs[k]]-values[k])<=step[dirs[k]]);
				}
				if(test)
				{

					for(int i=0;i<dim;i++)
					{
						fichier<<" "<< coordReelles[i]<<" ";
					}
					fichier<<"\n";
				}
			}
			pos++;
		}
		//cout<<" nbPoints ecrits="<<compte<<"\n";
		fichier.close();  // je referme le fichier
	}
	else  // sinon
		cerr << "Erreur � l'ouverture !" << endl;
	//cout<<"fichier fini\n";
}
void Grid_BitSet::saveCoupeBoundary(string nomFichier)
{
	boost::dynamic_bitset<> masqueFront;
	bool test;

	double  coordReelles[dim];
	unsigned long long int coordInt[dim];

	ofstream fichier(nomFichier.c_str());
	if(fichier)  // si l'ouverture a r�ussi
	{
		int compte=0;
		for(unsigned long long   int posX=0;posX<nbPointsTotalSubGrid;posX++)
		{
			this->numToIntAndDoubleCoords_dm1(posX, coordInt, coordReelles);
			masqueFront=this->analyseTrameMasque(posX);

			for(unsigned long long int k=0;k<longTrame;k++)
			{

				if(masqueFront[k])
				{
					coordReelles[dirTramage]=k*step[dirTramage]+limInf[dirTramage];
					test=true;
					for(int k=0;k<(int)dirs.size();k++)
					{
						test=test&(abs(coordReelles[dirs[k]]-values[k])<=step[dirs[k]]);
					}

					if(test)
					{
						compte++;
						// //cout<<"compt="<<compte<< "\n";
						//si le point appartient � l'ensemble
						//on  va enregistrer ses corrdonn�es r�elles

						for(int i=0;i<dim;i++)
						{
							fichier<<" "<< coordReelles[i]<<" ";

							////cout<<"i= "<<i<<"liminf "<< espace->getLimitesInf()[i] <<"pas "<<espace->getPasDiscretisation()[i]<<"coordrelle "<<CR[i]<<"\n";
						}
						fichier<<"\n";
					}
				}
			}
		}

		fichier.close();
		// je referme le fichier
	}
	else  // sinon
		cerr << "Erreur � l'ouverture !" << endl;
	//cout<<"fichier fini\n";
}
void Grid_BitSet::saveBoundary(string nomFichier )
{

	boost::dynamic_bitset<> masqueFront;
	double  coordReelles[dim];
	unsigned long long int coordInt[dim];

	ofstream fichier(nomFichier.c_str());

	if(fichier)  // si l'ouverture a r�ussi
	{
		int compte=0;
		for(  int posX=0;posX<(int)nbPointsTotalSubGrid;posX++)
		{
			this->numToIntAndDoubleCoords_dm1(posX, coordInt, coordReelles);
			masqueFront=this->analyseTrameMasque(posX);

			for(int k=0;k<(int)longTrame;k++)
			{

				if(masqueFront[k])
				{
					coordReelles[dirTramage]=k*step[dirTramage]+limInf[dirTramage];

					compte++;

					for(int i=0;i<dim;i++)
					{
						fichier<<" "<< coordReelles[i]<<" ";

						////cout<<"i= "<<i<<"liminf "<< espace->getLimitesInf()[i] <<"pas "<<espace->getPasDiscretisation()[i]<<"coordrelle "<<CR[i]<<"\n";
					}
					fichier<<"\n";
				}
			}
		}

		fichier.close();
		// je referme le fichier
	}
	else  // sinon
		cerr << "Erreur � l'ouverture !" << endl;
	//cout<<"fichier fini\n";
}
bool Grid_BitSet::isInSet(unsigned long long int * coords )
{
	unsigned long long int posX;


	this->intCoordsToNum_dm1(coords, &posX);


	return (*gridTab[posX])[coords[dirTramage]];
}

void savePointsList(string fileName)
{

}

void Grid_BitSet::computeSubGridShifts()
{

	//  cout<< " compute subgrid shifts \n";


	nbTotalCellsSub=1;
	for(int i=0; i<dim-1;i++)
	{
		nbTotalCellsSub*=(nbPointsSubGrid[i]-1);
		nbCellsSub[i]=(nbPointsSubGrid[i]-1);
	}

	/*
	 * On calcule les vecteurs  binaires  d�finissant les coordonn�es entieres
	 *  des sommets par rapport au sommet  en bas � gauche
	 */
	int numX;
	vector< boost::dynamic_bitset<> >bb(nbPointsCubeSub),bb1(2*(dim-1));
	boost::dynamic_bitset<> BBtemp(dim-1);

	for(int k=0;k<nbPointsCubeSub;k++)
	{
		bb.at(k)=boost::dynamic_bitset<>(dim-1, k);
		//cout << "bits(k) = " << bb.at(k) << std::endl;

		for(int i=0;i<dim-1;i++)
		{
			lesDecalagesCellSub[k][i]=(int)bb.at(k)[dim-2-i];
		}

		numX=lesDecalagesCellSub[k][0];

		for(int i=0;i<dim-2;i++)
		{
			numX=numX*nbPointsSubGrid[i+1]+lesDecalagesCellSub[k][i+1];
		}
		indicesDecalCellSub[k]=numX;


	}
	/* cout<<" SUB grid: les decalages  dans le cube  sont ";
  for(int ll=0;ll<nbPointsCubeSub;ll++)
    {
    for(int i=0;i<dim-1;i++)
      {
      cout<< " "<<lesDecalagesCellSub[ll][i];

      }
  //  cout<< " numero correspondant " <<indicesDecalCellSub[ll]<<endl;
    }
  cout<<endl;
	 */
	////system("pause");


	int tempInt;
	int *indiceTemp=new int[dim-1];
	// cout<< "  calcul des  voisins d'un point de grilles  \n";
	for(int p=0;p<pow3Sub;p++)
	{
		tempInt=p;
		for(int k=dim-2;k>-1;k--)
		{

			indiceTemp[k]=tempInt%3;
			tempInt=(tempInt-indiceTemp[k])/3;
		}
		numX=indiceTemp[0];
		for(int i=0;i<dim-2;i++)
		{
			numX=numX*nbPointsSubGrid[i+1]+indiceTemp[i+1];

		}
		indicesDecalSub[p]=numX;

	}
	tempInt=indicesDecalSub[(pow3Sub-1)/2];
	for(int p=0;p<pow3Sub;p++)
	{

		indicesDecalSub[p]-=tempInt;
		// cout<<" "<<indicesDecalSub[p]<<" ";

	}
	// cout<< "\n fini le voisinage  entier il reste les axes\n";

	int Pow2=1;
	for(int k=0;k<dim-1;k++)
	{
		bb1.at(k)=boost::dynamic_bitset<>(dim-1, Pow2);
		//	cout << "bits(k) = " << bb.at(k) << std::endl;

		for(int i=0;i<dim-1;i++)
		{
			lesDecalagesAxesSub[k][i]= (int)bb1.at(k)[dim-2-i];
		}
		numX=lesDecalagesAxesSub[k][0];

		for(int i=0;i<dim-2;i++)
		{
			numX=numX*nbPointsSubGrid[i+1]+lesDecalagesAxesSub[k][i+1];
		}

		(indicesDecalAxesSub)[k]=numX;
		// cout<<  "numX= "<<numX<<"\n";

		Pow2*=2;
	}

	/*//cout<< "v�rif des indices de decalage axe";
	for(int k=0;k<dim;k++)
	{
		//cout<<indicesDecalSubAxes[k]<<" ";

	}
	//cout<<endl;
	 */


	/*//cout<< "v�rif des indices de decalage du cube ";
		for(int k=0;k<nbPointsCube;k++)
		{
			//cout<<indicesDecalSubCube.at(k)<<" ";

		}
		//cout<<endl;
		//system("pause");
	 */
	/*time_t start,end;
			//cout<< " debut de cr�ation de vertex\n";
				time (&start);



				time (&endTimeTime);

					dif = difftime (endTime,start);
					printf (" temps d'execution   %.5f secondes.\n", dif );

	 */
}


