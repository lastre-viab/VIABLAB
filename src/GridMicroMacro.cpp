/*
 * GridHJB.cpp
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
 *
 *      Author: ANYA
 */

#include "../include/GridMicroMacro.h"

GridMicroMacro::GridMicroMacro() {
  // TODO Auto-generated constructor stub

}
GridMicroMacro::GridMicroMacro(   gridParams gp) {
  // TODO Auto-generated constructor stub

  //	 cout<< " ***************************************************************\n";
  if(DEV_PRINT) cout<< " cretion de grille micro-macro stockage tableau\n";

  string dataFileName;
  ostringstream os;

  filePrefix=gp.FILE_PREFIX;

  os <<"OUTPUT/"<<filePrefix<<"-grid_data.dat";
  dataFileName = os.str();
  os.str("");

  ofstream gridDataFile(dataFileName.c_str());
  nbOMPThreads=gp.OMP_THREADS;

  limInf=gp.LIMINF;
  limSup=gp.LIMSUP;
  dim=gp.DIM;


  nbPoints=gp.NBPOINTS;
  periodic=gp.PERIODIC;



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


  /*
   * gridType=0 <=> les points  sont les noeuds
   * gridType=1 <=> les points sont les centres des mailles
   */
  gridType=gp.GRID_TYPE;

  step=new double[dim];
  if(DEV_PRINT)  cout<< " le pas est \n";
  int d;
  maxStep=0;
  nbTotalPoints=1;
  nbTotalCells=1;
  nbCells=new unsigned long long int[dim];
  for(d=0;d<dim;d++)
    {
    step[d]=(limSup[d]-limInf[d])/(nbPoints[d]-(1-gridType)*1);
    if(DEV_PRINT) 	 cout<< " "<<step[d];
    maxStep=max(maxStep, step[d]);
    nbTotalPoints*=nbPoints[d];
    nbTotalCells*=(nbPoints[d]-1);
    nbCells[d]=(nbPoints[d]-1);

    }
  if(DEV_PRINT) cout<<endl;

  gridPtr=new double [nbTotalPoints];

  vectUnsigIntTemp=new unsigned long long int [dim];
  vectInt=new int unsigned long long [dim];


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

  if(DEV_PRINT) cout<< " construction du grille micro-macro finie:OK\n";

}


GridMicroMacro::~GridMicroMacro() {
  // TODO Auto-generated destructor stub
  //fermeture de la base
  delete [] gridPtr;
  cout<<" grid detruit \n";

}

void GridMicroMacro::loadSet(string fileName)
{
  string line;
  ifstream * userDataFile= new ifstream();
  double xCoords[dim], val;
  istringstream sstr;


  userDataFile->open(fileName);
  if (userDataFile->good())
    {
    for (unsigned long long int pos=0;pos<nbTotalPoints; pos++)
        {
        getline((*userDataFile), line);
        line.append(" ");
        sstr.str(line);

        for(int i=0;i<dim;i++)
          {
          sstr>>xCoords[i];
          }
        sstr>>val;
        gridPtr[pos]=val;

      }
    }
}

void GridMicroMacro::printGrid(void)
{
  cout<< " La grille ici est HJB!\n ********************\n\n ***************************\n";
  cout<< " nbPoints de grid ";
  printVector(nbPoints, dim);
  cout<< " la dimention est "<<dim<<endl;

  cout<< " lim inf ";
  printVector(limInf, dim);

  cout<< " lim  sup";
  printVector(limSup, dim);


  cout<< " La grille ici est HJB!\n ********************\n\n ***************************\n";

  //system("pause");
}


double * GridMicroMacro::getGridPtr()
{
  return gridPtr;
}

double * GridMicroMacro::getGridPtr_tmp()
{
  return gridPtr_tmp;
}


void GridMicroMacro::copyGrid(  double *grIn,  double *grOut)
{
#pragma omp parallel for num_threads(nbOMPThreads)  shared( grIn, grOut) default(none)
  for( unsigned long long int posX=0;posX< nbTotalPoints;posX++)
    {
    grOut[posX]=grIn[posX];
    }
}
bool GridMicroMacro::isInSet(unsigned long long int * coords )
{
  //cout<< " coucouc\n";
  return true;
}
void GridMicroMacro::savePointsList(string fileName)
{
  //cout<<"ecriture  de l'ensemble dans un fichier \n";
  unsigned long long int * x = new unsigned long long int[dim] ;
  double *xReel=new double[dim];

  ofstream fichierB(fileName.c_str());

  if(fichierB)  // si l'ouverture a r�ussi
    {
    int   compteB=0;
    for (unsigned long long int pos=0;pos<nbTotalPoints; pos++)
      {
      numToIntAndDoubleCoords(pos, x,xReel);
      if( (gridPtr[pos])<(double)PLUS_INF)
        {
        compteB++;

        for(int l1=0;l1<dim;l1++)
          {
          fichierB<<   xReel[l1]<<" ";
          }
        fichierB<<"\n";
        }
      }
    fichierB.close();
    }
  else  // sinon
    cerr << "Erreur de  l'ouverture !" << endl;
}


void GridMicroMacro::saveValOnGrid(string fileName)
{
  //cout<<"ecriture  de l'ensemble dans un fichier \n";
  unsigned long long int * x = new unsigned long long int[dim] ;
  double *xReel=new double[dim];

  ofstream fichierB(fileName.c_str());

  if(fichierB)  // si l'ouverture a r�ussi
    {
    int   compteB=0;
    for (unsigned long long int pos=0;pos<nbTotalPoints; pos++)
      {
      numToIntAndDoubleCoords(pos, x,xReel);

      compteB++;

      for(int l1=0;l1<dim;l1++)
        {
        fichierB<<   xReel[l1]<<" ";
        }
      fichierB<<gridPtr[pos]<<"\n";
      }

    fichierB.close();
    }
  else  // sinon
    cerr << "Erreur de  l'ouverture !" << endl;
}

void  GridMicroMacro::saveSubLevelset(double level, string fileName )
{
  //cout<<"ecriture  de l'ensemble dans un fichier \n";
  unsigned long long int * x = new unsigned long long int[dim] ;
  double *xReel=new double[dim];

  ofstream fichierB(fileName.c_str());

  if(fichierB)  // si l'ouverture a r�ussi
    {
    int   compteB=0;
    for (unsigned long long int pos=0;pos<nbTotalPoints; pos++)
      {
      if(gridPtr[pos]<=level)
        {
        numToIntAndDoubleCoords(pos, x,xReel);

        compteB++;

        for(int l1=0;l1<dim;l1++)
          {
          fichierB<<   xReel[l1]<<" ";
          }
        fichierB<<gridPtr[pos]<<"\n";
        }
      }

    fichierB.close();
    }
  else  // sinon
    cerr << "Erreur de  l'ouverture !" << endl;
}

void GridMicroMacro::computeMinMaxValues(double &minV, double & maxV)
{
  minV=PLUS_INF;
  maxV=-PLUS_INF;

  for (unsigned long long int pos=0;pos<nbTotalPoints; pos++)
    {
    maxV=max(maxV, gridPtr[pos]);
    minV=min(maxV, gridPtr[pos]);
    }
}


void GridMicroMacro::saveProjetion(string fileName, unsigned long long int * projection)
{
  //cout<<"projection de la fonction valeur  dans un fichier \n";
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
  double val;
  unsigned long long int * x = new unsigned long long int[dim] ;
  ofstream fichierB(fileName.c_str());

  if(fichierB)  // si l'ouverture a r�ussi
    {
    unsigned long long int pointNum;
    int   compteB=0;
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

      val=PLUS_INF;
      for(unsigned long long int k=0;k<nbPoints[projAxe];k++)
        {
        x[projAxe]=k;
        /*//cout<< " x= ";
				for(unsigned long long int m=0;m<dim;m++)
				{
					//cout<< x[m]<<" ";
				}
				//cout<< " \n";*/
        this->intCoordsToNum(x,&pointNum);
        ////cout<< " nmero "<< pointNum<<endl;
        val=min(val, gridPtr[pointNum]);
        if(val<PLUS_INF-1)
          {
          compteB++;

          for(int l1=0;l1<dim-1;l1++)
            {
            fichierB<< tempVect[l1]<<" ";
            }
          fichierB<<val<<"\n";
          }
        }
      //cout<< " nbPoints  d'ensemble  ecrits="<<compteB<<"\n";
      fichierB.close();
      // je referme le fichier
      }
    }
  else  // sinon
    cerr << "Erreur � l'ouverture !" << endl;
}

void GridMicroMacro::addPointToSet(unsigned long long int * coords, double val)
{

  unsigned long long int pos;

  Grid::intCoordsToNum(coords, &pos);
  gridPtr[pos]=val;
}

void GridMicroMacro::addPointToSet(unsigned long long int pos, double val)
{
  gridPtr[pos]=val;
}
