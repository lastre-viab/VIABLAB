/*
 * Grid.cpp
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
 *  Created on: 10 sept. 2013
 *      Author: Ana DESILLES
 */

#include "../include/Grid.h"

Grid::Grid()
    {

    }

Grid::~Grid()
    {
    // TODO Auto-generated destructor stub
    // cout<< " delete grid \n";
    //cout<< "  destructeur de grid generique \n";

    delete[] step; /*!< \brief  pas de discretisation pour X , tableau	 */

    delete[] vectUnsigIntTemp; /*!< \brief une m�moire tempon pour les calculs ,tableau	 */
    delete[] vectInt; /*!< \brief une m�moire tempon pour les calculs ,tableau	 */

    delete[] nbCells;

    delete[] indicesDecalCell;

    delete[] indicesDecalAxes;
    delete[] indicesDecal;

    int nbPointsCube = (int) pow(2.0, dim);

    for (int k = 0; k < nbPointsCube; k++)
	{
	delete[] lesDecalagesCell[k];
	}
    delete[] lesDecalagesCell;

    for (int k = 0; k < dim; k++)
	{
	delete[] lesDecalagesAxes[k];
	}
    delete[] lesDecalagesAxes;

    // cout<< "  fini  destructeur de grid generique \n";
    }

void Grid::printGrid()
    {

    }
bool Grid::isInSet(unsigned long long int *coords)
    {
    return false;
    }

unsigned long long int Grid::getNearestPointInSet(double *coords)
    {
    return 0;
    }

void Grid::savePointsList(string fileName)
    {

    }
void Grid::saveValOnGrid(string fileName)
    {

    }
void Grid::saveValOnGridLight(string fileName)
    {

    }
unsigned long long int Grid::getDim()
    {
    return dim;
    }
void Grid::numToIntCoords(unsigned long long int num,
	unsigned long long int *res)
    {

    unsigned long long int temp = num;

    for (int d = dim - 1; d >= 0; d--)
	{
	res[d] = temp % nbPoints[d];
	temp = (temp - res[d]) / nbPoints[d];

	}

    }
void Grid::numToIntAndDoubleCoords(unsigned long long int num,
	unsigned long long int *resI, double *resD)
    {
    bool print = false;  //(num == 5100);
    unsigned long long int temp = num;
    for (int d = (int) dim - 1; d >= 0; d--)
	{
	if (print)
	    cout << "  temp = " << temp << " nbPoints[d] = " << nbPoints[d];
	resI[d] = temp % nbPoints[d];  // coordonn�es enti�res du point
	if (print)
	    cout << " resI = " << resI[d];
	resD[d] = limInf[d] + step[d] * resI[d] + 0.5 * gridType * step[d];
	if (print)
	    cout << "  gridType = " << gridType << " resD = " << resD[d];
	temp = (temp - resI[d]) / nbPoints[d];
	if (print)
	    cout << "  temp = " << temp << endl;
	}

    }
/*
 * Transformation  de num�rotation alpha-num�rique :
 *  � partir du num�ro des coordonn�es enti�res du points dans la grille
 *  sont num�ro
 */
void Grid::intCoordsToNum(unsigned long long int *coords,
	unsigned long long int *res)
    {
    (*res) = coords[0];
    for (int i = 0; i < dim - 1; i++)
	{
	(*res) = (*res) * nbPoints[i + 1] + coords[i + 1];

	}
    }

/*!
 *  Cette  fonction d�termine la maille qui contient le point de coordonn�es r�elles donn�es
 *  Important: n suppose ici  que les coordonn�es se  trouvent bien dans les limites d�finissant la grille
 *  Il est donc n�cessaire de v�rifier avant l'apartenance du point � la grille
 * @param coords
 * @return num�ro de la cellule  qui contient le point r�el
 */
unsigned long long int Grid::localizePoint(double *coords)
    {
    unsigned long long int *indiceCell = new unsigned long long int[dim];

    //int *indiceCell =vectInt;
    // cout<< " indice Cell dans localize ";
    for (int k = 0; k < dim; k++)
	{
	indiceCell[k] = (unsigned long long int) (floor(
		(((coords)[k] - limInf[k]) / step[k])));

	if (indiceCell[k] == (nbPoints[k] - 1))
	    {
	    // cout<<  " coucouc"<<endl;
	    indiceCell[k]--;

	    }
	// cout<< " "<<indiceCell[k];
	}
    // cout<< endl;
    unsigned long long int numCell = indiceCell[0];

    for (int i = 0; i < dim - 1; i++)
	{
	numCell = numCell * (nbPoints[i + 1]) + indiceCell[i + 1];

	}
    delete[] indiceCell;
    return numCell;
    }

bool Grid::ArePointsInTheSameCell(double *coords1, double *coords2)
    {
    bool res = true;
    for (int k = 0; k < dim; k++)
	{
	res &= (abs(coords1[k] - coords2[k]) <= this->step[k]);
	}
    return res;
    }

void Grid::computeGridShifts()
    {

//  cout<< "  compute grid shifts dans Grid.cpp\n";
    nbTotalCells = 1;

    for (int i = 0; i < dim; i++)
	{
	nbTotalCells *= (nbPoints[i] - 1);
	nbCells[i] = (nbPoints[i] - 1);
	}

    /*
     * On calcule les vecteurs  binaires  d�finissant les coordonn�es entieres
     *  des sommets par rapport au sommet  en bas � gauche
     */
    unsigned long long int numX;
    vector<boost::dynamic_bitset<> > bb(nbPointsCube), bb1(2 * dim);
    boost::dynamic_bitset<> BBtemp(dim);

    for (int k = 0; k < nbPointsCube; k++)
	{
	bb.at(k) = boost::dynamic_bitset<>(dim, k);
	//cout << "bits(k) = " << bb.at(k) << std::endl;

	for (int i = 0; i < dim; i++)
	    {
	    lesDecalagesCell[k][i] = (int) bb.at(k)[dim - 1 - i];
	    }

	numX = lesDecalagesCell[k][0];

	for (int i = 0; i < dim - 1; i++)
	    {
	    numX = numX * nbPoints[i + 1] + lesDecalagesCell[k][i + 1];
	    }
	indicesDecalCell[k] = numX;

	}
    /* cout<< " nbPointsCube = "<<nbPointsCube<< " dim = "<<dim<<endl;
     cout<<" grid: les decalages  dans le cube  sont ";
     for(int ll=0;ll<nbPointsCube;ll++)
     {
     for(int i=0;i<dim;i++)
     {
     cout<< " "<<lesDecalagesCell[ll][i];

     }
     cout<< " numero correspondant " <<indicesDecalCell[ll]<<endl;
     }
     cout<<endl;*/

    ////system("pause");
    int tempInt;
    int *indiceTemp = new int[dim];
    ////cout<< "  calcul des  voisins d'un point de grilles  \n";
    for (int p = 0; p < pow3; p++)
	{
	tempInt = p;
	for (int k = dim - 1; k > -1; k--)
	    {

	    indiceTemp[k] = tempInt % 3;
	    tempInt = (tempInt - indiceTemp[k]) / 3;
	    }
	numX = indiceTemp[0];
	for (int i = 0; i < dim - 1; i++)
	    {
	    numX = numX * nbPoints[i + 1] + indiceTemp[i + 1];

	    }
	indicesDecal[p] = numX;

	}
    tempInt = indicesDecal[(pow3 - 1) / 2];
    for (int p = 0; p < pow3; p++)
	{

	indicesDecal[p] -= tempInt;
	//   cout<<" "<<indicesDecal[p]<<" ";

	}
//  cout<< "\n fini le voisinage  entier il reste les axes\n";

    int Pow2 = 1;
    for (int k = 0; k < dim; k++)
	{
	bb1.at(k) = boost::dynamic_bitset<>(dim, Pow2);
	//	cout << "bits(k) = " << bb.at(k) << std::endl;

	for (int i = 0; i < dim; i++)
	    {
	    lesDecalagesAxes[k][i] = (int) bb1.at(k)[dim - 1 - i];
	    }
	numX = lesDecalagesAxes[k][0];

	for (int i = 0; i < dim - 1; i++)
	    {
	    numX = numX * nbPoints[i + 1] + lesDecalagesAxes[k][i + 1];

	    }

	(indicesDecalAxes)[k] = numX;
	//cout<<  "numX= "<<numX<<"\n";

	Pow2 *= 2;
	}
    delete[] indiceTemp;

    /*//cout<< "v�rif des indices de decalage axe";
     for(int k=0;k<dim;k++)
     {
     //cout<<indicesDecalAxes[k]<<" ";

     }
     //cout<<endl;
     */

    // cout<< " decal shifts grid.cpp fini\n";
    /*time_t start,end;
     //cout<< " debut de cr�ation de vertex\n";
     time (&start);



     time (&endTimeTime);

     dif = difftime (endTime,start);
     printf (" temps d'execution   %.5f secondes.\n", dif );

     */
    }

/*!
 * cette fonction corrige les coordonn�es r�elles d'un vecteur
 * qui sont p�riodiques en les ramenant le cas �ch�ant  dans l'intervalle de la p�riode
 */
void Grid::periodizePoint(double *vect)
    {
    int i;
    double dist;
    bool testDepass;
    ////cout<< " peridisation  vect=";
    //printVector(vect,dim);

    for (i = 0; i < dim; i++)
	{
	/*
	 * On teste  si la variable est d�crar�e p�riodique
	 */
	////cout<< " i= "<< i<< " periodic["<<i<<"]= "<<periodic[i]<<endl;
	if (periodic[i])
	    {

	    testDepass = true;
	    /*
	     * La periode
	     */
	    dist = limSup[i] - limInf[i];
	    ////cout<< " dist= "<<dist<<endl;
	    /*
	     * Correction par p�riode enti�re
	     */
	    while (testDepass)
		{
		//	//cout<< " test depass="<<testDepass<<endl;
		if (vect[i] > limSup[i])
		    {
		    vect[i] -= dist;
		    }
		else
		    {
		    if (vect[i] < limInf[i])
			{
			vect[i] += dist;
			}
		    else
			{
			testDepass = false;
			}
		    }
		//		//cout<< " vect= ";
		//		printVector(vect,dim);
		}

	    }
	}
    //	//cout<< "Periodise a fini  vect= ";
    //				printVector(vect,dim);
    ////system("pause");
    }

/*!
 *  Cette  fonction v�rifie l'apartenance  d'un point  dont les coordonnees sont donn�es
 *  au pav�  definissant la grille
 * @param coords coordonnees reelles d'un point
 * @return isInGrid : l'indicateur booleen
 */
bool Grid::isPointInGrid(double *coords)
    {

    bool isInGrid = true;
    int i = 0;

    while (isInGrid & (i < dim))
	{
	isInGrid &= ((coords[i] <= limSup[i]) & (coords[i] >= limInf[i]));
	i++;
	}
    return isInGrid;

    }

bool Grid::isPointInGrid_fd(unsigned long long int *coords)
    {

    bool isInGrid = true;
    int i = 0;

    while (isInGrid & (i < dim))
	{
	isInGrid &= ((coords[i] < nbPoints[i]) & (coords[i] >= 0));
	i++;
	}
    return isInGrid;

    }

unsigned long long int Grid::getNbTotalCells()
    {
    return nbTotalCells;
    }

unsigned long long int Grid::getNbTotalPoints()
    {
    // cout<< " ici grid nb total points est "<<nbTotalPoints<<endl;
    return nbTotalPoints;
    }

double Grid::getMaxStep()
    {
    return maxStep;
    }

long long int* Grid::getIndicesDecalCell()
    {
    return indicesDecalCell;
    }

unsigned long long int* Grid::getNbPoints()
    {
    return nbPoints;
    }

bool Grid::isPointInGridWithConstr(double *coords)
    {

    bool isInGrid = true;
    int i = 0;

    while (isInGrid & (i < dim))
	{
	isInGrid &= (((coords[i] <= limSup[i]) | (sortieOKsup[i]))
		& ((coords[i] >= limInf[i]) | (sortieOKinf[i])));
	i++;
	}
    return isInGrid;

    }

