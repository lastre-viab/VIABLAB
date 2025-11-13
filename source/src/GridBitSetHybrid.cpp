/*
 * GridBitSetHybrid.cpp
 *
 *  Created on: 17 août 2025
 *      Author: adesi
 */

#include "../include/GridBitSetHybrid.h"

GridBitSetHybrid::~GridBitSetHybrid()
    {
    // TODO Auto-generated destructor stub
    }

GridBitSetHybrid::GridBitSetHybrid(const gridParams &gp) :
						Grid_BitSet(gp)
    {
    dim_hc = gp.DIM_HC;
    dim_hd = gp.DIM_HD;
    maxStep = 0;//recompute here the maxStep only on continuous subgrid
    for (int i = 0; i < dim_hc; i++)
	{
	maxStep = max(maxStep, step[i]);
	}
    nbCellPointsContinuousState = pow(2, dim_hc);
    continuousStateShifts = new unsigned long long int[nbCellPointsContinuousState];
    for (int k = 0; k < (int) nbCellPointsContinuousState; k++)
	{
	continuousStateShifts[k] = indicesDecalCell[k];
	}
    }

unsigned long long int* GridBitSetHybrid::GetContinuousStateShifts()
    {
    return continuousStateShifts;
    }

unsigned long long int GridBitSetHybrid::getNbCellsContinuousState()
    {
    return nbCellPointsContinuousState;
    }

bool GridBitSetHybrid::isHybridPointInGrid(double *continuousState, unsigned long long int *discreteState)
    {
    bool isInGrid = true;
    int i = 0;

    while (isInGrid & (i < dim_hc))
	{
	isInGrid &= ((continuousState[i] <= limSup[i]) & (continuousState[i] >= limInf[i]));
	i++;
	}
    if (isInGrid)
	{
	while (isInGrid & (i < dim))
	    {
	    isInGrid &= (discreteState[i - dim_hc] < nbPoints[i]);
	    i++;
	    }
	}

    return isInGrid;
    }

unsigned long long int GridBitSetHybrid::LocalizeHybridPoint(double *doubleCoords, unsigned long long int *intCoords)
    {
    unsigned long long int numCell = (unsigned long long int) (((doubleCoords)[0] - limInf[0]) / step[0]);
    if (numCell == (nbPoints[0] - 1))
	numCell--;

    for (int i = 1; i < dim_hc; i++)
	{
	unsigned long long int indiceCell = (unsigned long long int) (((doubleCoords)[i] - limInf[i]) / step[i]);
	if (indiceCell == (nbPoints[i] - 1))
	    {
	    indiceCell--;
	    }
	numCell = numCell * (nbPoints[i]) + indiceCell;
	}

    for (int i = dim_hc; i < dim; i++)
	{
	numCell = numCell * (nbPoints[i]) + intCoords[i - dim_hd];
	}
    return numCell;
    }

void GridBitSetHybrid::saveHybridCoupe(const string &nomFichier) const
    {

    unsigned long long int pos;
    double coordReelles[dim];
    unsigned long long int coordsInt[dim];
    ofstream fichier(nomFichier.c_str());
    if (fichier)  // si l'ouverture a r�ussi
	{
	int compte = 0;
	bool test;
	pos = 0;
	while (pos < nbTotalPoints)
	    {
	    numToIntAndDoubleCoords(pos, coordsInt, coordReelles);
	    if (this->isInSet(coordsInt))
		{
		compte++;
		//si le point appartient � l'ensemble
		//on  va enregistrer ses corrdonn�es r�elles

		test = true;
		for (int k = 0; k < (int) dirs.size(); k++)
		    {
		    if (dirs[k] < dim_hc)
			{
			test = test & (abs(coordReelles[dirs[k]] - values[k]) <= step[dirs[k]]);
			}
		    else
			{
			test = test & (coordsInt[dirs[k]] == values_fd[k]);
			}

		    }
		if (test)
		    {

		    for (int i = 0; i < dim; i++)
			{
			fichier << " " << coordReelles[i] << " ";
			}
		    fichier << "\n";
		    }
		}
	    pos++;
	    }
	//cout<<" nbPoints ecrits="<<compte<<"\n";
	fichier.close();  // je referme le fichier
	}
    else
	// sinon
	cerr << "Erreur � l'ouverture !" << endl;
    //cout<<"fichier fini\n";
    }

void GridBitSetHybrid::saveValOnGridHybrid(const string &fileName) const
    {

    unsigned long long int indice[dim - 1];
    double xCoordsDouble[dim];
    unsigned long long int xCoordsInt[dim];
    ofstream fichierB(fileName.c_str());

    if (fichierB)
	{

	for (unsigned long long int posX = 0; posX < nbPointsTotalSubGrid; posX++)
	    {
	    //cout<< " posx = "<<posX<<endl;
	    numToIntCoords_gen(posX, dim - 1, nbPointsSubGrid, indice);

	    for (int j = 0; j < dirTramage; j++)
		{
		xCoordsDouble[j] = limInf[j] + indice[j] * step[j];
		xCoordsInt[j] = indice[j];
		}

	    for (int j = dirTramage + 1; j < (int) dim; j++)
		{
		xCoordsDouble[j] = limInf[j] + indice[j - 1] * step[j];
		xCoordsInt[j] = indice[j-1];
		}

	    for (unsigned long long int k = 0; k < longTrame; k++)
		{
		xCoordsDouble[dirTramage] = limInf[dirTramage] + k * step[dirTramage];
		xCoordsInt[dirTramage] = k;
		if ((*gridTab[posX])[k])
		    {
		    for (int l1 = 0; l1 < dim_hc; l1++)
			{
			fichierB << xCoordsDouble[l1] << " ";
			}
		    for (int l1 = 0; l1 < dim_hd; l1++)
			{
			fichierB << xCoordsInt[dim_hc + l1] << " ";
			}
		    fichierB << "1 \n";
		    }
		else
		    {
		    for (int l1 = 0; l1 < dim_hc; l1++)
			{
			fichierB << xCoordsDouble[l1] << " ";
			}
		    for (int l1 = 0; l1 < dim_hd; l1++)
			{
			fichierB << xCoordsInt[dim_hc + l1] << " ";
			}
		    fichierB << "0 \n";
		    }
		}
	    }  // fin de for de parcours de la trame
	fichierB.close();
	// je referme le fichier
	}
    else
	// sinon
	cerr << "Erreur � l'ouverture !" << endl;
    }

void GridBitSetHybrid::saveValOnGridLightHybrid(const string &fileName) const
    {
    unsigned long long int indice[dim - 1];
    double xCoordsDouble[dim];
    unsigned long long int xCoordsInt[dim];
    ofstream fichierB(fileName.c_str());

    if (fichierB)
	{

	for (unsigned long long int posX = 0; posX < nbPointsTotalSubGrid; posX++)
	    {
	    //cout<< " posx = "<<posX<<endl;
	    numToIntCoords_gen(posX, dim - 1, nbPointsSubGrid, indice);

	    for (int j = 0; j < dirTramage; j++)
		{
		xCoordsDouble[j] = limInf[j] + indice[j] * step[j];
		xCoordsInt[j] = indice[j];
		}

	    for (int j = dirTramage + 1; j < (int) dim; j++)
		{
		xCoordsDouble[j] = limInf[j] + indice[j - 1] * step[j];
		xCoordsInt[j] = indice[j-1];
		}

	    for (unsigned long long int k = 0; k < longTrame; k++)
		{
		xCoordsDouble[dirTramage] = limInf[dirTramage] + k * step[dirTramage];
		xCoordsInt[dirTramage] = k;
		if ((*gridTab[posX])[k])
		    {
		    for (int l1 = 0; l1 < dim_hc; l1++)
			{
			fichierB << xCoordsDouble[l1] << " ";
			}
		    for (int l1 = 0; l1 < dim_hd; l1++)
			{
			fichierB << xCoordsInt[dim_hc + l1] << " ";
			}
		    fichierB << "1 \n";
		    }
		}
	    }  // fin de for de parcours de la trame
	fichierB.close();
	// je referme le fichier
	}
    else
	// sinon
	cerr << "Erreur � l'ouverture !" << endl;
    }

void GridBitSetHybrid::loadSetHybrid(const string &fileName)
    {
    string line;
    ifstream userDataFile(fileName);
    double xCoords[dim_hc];
    bool val;
    unsigned long long int xCoordsInt[dim];
    istringstream sstr;
    int d = 0;

    if (!userDataFile)
	{
	spdlog::error("Error loading set, could not open set file : {}", fileName);
	}
    else
	{
	while (getline(userDataFile, line))
	    {
	    sstr.str(line);
	    for (d = 0; d < dim_hc; d++)
		{
		sstr >> xCoords[d];
		// Transformation inverse des coordonnées en double vers les coordonnées entières
		const double xCoordsIntI = (xCoords[d] - limInf[d]) / step[d] - 0.5 * arePointsGridCenters;
		xCoordsInt[d] = llround(xCoordsIntI);
		}
	    for(d = dim_hc; d< dim; d++)
		{
		sstr >> xCoordsInt[d];
		}
	    sstr >> val;

	    unsigned long long int posX = this->intCoordsToNum_dm1(xCoordsInt);
	    (*gridTab[posX])[xCoordsInt[dirTramage]] = val;
	    }
	}
    }
