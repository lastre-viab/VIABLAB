/*
 * GridMicroMacro.cpp
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
 *  Created on:  2018
 *      Author: Anna DESILLES
 */

#include "../include/GridMicroMacro.h"

GridMicroMacro::GridMicroMacro()
    {
    // TODO Auto-generated constructor stub

    }
GridMicroMacro::GridMicroMacro(const gridParams &gp)
    {

    string dataFileName;
    ostringstream os;

    filePrefix = gp.FILE_PREFIX;

    os << "../OUTPUT/" << filePrefix << "-grid_data.dat";
    dataFileName = os.str();
    os.str("");

    ofstream gridDataFile(dataFileName.c_str());
    nbOMPThreads = gp.OMP_THREADS;

    limInf = gp.LIMINF;
    limSup = gp.LIMSUP;
    dim = gp.DIM;

    nbPoints = gp.NBPOINTS;
    periodic = gp.PERIODIC;

    dirs.clear();
    values.clear();
    values_fd.clear();
    for (int k = 0; k < dim; k++)
	{
	if (gp.SLICE_DIRS[k])
	    {
	    dirs.push_back(k);
	    values.push_back(gp.SLICE_VALUES[k]);
	    values_fd.push_back(gp.SLICE_VALUES_FD[k]);
	    }
	}

    sortieOKinf = gp.SORTIE_OK_INF;
    sortieOKsup = gp.SORTIE_OK_SUP;

    int d = 0;
    unboundedDomain = false;

    while (d < dim && !unboundedDomain)
	{
	unboundedDomain |= (sortieOKinf[d] == 1) | (sortieOKsup[d] == 1);
	d++;
	}
    spdlog::debug("[GridMicroMacro] : check unbounded domain : {}", unboundedDomain);

    gridDataFile << dim << endl;
    for (int k = 0; k < dim; k++)
	{
	gridDataFile << limInf[k] << endl;
	}
    for (int k = 0; k < dim; k++)
	{
	gridDataFile << limSup[k] << endl;
	}
    for (int k = 0; k < dim; k++)
	{
	gridDataFile << nbPoints[k] << endl;
	}

    for (int k = 0; k < dim; k++)
	{
	gridDataFile << gp.SLICE_DIRS[k] << endl;
	}

    for (int k = 0; k < dim; k++)
	{
	gridDataFile << gp.SLICE_VALUES[k] << endl;
	}
    gridDataFile.close();

    /*
     * arePointsGridCenters=0 <=> les points  sont les noeuds
     * arePointsGridCenters=1 <=> les points sont les centres des mailles
     */
    arePointsGridCenters = gp.GRID_TYPE;

    step = new double[dim];

    d = 0;
    maxStep = 0;
    nbTotalPoints = 1;
    nbTotalCells = 1;
    nbCells = new unsigned long long int[dim];
    for (d = 0; d < dim; d++)
	{
	step[d] = (limSup[d] - limInf[d]) / (nbPoints[d] - (1 - arePointsGridCenters) * 1);

	maxStep = max(maxStep, step[d]);
	nbTotalPoints *= nbPoints[d];
	nbTotalCells *= (nbPoints[d] - 1);
	nbCells[d] = (nbPoints[d] - 1);

	}
    logVector("[GridMicroMacro] : grid step : ", step, dim);

    try
	{
	gridPtr = new double[nbTotalPoints];
	}
    catch (...)
	{
	cout << " error while allocating ";
	throw;
	}

    if (nbOMPThreads > 1)
	{
	try
	    {
	    gridPtr_tmp = new double[nbTotalPoints];
	    }
	catch (...)
	    {
	    cout << " error while allocating ";
	    throw;
	    }
	}

    vectUnsigIntTemp = new unsigned long long int[dim];
    vectInt = new int unsigned long long[dim];

    nbPointsCube = (int) pow(2.0, dim); //pow(2.0, dim);

    indicesDecalCell = new long long int[nbPointsCube];

    lesDecalagesCell = new unsigned long long int*[nbPointsCube];

    for (int k = 0; k < nbPointsCube; k++)
	{
	lesDecalagesCell[k] = new unsigned long long int[dim];
	}

    pow3 = 1;
    for (int p = 0; p < dim; p++)
	{
	pow3 *= 3;
	}
    indicesDecal = new long long int[pow3];

    lesDecalagesAxes = new unsigned long long int*[dim];
    indicesDecalAxes = new unsigned long long int[dim];
    for (int k = 0; k < dim; k++)
	{
	lesDecalagesAxes[k] = new unsigned long long int[dim];
	}

    computeGridShifts();

    spdlog::info("[GridMicorMacro] : grid initialization finished");

    }

GridMicroMacro::~GridMicroMacro()
    {
    delete[] gridPtr;
    if (nbOMPThreads > 1)
	{
	delete[] gridPtr_tmp;
	}
    spdlog::debug("[GridMicroMacro] : grid destructor finished");

    }

void GridMicroMacro::loadSet(const string &fileName)
    {
    string line;
    ifstream userDataFile(fileName);
    double val;
    istringstream sstr;
    double *xCoords = new double[dim];

    if (userDataFile.good())
	{
        while (getline(userDataFile, line)) {
            sstr.str(line);

            unsigned long long int pos = 0;
            for (int d = 0; d < dim; ++d) {
                sstr >> xCoords[d];

                // Transformation double vers int...
                const double invDoubleCoord =
                    (xCoords[d] - limInf[d])/step[d] - 0.5*arePointsGridCenters;                
                const unsigned long long int coordInt = llround(invDoubleCoord);
                // ...puis int vers num
                pos = pos * nbPoints[d] + coordInt;
            }            
            sstr >> val;
            gridPtr[pos] = val;
        }
	}
    delete[] xCoords;
    }

void GridMicroMacro::printGrid(void) const
    {
    spdlog::debug("[GridMicorMacro] : Grid type : Micro-Macro");
    spdlog::debug("[GridMicorMacro] : Dimension = {}", dim);
    logVector("[GridMicorMacro] : number of points per axis : ", nbPoints, dim);
    logVector("[GridMicorMacro] : inf limits : ", limInf, dim);
    logVector("[GridMicorMacro] : sup limits : ", limSup, dim);
    }

const double* GridMicroMacro::getGridPtr() const
    {
    return gridPtr;
    }

double* GridMicroMacro::getGridPtr()
    {
    return gridPtr;
    }

const double* GridMicroMacro::getGridPtr_tmp() const
    {
    return gridPtr_tmp;
    }

void GridMicroMacro::copyGrid(const double *grIn, double *grOut) const
    {
    unsigned long long int posX = 0;
#pragma omp parallel for num_threads(nbOMPThreads) private (posX) shared( grIn, grOut) default(none)
    for (posX = 0; posX < nbTotalPoints; posX++)
	{
	grOut[posX] = grIn[posX];
	}
    }
bool GridMicroMacro::isInSet(const unsigned long long int *coords) const
    {
    unsigned long long int posX;
    posX = this->intCoordsToNum(coords);
    return ((gridPtr[posX]) < (double) PLUS_INF);
    }

unsigned long long int GridMicroMacro::getNearestPointInSet(const double *coords) const
    {
    unsigned long long int nearest = this->nbTotalPoints + 1; // means that teh is no near points that are viable
    double minDist = PLUS_INF;
    unsigned long long int cellNum = localizePoint(coords);
    unsigned long long int posTemp;
    unsigned long long int testI[dim];
    double testV[dim];
    for (int ii = 0; ii < nbPointsCube; ii++)
	{
	posTemp = cellNum + indicesDecalCell[ii];

	numToIntAndDoubleCoords(posTemp, testI, testV);
	if (gridPtr[posTemp] < PLUS_INF)
	    {
	    double dist = 0.0;
	    for (int k = 0; k < dim; k++)
		{
		dist = max(dist, abs(testV[k] - coords[k]));
		}
	    if (dist < minDist)
		{
		minDist = dist;
		nearest = posTemp;
		}
	    }
	}
    return nearest;
    }

unsigned long long int GridMicroMacro::getBestNearPointInSet(const double *coords) const
    {
    unsigned long long int nearest = this->nbTotalPoints + 1; // means that teh is no near points that are viable
    double minVal = PLUS_INF;
    unsigned long long int cellNum = localizePoint(coords);
    unsigned long long int posTemp;
    for (int ii = 0; ii < nbPointsCube; ii++)
	{
	posTemp = cellNum + indicesDecalCell[ii];
	if (gridPtr[posTemp] < minVal)
	    {
	    minVal = gridPtr[posTemp];
	    nearest = posTemp;
	    }
	}
    return nearest;
    }

double GridMicroMacro::getOptimalValue(const double *coords) const
    {
    double minCellVal = PLUS_INF;
    unsigned long long int cellNum = localizePoint(coords);
    for (int ii = 0; ii < nbPointsCube; ii++)
	{
	minCellVal = min(minCellVal, gridPtr[cellNum + indicesDecalCell[ii]]);
	}
    return minCellVal;
    }

void GridMicroMacro::savePointsList(const string &fileName) const
    {
    //cout<<"ecriture  de l'ensemble dans un fichier \n";
    unsigned long long int *x = new unsigned long long int[dim];
    double *xReel = new double[dim];

    ofstream fichierB(fileName.c_str());

    if (fichierB)  // si l'ouverture a r�ussi
	{
	int compteB = 0;
	for (unsigned long long int pos = 0; pos < nbTotalPoints; pos++)
	    {
	    numToIntAndDoubleCoords(pos, x, xReel);
	    if ((gridPtr[pos]) < (double) PLUS_INF)
		{
		compteB++;

		for (int l1 = 0; l1 < dim; l1++)
		    {
		    fichierB << xReel[l1] << " ";
		    }
		fichierB << "\n";
		}
	    }
	fichierB.close();
	}
    else
	{
	spdlog::error("[GridMicorMacro] : Error while open file {}", fileName);
	}

    delete[] x;
    delete[] xReel;
    }

void GridMicroMacro::saveValOnGrid(const string &fileName) const
    {
    //cout<<"ecriture  de l'ensemble dans un fichier \n";
    unsigned long long int *x = new unsigned long long int[dim];
    double *xReel = new double[dim];

    ofstream fichierB(fileName.c_str());

    if (fichierB)  // si l'ouverture a r�ussi
	{
	int compteB = 0;
	for (unsigned long long int pos = 0; pos < nbTotalPoints; pos++)
	    {
	    numToIntAndDoubleCoords(pos, x, xReel);

	    compteB++;

	    for (int l1 = 0; l1 < dim; l1++)
		{
		fichierB << xReel[l1] << " ";
		}
	    fichierB << gridPtr[pos] << "\n";
	    }

	fichierB.close();
	}
    else
	{
	spdlog::error("[GridMicorMacro] : Error while open file {}", fileName);
	}
    delete[] x;
    delete[] xReel;
    }

void GridMicroMacro::saveValOnGridLight(const string &fileName) const
    {
    cout << "ecriture  de l'ensemble light  dans un fichier \n";
    unsigned long long int *x = new unsigned long long int[dim];
    double *xReel = new double[dim];

    ofstream fichierB(fileName.c_str());

    if (fichierB)  // si l'ouverture a r�ussi
	{
	int compteB = 0;
	for (unsigned long long int pos = 0; pos < nbTotalPoints; pos++)
	    {
	    if (gridPtr[pos] < PLUS_INF)
		{
		numToIntAndDoubleCoords(pos, x, xReel);

		compteB++;

		for (int l1 = 0; l1 < dim; l1++)
		    {
		    fichierB << xReel[l1] << " ";
		    }
		fichierB << gridPtr[pos] << "\n";
		}
	    }

	fichierB.close();
	}
    else
	{
	spdlog::error("[GridMicorMacro] : Error while open file {}", fileName);
	}
    delete[] x;
    delete[] xReel;
    }

void GridMicroMacro::saveCoupeBoundary(const string &fileName) const
    {
    bool test;

    //cout<<"ecriture  de l'ensemble dans un fichier \n";
    unsigned long long int *x = new unsigned long long int[dim];
    double *xReel = new double[dim];

    ofstream fichierB(fileName.c_str());
    if (fichierB)  // si l'ouverture a r�ussi
	{
	int compteB = 0;
	for (unsigned long long int pos = 0; pos < nbTotalPoints; pos++)
	    {
	    numToIntAndDoubleCoords(pos, x, xReel);
	    test = true;
	    for (int k = 0; k < (int) dirs.size(); k++)
		{
		test = test & (abs(xReel[dirs[k]] - values[k]) <= step[dirs[k]]);
		}
	    if (test)
		{
		compteB++;

		for (int l1 = 0; l1 < dim; l1++)
		    {
		    fichierB << xReel[l1] << " ";
		    }
		fichierB << gridPtr[pos] << "\n";
		}

	    }

	fichierB.close();
	}
    else
	{
	spdlog::error("[GridMicorMacro] : Error while open file {}", fileName);
	}
    delete[] x;
    delete[] xReel;
    }

void GridMicroMacro::saveCoupeBoundary_DD(const string &fileName) const
    {
    bool test;

    //cout<<"ecriture  de l'ensemble dans un fichier \n";
    unsigned long long int *x = new unsigned long long int[dim];
    double *xReel = new double[dim];

    ofstream fichierB(fileName.c_str());
    if (fichierB)  // si l'ouverture a r�ussi
	{
	int compteB = 0;
	for (unsigned long long int pos = 0; pos < nbTotalPoints; pos++)
	    {
	    numToIntAndDoubleCoords(pos, x, xReel);
	    test = true;
	    for (int k = 0; k < (int) dirs.size(); k++)
		{
		test = test & (x[dirs[k]] == values_fd[k]);
		}
	    if (test)
		{
		compteB++;

		for (int l1 = 0; l1 < dim; l1++)
		    {
		    fichierB << x[l1] << " ";
		    }
		fichierB << gridPtr[pos] << "\n";
		}

	    }

	fichierB.close();
	}
    else
	{
	spdlog::error("[GridMicorMacro] : Error while open file {}", fileName);
	}
    }

void GridMicroMacro::saveValOnGrid_DD(const string &fileName) const
    {
    unsigned long long int *xInt = new unsigned long long int[dim];
    double *xReel = new double[dim];

    ofstream fichierB(fileName.c_str());

    if (fichierB)
	{
	int compteB = 0;
	for (unsigned long long int pos = 0; pos < nbTotalPoints; pos++)
	    {
	    numToIntAndDoubleCoords(pos, xInt, xReel);

	    compteB++;

	    for (int l1 = 0; l1 < dim; l1++)
		{
		fichierB << xInt[l1] << " ";
		}
	    fichierB << gridPtr[pos] << "\n";
	    }

	fichierB.close();
	}
    else
	{
	spdlog::error("[GridMicorMacro] : Error while open file {}", fileName);
	}
    }

void GridMicroMacro::saveSubLevelset(double level, const string &fileName) const
    {
    unsigned long long int *x = new unsigned long long int[dim];
    double *xReel = new double[dim];

    ofstream fichierB(fileName.c_str());

    if (fichierB)
	{
	int compteB = 0;
	for (unsigned long long int pos = 0; pos < nbTotalPoints; pos++)
	    {
	    if (gridPtr[pos] <= level)
		{
		numToIntAndDoubleCoords(pos, x, xReel);

		compteB++;

		for (int l1 = 0; l1 < dim; l1++)
		    {
		    fichierB << xReel[l1] << " ";
		    }
		fichierB << gridPtr[pos] << "\n";
		}
	    }

	fichierB.close();
	}
    else
	{
	spdlog::error("[GridMicorMacro] : Error while open file {}", fileName);
	}
    }
void GridMicroMacro::saveSubLevelset_DD(double level, const string &fileName) const
    {
    unsigned long long int *xInt = new unsigned long long int[dim];
    double *xReel = new double[dim];

    ofstream fichierB(fileName.c_str());

    if (fichierB)
	{
	int compteB = 0;
	for (unsigned long long int pos = 0; pos < nbTotalPoints; pos++)
	    {
	    if (gridPtr[pos] <= level)
		{
		numToIntAndDoubleCoords(pos, xInt, xReel);

		compteB++;

		for (int l1 = 0; l1 < dim; l1++)
		    {
		    fichierB << xInt[l1] << " ";
		    }
		fichierB << gridPtr[pos] << "\n";
		}
	    }

	fichierB.close();
	}
    else
	{
	spdlog::error("[GridMicorMacro] : Error while open file {}", fileName);
	}
    }

void GridMicroMacro::computeMinMaxValues(double &minV, double &maxV) const
    {
    minV = PLUS_INF;
    maxV = -PLUS_INF;

    for (unsigned long long int pos = 0; pos < nbTotalPoints; pos++)
	{
	maxV = max(maxV, gridPtr[pos]);
	minV = min(maxV, gridPtr[pos]);
	}
    }

void GridMicroMacro::saveProjection(const string &fileName, const unsigned long long int *projection) const
    {
    unsigned long long int *reducedNbPoints = new unsigned long long int[dim - 1];
    double *reducedLimInf = new double[dim - 1];
    double *reducedStep = new double[dim - 1];
    unsigned long long int *tempVect = new unsigned long long int[dim - 1];

    int i = 0;
    int j = 0;
    int projAxe = 0;
    while (j < dim)
	{
	if (!projection[j])
	    {
	    reducedNbPoints[i] = nbPoints[j];
	    reducedLimInf[i] = limInf[j];
	    reducedStep[i] = step[j];
	    i++;
	    j++;
	    }
	else
	    {
	    projAxe = j;
	    j++;
	    }
	}

    unsigned long long int reducedTotalPoints = 1;
    for (i = 0; i < dim - 1; i++)
	{
	reducedTotalPoints *= reducedNbPoints[i];
	}
    double val;
    unsigned long long int *x = new unsigned long long int[dim];
    ofstream fichierB(fileName.c_str());

    if (fichierB)
	{
	unsigned long long int pointNum;
	int compteB = 0;
	for (unsigned long long int posX = 0; posX < reducedTotalPoints; posX++)
	    {

	    numToIntCoords_gen(posX, dim - 1, reducedNbPoints, tempVect);
	    for (j = 0; j < projAxe; j++)
		{
		x[j] = tempVect[j];
		}
	    for (j = projAxe + 1; j < dim; j++)
		{
		x[j] = tempVect[j - 1];
		}

	    val = PLUS_INF;
	    for (unsigned long long int k = 0; k < nbPoints[projAxe]; k++)
		{
		x[projAxe] = k;
		pointNum = this->intCoordsToNum(x);
		val = min(val, gridPtr[pointNum]);
		if (val < PLUS_INF - 1)
		    {
		    compteB++;

		    for (int l1 = 0; l1 < dim - 1; l1++)
			{
			fichierB << tempVect[l1] << " ";
			}
		    fichierB << val << "\n";
		    }
		}

	    fichierB.close();
	    }
	}
    else
	{
	spdlog::error("[GridMicorMacro] : Error while open file {}", fileName);
	}
    }

void GridMicroMacro::addPointToSet(const unsigned long long int *coords, double val)
    {

    unsigned long long int pos;

    pos = Grid::intCoordsToNum(coords);
    gridPtr[pos] = val;
    }

void GridMicroMacro::addPointToSet(unsigned long long int pos, double val)
    {
    gridPtr[pos] = val;
    }
