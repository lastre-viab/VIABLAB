/*
 * ControlGrid.cpp
 *
 *  Created on: 18 août 2025
 *      Author: adesi
 */

#include "../include/ControlGrid.h"

ControlGrid::ControlGrid()
    {
    dimC = 0;
    stepC = new double[1];
    limSupC = new double[1];
    limInfC = new double[1];
    nbPointsC = new unsigned long long int[1];
    totalNbPointsC = 1;
    controlCoords = new double*[1];
    controlIntCoords = new unsigned long long int*[1];
    }

ControlGrid::~ControlGrid()
    {
    if(totalNbPointsC > 0)
	{
	for (unsigned long long int k = 0; k < totalNbPointsC; k++)

	    {
	    delete[] controlCoords[k];
	    delete[] controlIntCoords[k];
	    }
	}
    delete[] controlCoords;
    delete[] controlIntCoords;
    delete[] stepC;
    delete[] limInfC;
    delete[] limSupC;
    delete[] nbPointsC;
    }
/*!
 *  attention! Pour l'instant, la convention par d�faut
 *  pour la grille de contr�les  est la suivante:
 *  les  points sont les centres de mailles
 *  avec cette convention on peut m�me avoit un point  d controle par axe
 *  sans problemes
 *  Donc le nombre de points=nombre d'intervalles !!!!!
 */
ControlGrid::ControlGrid(int dim, double * limInf, double * limSup, unsigned long long int * nbPoints)
    {
    dimC = dim;
	int tabsize = max(dim, 1);
    stepC = new double[tabsize];
    limSupC = new double[tabsize];
    limInfC = new double[tabsize];
    nbPointsC = new unsigned long long int[tabsize];
if (dimC > 0)
{
	for(int k = 0; k < dim; k++)
	{
		limSupC[k] = limSup[k];
		limInfC[k] = limInf[k];
		nbPointsC[k] = nbPoints[k];
	}

	spdlog::debug("[System] : Control dimension {}", dimC);
	totalNbPointsC = 1;
	for (int dc = 0; dc < dimC; dc++)
	{
		totalNbPointsC *= nbPointsC[dc];
		stepC[dc] = (limSupC[dc] - limInfC[dc]) / (nbPointsC[dc] - 1);
	}

	spdlog::debug("[System] : nbPoints controle total  {}", totalNbPointsC);
	logVector("[System]: control discretization step :", stepC, dimC);
	logVector("[System]: control inf limits :", limInfC, dimC);
	logVector("[System]: control sup limits :", limSupC, dimC);
	controlCoords = new double*[totalNbPointsC];
	controlIntCoords = new unsigned long long int*[totalNbPointsC];
	unsigned long long int *coordsIntC = new unsigned long long int[dimC];

	for (unsigned long long int k = 0; k < totalNbPointsC; k++)
	{
		controlCoords[k] = new double[dimC];
		controlIntCoords[k] = new unsigned long long int[dimC];
		numToIntCoords_gen(k, dimC, nbPointsC, coordsIntC);
		for (int dc = 0; dc < dimC; dc++)
		{
			controlCoords[k][dc] = limInfC[dc] + stepC[dc] * coordsIntC[dc]; //+0.5*stepC[dc];
			controlIntCoords[k][dc] = coordsIntC[dc];
		}
	}
}
	else
	{

			limSupC[0] = 0.0;
			limInfC[0] = 0.0;
			nbPointsC[0] = 1;

		spdlog::debug("[System] : Control dimension {}", dimC);
		totalNbPointsC = 1;
			stepC[0] = 1.0;

		spdlog::debug("[System] : nbPoints controle total  {}", totalNbPointsC);
		logVector("[System]: control discretization step :", stepC, dimC);
		logVector("[System]: control inf limits :", limInfC, dimC);
		logVector("[System]: control sup limits :", limSupC, dimC);
		controlCoords = new double*[totalNbPointsC];
		controlIntCoords = new unsigned long long int*[totalNbPointsC];
		unsigned long long int *coordsIntC = new unsigned long long int[1];

		for (unsigned long long int k = 0; k < totalNbPointsC; k++)
		{
			controlCoords[k] = new double[1];
			controlIntCoords[k] = new unsigned long long int[1];
				controlCoords[k][0] = 0.0; //+0.5*stepC[dc];
				controlIntCoords[k][0] = 0;
		}
	}
    }

double ** ControlGrid::GetControlCoords()
    {
    return controlCoords;
    }

unsigned long long int **ControlGrid::GetControlIntCoords()
    {
    return controlIntCoords;
    }

unsigned long long int ControlGrid::GetTotalNbPoints()
    {
    return totalNbPointsC;
    }

unsigned long long int * ControlGrid::GetNbPoints()
    {
    return nbPointsC;
    }

double * ControlGrid::GetLimSup()
    {
    return limSupC;
    }
double * ControlGrid::GetLimInf()
    {
    return limInfC;
    }
double * ControlGrid::GetStep()
    {
    return stepC;
    }
int ControlGrid::GetDim()
    {
    return dimC;
    }

