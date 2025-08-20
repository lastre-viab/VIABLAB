/*
 * GridBitSetHybrid.cpp
 *
 *  Created on: 17 août 2025
 *      Author: adesi
 */

#include "../include/GridBitSetHybrid.h"


GridBitSetHybrid::~GridBitSetHybrid() {
	// TODO Auto-generated destructor stub
}

GridBitSetHybrid::GridBitSetHybrid(const gridParams &gp) : Grid_BitSet(gp)
{
	dim_hc = gp.DIM_HC;
	dim_hd = gp.DIM_HD;
	nbCellPointsContinuousState = pow(2, dim_hc);
	continuousStateShifts = new unsigned long long int[nbCellPointsContinuousState];
	for(int k = 0; k < (int)nbCellPointsContinuousState; k++)
	{
		continuousStateShifts[k] = indicesDecalCell[k];
	}
}

unsigned long long int * GridBitSetHybrid::GetContinuousStateShifts()
{
	return continuousStateShifts;
}

unsigned long long int GridBitSetHybrid::getNbCellsContinuousState()
{
	return nbCellPointsContinuousState;
}


unsigned long long int GridBitSetHybrid::LocalizeHybridPoint(double * doubleCoords, unsigned long long int * intCoords)
{
    unsigned long long int numCell = (unsigned long long int) (((doubleCoords)[0] - limInf[0]) / step[0]);
	if (numCell == (nbPoints[0] - 1)) numCell--;

    for (int i = 1; i < dim_hc; i++)
	{
        unsigned long long int indiceCell = (unsigned long long int) (((doubleCoords)[i] - limInf[i]) / step[i]);
        if (indiceCell == (nbPoints[i] - 1)) {
            indiceCell--;
        }
        numCell = numCell * (nbPoints[i]) + indiceCell;
	}

    for (int i = dim_hc; i < dim; i++)
    	{
            numCell = numCell * (nbPoints[i]) + intCoords[i];
    	}
    return numCell;
}
