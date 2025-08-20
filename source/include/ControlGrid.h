/*
 * ControlGrid.h
 *
 *  Created on: 18 août 2025
 *      Author: adesi
 */

#ifndef INCLUDE_CONTROLGRID_H_
#define INCLUDE_CONTROLGRID_H_

#include "defs.h"

class ControlGrid final
    {
public:
    ControlGrid();
    ControlGrid(int, double *, double *, unsigned long long int *);
    virtual ~ControlGrid();

    ControlGrid(const ControlGrid&) = delete;
    ControlGrid(const ControlGrid&&) = delete;
    ControlGrid& operator=(const ControlGrid&) = delete;
    ControlGrid& operator=(ControlGrid &&data);

    //Getters
    int GetDim();
    double * GetLimInf();
    double * GetLimSup();
    double * GetStep();
    unsigned long long int GetTotalNbPoints();
    unsigned long long int * GetNbPoints();
    double **GetControlCoords();
    unsigned long long int **GetControlIntCoords();

private:
    int dimC;
    double *limInfC;
    double *limSupC;
    double *stepC;
    unsigned long long int *nbPointsC;
    double **controlCoords;
    unsigned long long int **controlIntCoords;
    unsigned long long int totalNbPointsC;
    };

#endif /* INCLUDE_CONTROLGRID_H_ */
