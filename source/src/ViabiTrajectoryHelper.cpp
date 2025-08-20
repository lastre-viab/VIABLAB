/*
 * ViabiTrajectoryHelper.cpp
 *
 *  Created on: 19 août 2025
 *      Author: adesi
 */

#include "../include/ViabiTrajectoryHelper.h"

ViabiTrajectoryHelper::ViabiTrajectoryHelper()
    {
    // TODO Auto-generated constructor stub

    }

ViabiTrajectoryHelper::~ViabiTrajectoryHelper()
    {
    // TODO Auto-generated destructor stub
    }

ViabiTrajectoryHelper::ViabiTrajectoryHelper(SysDyn *ds, TrajectoryParametersManager *tpm) :
    tpm(tpm),
    dynsys(ds) {
    dim = dynsys->getDim();
    dimC = dynsys->getDimC();

    const trajectoryParams *tp = tpm->getTrajectoryParameters();
        TypeTraj typeTraj = tp->TRAJECTORY_TYPE;

        unsigned long long int nbTotalC = dynsys->getTotalNbPointsC();
        preferedControlIndexes = new int[nbTotalC];
        std::iota(preferedControlIndexes, preferedControlIndexes + nbTotalC, 0);
        controlWeight = tp->CONTROL_WEIGHT;
        sortIndexes = tp->SORT_INDEXES;
        trajIndex = tpm->getTrajectoryIndex();
}

bool ViabiTrajectoryHelper::isViableControl(const double *currentPos, const double *controlCoord, double *imageVect, double rho) const {
    if (dynsys->constraintsXU(currentPos, controlCoord) >= PLUS_INF) {
	return false;
    }
    else {
	std::invoke(dynsys->discretDynamics, dynsys, currentPos, controlCoord, imageVect, rho);
	return checkKernelRelation(imageVect) == VALID_TRAJECTORY_POINT;
    }
}

bool ViabiTrajectoryHelper::isViableControl_tych(const double *currentPos, const double *controlCoord, const double *tycheCoord, double *imageVect, double rho) const {
    if (dynsys->constraintsXU(currentPos, controlCoord) >= PLUS_INF || dynsys->constraintsXV_tych(currentPos, tycheCoord) >= PLUS_INF) {
	return false;
    }
    else {
	std::invoke(dynsys->discretDynamics_tych, dynsys, currentPos, controlCoord, tycheCoord, imageVect, rho);
	return checkKernelRelation(imageVect) == VALID_TRAJECTORY_POINT;
    }
}

bool ViabiTrajectoryHelper::isViableGuaranteedControl(const double *currentPos, const double *controlCoord, double rho) const {

    if (dynsys->constraintsXU(currentPos, controlCoord) >= PLUS_INF) {
	return false;
    }

    bool allTychOk = true;
    unsigned long long int i = 0;
    double *imageVect = new double[dim];
    unsigned long long int totalNbPointsC = dynsys->getTotalNbPointsC();
    unsigned long long int totalNbPointsTych = dynsys->getTotalNbPointsTy();

    double ** controlCoords = dynsys->getControlCoords();
    double ** tychCoords = dynsys->getTychCoords();
    while (allTychOk && i < totalNbPointsTych) {
	const double *tychCoord = tychCoords[i];
	if (dynsys->constraintsXV_tych(currentPos, tychCoord) < PLUS_INF) {
	    std::invoke(dynsys->discretDynamics_tych, dynsys, currentPos, controlCoord, tychCoord, imageVect, rho);
	    allTychOk = (checkKernelRelation(imageVect) == VALID_TRAJECTORY_POINT);
	}
	++i;
    }
    delete [] imageVect;
    return allTychOk;
}

int ViabiTrajectoryHelper::getClosestControlTo(const double *u)
    {
    double **controlCoords = dynsys->getControlCoords();
    unsigned long long nbCTotal = dynsys->getTotalNbPointsC();

    double minDist = PLUS_INF;
    int argMinDist = -1;

    for (unsigned long long i = 0; i < nbCTotal; ++i)
	{
	double euclideanDistanceSquared = 0;
	double *uCandidate = controlCoords[i];
	for (int d = 0; d < dimC; ++d)
	    {
	    euclideanDistanceSquared += (u[d] - uCandidate[d]) * (u[d] - uCandidate[d]);
	    }
	if (euclideanDistanceSquared < minDist)
	    {
	    minDist = euclideanDistanceSquared;
	    argMinDist = i;
	    }
	}
    return argMinDist;
    }

SysDyn * ViabiTrajectoryHelper::GetDynSys()
    {
    return dynsys;
    }
