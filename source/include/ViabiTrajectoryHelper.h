/*
 * ViabiTrajectoryHelper.h
 *
 *  Created on: 19 août 2025
 *      Author: adesi
 */

#ifndef INCLUDE_VIABITRAJECTORYHELPER_H_
#define INCLUDE_VIABITRAJECTORYHELPER_H_
#include "SysDyn.h"
#include "ParametersManager.h"

class ViabiTrajectoryHelper
    {
public:
    ViabiTrajectoryHelper();
    ViabiTrajectoryHelper(SysDyn *ds, TrajectoryParametersManager *tpm);
    virtual ~ViabiTrajectoryHelper();
    enum PointStatus : unsigned char
	{
	VALID_TRAJECTORY_POINT,
	OUTSIDE_DOMAIN,
	OUTSIDE_CONSTRAINTS,
	OUTSIDE_GRID,
	};
    virtual PointStatus checkKernelRelation(double *point) const = 0;

    bool isViableControl(const double *x, const double *u, double *image, double rho) const;
    bool isViableControl_tych(const double *x, const double *u, const double *v, double *image, double rho) const;
    bool isViableGuaranteedControl(const double *x, const double *u, double rho) const;
    /**
         * Trouve le contrôle "le plus proche" de u.
         * La distance choisie est ici (arbitrairement) la distance euclidienne.
         */
        int getClosestControlTo(const double *u);
    SysDyn * GetDynSys();
    indexSorter_t sortIndexes;
protected:
    int dim;
    int dimC;

    TrajectoryParametersManager *tpm;
    SysDyn *dynsys;


    int *preferedControlIndexes;
        controlWeight_t controlWeight;
        int trajIndex;
    };

#endif /* INCLUDE_VIABITRAJECTORYHELPER_H_ */
