/*
 * ViabiMicroMacroTrajectoryHelper.h
 *
 *  Created on: 7 juil. 2024
 *      Author: adesi
 */

#ifndef SRC_VIABIMICROMACROTRAJECTORYHELPER_H_
#define SRC_VIABIMICROMACROTRAJECTORYHELPER_H_
#include "GridMicroMacro.h"
#include "SysDyn.h"

class ViabiMicroMacroTrajectoryHelper {
public:
	ViabiMicroMacroTrajectoryHelper();

	ViabiMicroMacroTrajectoryHelper(GridMicroMacro * gr, SysDyn * ds, int type);

	virtual ~ViabiMicroMacroTrajectoryHelper();

	double computeOptimalTrajectory_tmin(double *initPosition, string fileName, bool &succes);
	double computeOptimalTrajectory_new(double *initPosition, string fileName, bool &succes);

	int findOptiControl_tmin(double *currentPos, unsigned long long int optimDiscreteSuccessor,
			double &dt,
			int nbStepIter,
			double stepCoeff,
			double *resPos,
			bool &succes );

	unsigned long long int  findOptimalDiscreteSuccessor_tmin(unsigned long long int pos, double dt);


	double computeOptimalTrajectory_Lmin(double *initPosition, string fileName, bool &succes);

	int findOptiControl_Lmin(double budget, double *currentPos,
			double &dt,
			int nbStepIter,
			double stepCoeff,
			double *resPos,
			double & newBudget,
			bool &succes );

	double computeViableTrajectory(double  *initPosition,  double initValue, string fileName, bool &succes);
	double computeViableTrajectory_DD(unsigned long long int *initPosition,  double initValue, string fileName, bool &succes);
	double computeViableTrajectory_tych_DD(unsigned long long int *initPosition,  double initValue, string fileName, bool &succes);

	unsigned long long int (ViabiMicroMacroTrajectoryHelper::*findfViableControl_DD)(double budget, unsigned long long int *currentPos,
			unsigned long long int currentControl,
			unsigned long long int *resPos,
			double & newBudget,
			bool &succes );
	unsigned long long int  findViabControlDiffControl_DD(double budget, unsigned long long int *currentPos,
			unsigned long long int currentControl,
			unsigned long long int *resPos,
			double & newBudget,
			bool &succes );
	unsigned long long int  findViabControlDefault_DD(double budget, unsigned long long int *currentPos,
			unsigned long long int currentControl,
			unsigned long long int *resPos,
			double & newBudget,
			bool &succes );
	unsigned long long int  findViabControlDefault_tych_DD(double budget, unsigned long long int *currentPos,
			unsigned long long int currentControl, unsigned long long int currentTych,
			unsigned long long int *resPos,
			double & newBudget,
			bool &succes );
	unsigned long long int  findViabControlMinValue_DD(double budget, unsigned long long int *currentPos,
			unsigned long long int currentControl,
			unsigned long long int *resPos,
			double & newBudget,
			bool &succes );

private :

	GridMicroMacro * grid;
	/*!
	 *  \brief  Copie pour raisons de rapidté  de la valeur de dimension d'état
	 */
	int dim;
	/*!
	 *  \brief  Copie pour raisons de rapidté  de la valeur de dimension de contrôle
	 */
	int dimC;
	/*!
	 *  \brief Pointeur sur la base de données servant à enregister la rétroaction optimale
	 */
	double * vTab;
	int typeTraj;
	string filePrefix;
	SysDyn* dynsys;
};

#endif /* SRC_VIABIMICROMACROTRAJECTORYHELPER_H_ */
