/*
 * GaiaModelManager.h
 *
 *  Created on: 6 ao�t 2021
 *      Author: adesi
 */

#ifndef GAIAMODELMANAGER_H_
#define GAIAMODELMANAGER_H_

#include "defs.h"
#include "CycleData.h"

class FarmData_Gardener {
public:
	FarmData_Gardener();
	FarmData_Gardener( string usedCropsFile);

	void readAllCropsFile( string usedCropsFile);
	void readFarmFile( string usedCropsFile);
	void readParamsFile( string usedCropsFile);

	unsigned long long int  getItransition ( unsigned long long int In, unsigned long long int t, unsigned long long int un);
	unsigned long long int  getItransitionPrevious ( unsigned long long int In, unsigned long long int t, unsigned long long int un);

	unsigned long long int getTimeTransition (  unsigned long long int t, unsigned long long int un);
	unsigned long long int getCycleDuration (  unsigned long long int un);
	double getCycleBalance(unsigned long long int In, unsigned long long int t, unsigned long long int un);
	double getCycleMaxDeficit(unsigned long long int In, unsigned long long int t, unsigned long long int un);
	double getCycleSeasonality( unsigned long long int un, unsigned long long int t);
    double getAdmissibleDeficit();
    double getCycleBalanceForGivenPeriod(unsigned long long int In, unsigned long long int mStart, unsigned long long int mEnd, unsigned long long int un);
    double getCycleMaxDeficitForGivenPeriod(unsigned long long int In, unsigned long long int mStart, unsigned long long int mEnd, unsigned long long int un);
    double getCycleMaxDeficitForGivenPeriodMultiParcel(unsigned long long int *In, unsigned long long int *mStart, unsigned long long int nbMonths, unsigned long long int *un);

     void saveFarmDataForMatlab();
	virtual ~FarmData_Gardener();

	unsigned long long int getMaxCycleDuration();
	double * getIgrid();
	unsigned long long int getNbControls();
	unsigned long long int getNbIPoints();
	double getIStar();
	unsigned long long int getTimeHorizon();

	int getNbTrajs();
	unsigned long long int * getInitPointsCoords();
	double * getInitValues();
	double * getParcelsSurfaces();
	int  getNbParcels();

 string getPrefix();
private :
	map<string, CropData*> allCrops;
	CycleData ** cropsPortfolio;
	unsigned long long int ** IBQSTransistionMatrix;
	unsigned long long int ** IBQSTransistionMatrixPrevious;
	unsigned long long int *  cycleDurations;
	double ** cycleCompleteBalance;
	double ** cycleCompleteMaxDeficit;
	int ** seasonalities;
	double * Igrid;
	string prefix;


	int nbAllSpecs;
	int nbCycles;
	unsigned long long int nbPointsIBQS;
	unsigned long long int nbPointsValeur;

	double IBQS_min;
	double IBQS_max;
	double IBQS_target;
	double *surfaces;
	int nbParcels;
	double deficitMax;
	double fixedCostPerHa;
	unsigned long long int maxCycleDuration;
	unsigned long long int timeHorizon;

	int nbTrajectoires;
	unsigned long long int * initPointsCoords;
	double * initValues;







};

#endif /* GAIAMODELMANAGER_H_ */
