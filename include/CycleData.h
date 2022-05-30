/*
 * CycleData.h
 *
 *  Created on: 7 août 2021
 *      Author: adesi
 */

#ifndef SRC_CYCLEDATA_H_
#define SRC_CYCLEDATA_H_
#include "defs.h"
#include "CropData.h"
class CycleData {
public:
	CycleData();
	string cycleName;
	int cycleType; // 0 = single, 1= mixed
	string longCropKey;
	string * shortCropKeys;
	int nbShortCrops;
	CycleData( int type, string name, string longCrop, int nbShortCrops, string * shorts, double w, double lc, double ll, map<string, CropData*> * allCrops);
	double weightLong;
	double LERLong;
	double LERShort;
	unsigned long long int GetDuration();
	double calculBilanCycle(double I,int t, int T);
	double calculMaxDeficitCycle(double I,int t, int T);
	double calculINextMixedCrop(double initI,int t, int T);
	double computeNextIBQSValue(double initI,int t, int T);
	void setFixedCost(double fc);
	void computeCyclePartialData(string prefix, double * Igrid, unsigned long long int nbIpoints, unsigned long long int T);
	 double calculBilanMoisCycle(unsigned long long int t, double I, double fixedCost);
	unsigned long long int  getIPartialtransition ( unsigned long long int In, unsigned long long int t);
	unsigned long long int  getIPartialtransitionPrevious ( unsigned long long int In, unsigned long long int t);

	double getCycleBalancePerMonth(unsigned long long int In, unsigned long long int m);
	double getCyclePartialBalance(unsigned long long int In, unsigned long long int t);
	double getCyclePartialMaxDeficit(unsigned long long int In, unsigned long long int t);
	void saveCyclePartialData(string prefix, unsigned long long int nbIPoints, double * Igrid);
int getCycleSeasonality(int t);
	virtual ~CycleData();

private :
	double fixedCost;
	map<string, CropData*> * allCropsData;
	CropData * longCropData;
	CropData ** shortCropsData;
	double ** cyclePartialMaxDeficit;
	double ** cyclePartialBalance;
	double ** cycleBalancePerMonth;
	unsigned long long int ** cyclePartialIBQSTransitionMatrix;
	unsigned long long int ** cyclePartialIBQSTransitionMatrixPrevious;

	unsigned long long int globalCycleDuration;
	unsigned long long int longCropDuration;
	unsigned long long int * shortCropDurations;
	double calculBilanCycleSingleCrop(double I,int t, int T);
	double calculBilanCycleMixedCrop(double I,int t, int T);
	double calculMaxDeficitSingle(double I,int t, int T);
	double calculMaxDeficitMixte(double I,int t, int T);
	double calculBilanMoisMixte(double Ilong, double Ishort, unsigned long long int mLong, unsigned long long int mShort, CropData * shortCrop);


};

#endif /* SRC_CYCLEDATA_H_ */
