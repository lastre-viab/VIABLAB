/*
 * CropData.h
 *
 *  Created on: 7 août 2021
 *      Author: adesi
 */

#ifndef INCLUDE_CROPDATA_H_
#define INCLUDE_CROPDATA_H_

#include "defs.h"
class CropData {
public:
	CropData();
	CropData(string fileMane);
	string name;
	int duration;
	int * admissibleSeason;
	double sensi;
	double meanYield;
	double soilDegradationRate;
	double practiceImpactRate;
	int practice;
	double absolutePracticeImpact;
	double activityGrantPerMonth;
	double productionGrantPerMonth;
	double workCostPerMonth;
	double startingConst;
	double otherCosts;
	double pricePerTonn;
	int firstCroppingMonth;
	int croppingPeriod;
	double croppingCost;
	double inputCosts;
	double climateSensi;
	double restoreCost;
	double calculInext(double I, unsigned long long int t, unsigned long long int T);
	unsigned long long int calculINextIndex(double I, double stepI, double IMin, unsigned long long int nbPointsI);
	double calculInextAssocie(double I, double Rm, unsigned long long int t, unsigned long long int T);
		unsigned long long int calculINextIndexAssocie(double I, double Rm, double stepI, double IMin, unsigned long long int nbPointsI);

	double calculRendementAssocie(  double I, double Rm);
	double calculBilanMois(unsigned long long int t, double I, double fixedCost);
	double calculBilanMoisAssocie(unsigned long long int t, double I, double fixedCost, double wSurface, double lerPartiel);
	virtual ~CropData();

private :


	double calculRendement(  double I);
	double calculDeltaIbiomasse(  double I);
	double calculDeltaIbiomasseAssocie(  double I, double Rm);
	double calculDeltaIpratique(  double I);
	//unsigned long long int calculInext(double I, double stepI, double IMin, unsigned long long int nbPointsI);



};

#endif /* INCLUDE_CROPDATA_H_ */
