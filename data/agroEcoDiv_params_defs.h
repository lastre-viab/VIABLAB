/*
 * gaia_params_defs.h
 *
 *  Created on: 3 ao�t 2015
 *      Author: anna_desilles
 */

#ifndef GAIA_PARAMS_DEFS_H_
#define GAIA_PARAMS_DEFS_H_


#include "../include/FarmDataGardener.h"

FarmData_Gardener * modelData;
unsigned long long int ***IB, **Delta, **dpr, **perRec;
double  ***bilan, ****bilan_ty,   **sensiR,**ct, **ac, **cIn,**cr,**cie, **rendM, ****maxDeficit, *minMaxDeficit, **sp, **sa, **pv;
string prefix;
double globalMaxDeficit;
double IStar=0.5;
double Imin=0.0;
double Imax=1.0;
string* specNames;
string* pratNames;
double alphaJachere, betaJachere;
 int nbValPoints;
  double initValue;
 double  minVal;
		    double maxVal;
		    double stepVal;
int ***climatSensi;
double **cycloneTauxPerte;
double **coutCyclone;

double ** coeffsIncertPrix;

unsigned long long int cyclonMinPeriode=1;
unsigned long long int cyclonMaxPeriode=1;
unsigned long long int cyclonMaxNb=1;
string ** specPratFileNames;
double * IRealValues;
double deficitMax=500.0;
unsigned long long int nbSaisons;
int excludeSomeControls;
int ** controlExclusions;
double * realIvalues;
double minRichesse;

int nbTrajs;


int nbSpec=13;
int nbPr=2;

unsigned long long int nbIb=11;
unsigned long long int DeltaMax=18;
unsigned long long int nbLignes;

int computeViabSousJacent;
int computeViabTrajSousJacent;
int computeBestPerfSansConstr;
int computeOptimTrajSansConstr;
int computeViabConstrEco;
int computeBestPerfWithConstr;
int computeViabTrajWithConstr;
int computeOptimTrajWithConstr;

double richesseCible;
int useControlPrefs;
int computeTrajectory=1;

int computeWithControlsSubset;
int computeWithControlsPrefs;
int trajReconstructionType=0;
int computeParcelsKernels=0;
int computeAgregRetroDyn=0;
int optimCostType=0;//0= W1 seul (implique W2=0), 1= W2 seul (implique W2=1), 2= W1 sous contrainte de W2
int W2; // 1= crit�re W2, valaeur epigrahique, 0=crit�re W1, valeur hypographique
int computeW2andW1=0;

int *** saisonalite;
#endif /* GAIA_PARAMS_DEFS_H_ */
