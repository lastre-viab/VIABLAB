/*
 * CycleData.cpp
 *
 *  Created on: 7 août 2021
 *      Author: adesi
 */

#include "../include//CycleData.h"

CycleData::CycleData() {
	// TODO Auto-generated constructor stub

}

CycleData::CycleData( int type, string name,  string longCrop, int nb, string * shorts, double w, double lc, double ll, map<string, CropData*> * allCrops) {

	cycleType = type;
	cycleName = name;
	longCropKey = longCrop;
	shortCropKeys = shorts;
	nbShortCrops = nb;
	weightLong = w;
	LERLong = ll;
	LERShort = lc;
	allCropsData = allCrops;
	shortCropDurations = new unsigned long long int [nbShortCrops];
	shortCropsData = new CropData*[nbShortCrops];
	longCropData = (*allCrops)[longCropKey];
	for(int k=0;k<nbShortCrops; k++)
	{
		shortCropsData[k]=(*allCrops)[shortCropKeys[k]];
		shortCropDurations[k]=shortCropsData[k]->duration;
	}
	longCropDuration = longCropData->duration;
	if(cycleType == 1)
	{
		unsigned long long int tempDuration = 0;
		int i=0;
		while(tempDuration < longCropDuration)
		{
			while(i<nbShortCrops && tempDuration < longCropDuration)
			{
				tempDuration+=shortCropsData[i]->duration;
				i++;
			}
			i=0;
		}
		globalCycleDuration = tempDuration;
	}
	else
	{
		globalCycleDuration = longCropDuration;
	}
	this->fixedCost = 0.0;

}

void CycleData::computeCyclePartialData(string prefix, double * Igrid, unsigned long long int nbIPoints, unsigned long long int T)
{
	bool print = false;//this->cycleType == 1;
	unsigned long long int iNextIndex;
	cyclePartialMaxDeficit = new double *[nbIPoints];
	cyclePartialBalance = new double *[nbIPoints];
	cycleBalancePerMonth = new double *[nbIPoints];
	cyclePartialIBQSTransitionMatrix = new unsigned long long int *[nbIPoints];
	cyclePartialIBQSTransitionMatrixPrevious = new unsigned long long int *[nbIPoints];
	for(int k=0;k<nbIPoints;k++)
		{
			cyclePartialIBQSTransitionMatrixPrevious[k]=new unsigned long long int [globalCycleDuration-1];
			for(int t=0; t<globalCycleDuration-1;t++)
			{
				cyclePartialIBQSTransitionMatrixPrevious[k][t]=k;
			}
		}

	for(int k=0;k<nbIPoints;k++)
	{
		cyclePartialMaxDeficit[k]=new double [globalCycleDuration-1];
		cyclePartialBalance[k]=new double [globalCycleDuration-1];
		cycleBalancePerMonth[k]=new double [globalCycleDuration];
		cyclePartialIBQSTransitionMatrix[k]=new unsigned long long int [globalCycleDuration-1];
		for(int t=0; t<globalCycleDuration-1;t++)
		{
			cyclePartialMaxDeficit[k][t]=this->calculMaxDeficitCycle(Igrid[k], T-globalCycleDuration+1+t, T);
			cyclePartialBalance[k][t]=this->calculBilanCycle(Igrid[k], T-globalCycleDuration+1+t, T);
			cycleBalancePerMonth[k][t] = this->calculBilanMoisCycle(t+1, Igrid[k], this->fixedCost);
					double tempI = this->computeNextIBQSValue(Igrid[k], T-globalCycleDuration+1+t, T);
			if(tempI<Igrid[0])
			{
				cyclePartialIBQSTransitionMatrix[k][t]=2*nbIPoints;
			}
			else
			{
				iNextIndex = (unsigned long long int)max(0,(int)floor((tempI-Igrid[0])/(Igrid[1]-Igrid[0])));
				cyclePartialIBQSTransitionMatrix[k][t]=iNextIndex;
				cyclePartialIBQSTransitionMatrixPrevious[iNextIndex][t]=k;
			}

		}
		cycleBalancePerMonth[k][globalCycleDuration-1] = calculBilanMoisCycle(globalCycleDuration, Igrid[k], this->fixedCost);
	}
	saveCyclePartialData(prefix, nbIPoints, Igrid);
}


double CycleData::calculBilanMoisCycle(unsigned long long int m, double initI, double fixedCost)
{
	double res = 0.0;
	if(cycleType==0)
	{
		res = longCropData->calculBilanMois(m, initI, fixedCost);
	}
	else
	{
		unsigned long long int time=0;
		int numCycle=1;
		double revenuCumul=0.0;

		double I=initI;
		// cout<< " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		// cout<< "Calcul bilan  du cycle " << cycleName <<  " I initial " << initI << endl;

		unsigned long long int currentLongM=1;
		unsigned long long int duree=m;

		int shortCropIndex =0;
		double longRm = longCropData->meanYield;
		double bilanTemp = 0.0;
		unsigned long long int longCropDuration = longCropData->duration;
		double InextLong = longCropData->calculInextAssocie(I,longRm*LERLong, 0,2*longCropDuration);
		double longCropIStep = (InextLong - I)/longCropDuration;
		while(currentLongM<=duree)
		{
			//cout << " currentLongM = "<< currentLongM << " duree "<< duree << endl;
			//cout << " shortCropIndex = "<< shortCropIndex << "  nb shortCrop  "<< this->nbShortCrops<< endl;

			while((shortCropIndex < this->nbShortCrops)  && (currentLongM<=duree))
			{

				CropData * currentShortCrop = this->shortCropsData[shortCropIndex];
				double shortCropDuration = currentShortCrop->duration;
				unsigned long long int currentShortM=1;
				while((currentShortM <= shortCropDuration) && (currentLongM<=duree))
				{
					//cout << " currentLongM = "<< currentLongM << " duree "<< duree << endl;
					//cout << " currentShortM = "<< currentShortM << " nb shortCrop  "<< shortCropsData[shortCropIndex]->name << endl;
					bilanTemp = this->calculBilanMoisMixte(initI, I, currentLongM, currentShortM, currentShortCrop);
					//cout<< " \n le bilan du mois " << bilanTemp << " revenu cumule "<<revenuCumul<< endl;
					currentLongM++;
					currentShortM++;
				}
				double ILongEstimation = initI+min(currentLongM-1, longCropDuration)*longCropIStep;//estimation du IBQS intermédiaire de la spec longue
				double shortRm = currentShortCrop->meanYield;
				double InextShort = currentShortCrop->calculInextAssocie(I, shortRm*LERShort, 0,2*longCropDuration);

				//cout<< " ----------------------------\n";
				//cout<< "ILongEstimation = "<<ILongEstimation<< " Inext short "<< InextShort<< endl;

				//cout<< " Coefs LER :  COURT "<< LERShort <<   " LONG "<< LERLong << endl;
				//cout<< " ----------------------------\n";
				I= min(1.0, max(0.0, LERLong*ILongEstimation+LERShort*InextShort) );
				 //cout<< " nouveau Icourt "<< I <<endl;
				shortCropIndex++;
			}
			shortCropIndex =0;

		}
		res = bilanTemp;
					//cout<< " ICI res est "<<res << " long M = "<< currentLongM << endl ;
	}
	 //cout<< " fini  res = "<<res<<endl;
	return res;
}
double CycleData::calculBilanCycleSingleCrop(double I,int t, int T)
{
	double bl=0.0;
	unsigned long long int duree;
	if(t+globalCycleDuration<=T)
	{
		duree=this->globalCycleDuration;
	}
	else
	{
		duree=T-t;
	}
	for(unsigned long long int  k=1;k<=duree;k++)
	{
		bl=bl+longCropData->calculBilanMois(k, I, fixedCost);
		// cout<< " mois "<<k<< " bilan cumulÃ© "<<bl<<endl;
	}
	//cout<< "=======================\n";
	return bl;
}

double CycleData::calculMaxDeficitSingle(double I,int t, int T)
{
	double df=0.0;
	double dfMax=-PLUS_INF;
	//cout<< " deficit max  de la periode "<<d<<endl;
	unsigned long long int duree;
	if(t+globalCycleDuration<=T)
	{
		duree=this->globalCycleDuration;
	}
	else
	{
		duree=T-t;
	}
	for(unsigned long long int  k=1;k<=duree;k++)
	{
		df=df-longCropData->calculBilanMois(k, I, fixedCost);
		if(df>dfMax)
		{
			dfMax=df;
		}
		//cout<< "  dfMax = "<<dfMax<<endl;
	}

	// cout<< "  dfMax = "<<dfMax<<endl;
	return dfMax;
}


double CycleData::calculBilanCycleMixedCrop(double initI,int t, int T)
{

	unsigned long long int time=0;
	int numCycle=1;
	double revenuCumul=0.0;

	double I=initI;
	//cout<< " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	//cout<< "Calcul bilan  du cycle " << cycleName <<  " I initial " << initI << endl;

	unsigned long long int currentLongM=1;
	unsigned long long int duree;
	if(t+globalCycleDuration<=T)
	{
		duree=this->globalCycleDuration;
	}
	else
	{
		duree=T-t;
	}
	int shortCropIndex =0;
	double longRm = longCropData->meanYield;
	double InextLong = longCropData->calculInextAssocie(I,longRm*LERLong, t,T);
	unsigned long long int longCropDuration = longCropData->duration;
	double longCropIStep = (InextLong - I)/longCropDuration;
	while(currentLongM<=duree)
	{
		//cout << " currentLongM = "<< currentLongM << " duree "<< duree << endl;
		while((shortCropIndex < this->nbShortCrops)  && (currentLongM<=duree))
		{

			CropData * currentShortCrop = this->shortCropsData[shortCropIndex];
			double shortCropDuration = currentShortCrop->duration;
			unsigned long long int currentShortM=1;
			while((currentShortM <= shortCropDuration) && (currentLongM<=duree))
			{
				//cout << " currentLongM = "<< currentLongM << " duree "<< duree << endl;
				//cout << " currentShortM = "<< currentShortM << " shortCrop  "<< shortCropsData[shortCropIndex]->name << endl;
				double bilanTemp = this->calculBilanMoisMixte(initI, I, currentLongM, currentShortM, currentShortCrop);
				revenuCumul+= bilanTemp;
				//cout<< " \n le bilan du mois " << bilanTemp << " revenu cumule "<<revenuCumul<< endl;
				currentLongM++;
				currentShortM++;
			}
			double ILongEstimation = initI+min(currentLongM-1, longCropDuration)*longCropIStep;//estimation du IBQS intermédiaire de la spec longue
			double shortRm = currentShortCrop->meanYield;
			double InextShort = currentShortCrop->calculInextAssocie(I, shortRm*LERShort, t+currentLongM-1,T);

			//cout<< " ----------------------------\n";
			//cout<< "ILongEstimation = "<<ILongEstimation<< " Inext short "<< InextShort<< endl;

			//cout<< " Coefs LER :  COURT "<< LERShort <<   " LONG "<< LERLong << endl;
			//cout<< " ----------------------------\n";
			I= min(1.0, max(0.0, LERLong*ILongEstimation+LERShort*InextShort) );
			//cout<< " nouveau Icourt "<< I <<endl;
			shortCropIndex++;
		}
		shortCropIndex =0;

	}
	//cout<< " dernier Icourt "<< I << " revenu total du cycle "<< revenuCumul<<endl;
	return revenuCumul;
	//cout<< " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}


double CycleData::calculMaxDeficitMixte(double initI,int t, int T)
{

	unsigned long long int time=0;
	int numCycle=1;
	double revenuCumul=0.0;

	double I=initI;

	double df=0.0;
	double dfMax=-PLUS_INF;

	unsigned long long int currentLongM=1;
	unsigned long long int duree;
	if(t+globalCycleDuration<=T)
	{
		duree=this->globalCycleDuration;
	}
	else
	{
		duree=T-t;
	}
	int shortCropIndex =0;
	double longRm = longCropData->meanYield;
	double InextLong = longCropData->calculInextAssocie(I,longRm*LERLong, t,T);
	unsigned long long int longCropDuration = longCropData->duration;
	double longCropIStep = (InextLong - I)/longCropDuration;
	while(currentLongM<=duree)
	{

		while((shortCropIndex < this->nbShortCrops)  && (currentLongM<=duree))
		{
			CropData * currentShortCrop = this->shortCropsData[shortCropIndex];
			double shortCropDuration = currentShortCrop->duration;
			unsigned long long int currentShortM=1;
			while((currentShortM <= shortCropDuration) && (currentLongM<=duree))
			{
				df=df-this->calculBilanMoisMixte(initI, I, currentLongM, currentShortM, currentShortCrop);
				if(df>dfMax)
				{
					dfMax=df;
				}
				currentLongM++;
				currentShortM++;
			}
			double ILongEstimation = initI+min(currentLongM-1, longCropDuration)*longCropIStep;//estimation du IBQS intermédiaire de la spec longue
			double shortRm = currentShortCrop->meanYield;
			double InextShort = currentShortCrop->calculInextAssocie(I, shortRm*LERShort, t+currentLongM-1,T);
			I=LERLong*ILongEstimation+LERShort*InextShort;

			shortCropIndex++;
		}
		shortCropIndex =0;

	}
	return dfMax;
}

double CycleData::computeNextIBQSValue(double I,int t, int T)
{
	return (cycleType==0) ? longCropData->calculInext(I, t, T) : this->calculINextMixedCrop(I, t, T);
}

double CycleData::calculINextMixedCrop(double initI, int t, int T)
{
	//cout<< "Calcul du Inext mixed for "<< cycleName<<endl;
	double I=initI;
	unsigned long long int currentLongM=1;
	unsigned long long int duree;
	if(t+globalCycleDuration<=T)
	{
		duree=this->globalCycleDuration;
	}
	else
	{
		duree=T-t;
		//		cout<< " t = "<<t<<" T= "<< T <<  " Duree = "<< duree << endl;
	}
	int shortCropIndex =0;
	double longRm = longCropData->meanYield;
	double InextLong = longCropData->calculInextAssocie(I, longRm*LERLong, t,T);
	unsigned long long int longCropDuration = longCropData->duration;
	double longCropIStep = (InextLong - I)/longCropDuration;
	//cout<< " Inext longue "<< InextLong << "  step " << longCropIStep<<endl;
	while(currentLongM<=duree)
	{

		while((shortCropIndex < this->nbShortCrops)  && (currentLongM<=duree))
		{
			CropData * currentShortCrop = this->shortCropsData[shortCropIndex];
			double shortCropDuration = currentShortCrop->duration;
			unsigned long long int currentShortM=1;
			while((currentShortM <= shortCropDuration) && (currentLongM<=duree))
			{
				currentLongM++;
				currentShortM++;
			}

			// cout<<" fin de cycle court "<< shortCropsData[shortCropIndex]->name<< " mois long "<<currentLongM<< " mois court "<< currentShortM<< endl;
			double ILongEstimation = initI+min(currentLongM-1, longCropDuration)*longCropIStep;//estimation du IBQS intermédiaire de la spec longue
			double shortRm = currentShortCrop->meanYield;
			double InextShort = currentShortCrop->calculInextAssocie(I,shortRm*LERShort, t+currentLongM-1, T);
			I=min(1.0, LERLong*ILongEstimation+LERShort*InextShort);
			//cout<< " Inext longue "<< ILongEstimation << "  I next short  " << InextShort<< "  nouveau I "<<I<<endl;
			shortCropIndex++;
		}
		shortCropIndex =0;
	}
	//system("pause");
	return I;
}


double CycleData::calculBilanMoisMixte(double Ilong, double Ishort, unsigned long long int mLong, unsigned long long int mShort, CropData * shortCrop)
{
	//cout<< " fixedCost = "<<fixedCost<<endl;
	double bilanL=  longCropData->calculBilanMoisAssocie(mLong, Ilong, fixedCost, weightLong, LERLong);
	double bilanC= shortCrop->calculBilanMoisAssocie(mShort, Ishort, fixedCost, 1.0-weightLong, LERShort);
	//cout<< " bilan mois mixte mois court :  "<< bilanC << " mois long "<< bilanL<< endl;
	return bilanL+bilanC;

}


double CycleData::calculBilanCycle(double I,int t, int T)
{
	return (cycleType==0)? this->calculBilanCycleSingleCrop(I, t, T):this->calculBilanCycleMixedCrop(I, t, T);
}
double CycleData::calculMaxDeficitCycle(double I,int t, int T)
{
	return (cycleType==0)? this->calculMaxDeficitSingle(I, t, T):this->calculMaxDeficitMixte(I, t, T);
}


unsigned long long int CycleData::GetDuration()
{
	return globalCycleDuration;
}

void CycleData::setFixedCost(double fc)
{
	fixedCost = fc;
}


unsigned long long int  CycleData::getIPartialtransition ( unsigned long long int In, unsigned long long int t)
{
	return cyclePartialIBQSTransitionMatrix[In][t];
}

unsigned long long int  CycleData::getIPartialtransitionPrevious ( unsigned long long int In, unsigned long long int t)
{
	return cyclePartialIBQSTransitionMatrixPrevious[In][t];
}

double CycleData::getCyclePartialBalance(unsigned long long int In, unsigned long long int t)
{
	return cyclePartialBalance[In][t];
}

double CycleData::getCycleBalancePerMonth(unsigned long long int In, unsigned long long int m)
{
	//cout<< " getting cycle month for In = "<<In<< " month = "<<m<<endl;
	return cycleBalancePerMonth[In][m-1];
}
double CycleData::getCyclePartialMaxDeficit(unsigned long long int In, unsigned long long int t)
{
	return  cyclePartialMaxDeficit[In][t];
}

int CycleData::getCycleSeasonality(int t)
{
	return longCropData->admissibleSeason[t];
}

void CycleData::saveCyclePartialData(string prefix, unsigned long long int nbIPoints, double * Igrid)
{
	ostringstream os;
	string fileName;
	FILE * fi;

	os<<"../OUTPUT/"<<prefix<<"-"<<cycleName<<"-BilansCycle.dat";
	fileName=os.str();
	os.str("");

	fi = fopen( fileName.c_str(),"w");
	if(fi==NULL){
		printf("** error: impossible to open the file %s.\n", fileName.c_str());

	}
	else
	{
		for(int t=0; t<globalCycleDuration;t++)
		{
			fprintf(fi,  "%d " ,    t);
		}
		fprintf(fi,  "\n");
		for(int l1=0;l1<nbIPoints;l1++)
		{
			double bilanComplet = calculBilanCycle(Igrid[l1], 0, globalCycleDuration+100);

			fprintf(fi,  "%15.8f " ,    bilanComplet);
			for(int t=0; t<globalCycleDuration-1;t++)
			{
				fprintf(fi,  "%15.8f " ,    cyclePartialBalance[l1][t]);
			}

			fprintf(fi,  "\n");
		}

		fclose(fi);
	}

	os<<"../OUTPUT/"<<prefix<<"-"<<cycleName<<"-BilansCyclePerMonth.dat";
		fileName=os.str();
		os.str("");

		fi = fopen( fileName.c_str(),"w");
		if(fi==NULL){
			printf("** error: impossible to open the file %s.\n", fileName.c_str());

		}
		else
		{
			for(int t=0; t<globalCycleDuration;t++)
			{
				fprintf(fi,  "%d " ,    t);
			}
			fprintf(fi,  "\n");
			for(int l1=0;l1<nbIPoints;l1++)
			{

				for(int t=0; t<globalCycleDuration;t++)
				{
					fprintf(fi,  "%15.8f " ,    cycleBalancePerMonth[l1][t]);
				}

				fprintf(fi,  "\n");
			}

			fclose(fi);
		}

	os<<"../OUTPUT/"<<prefix<<"-"<<cycleName<<"-PartialIBQS.dat";
	fileName=os.str();
	os.str("");

	fi = fopen( fileName.c_str(),"w");
	if(fi==NULL){
		printf("** error: impossible to open the file %s.\n", fileName.c_str());

	}
	else
	{
		for(int t=0; t<globalCycleDuration;t++)
		{
			fprintf(fi,  "%d " ,    t);
		}
		fprintf(fi,  "\n");
		for(int l1=0;l1<nbIPoints;l1++)
		{
			double IbqsComplet = computeNextIBQSValue(Igrid[l1], 0, globalCycleDuration+100);
			unsigned long long int index =(unsigned long long int)max(0,(int)floor((IbqsComplet-Igrid[0])/(Igrid[1]-Igrid[0])));
			fprintf(fi,  "%d " ,    index);
			for(int t=0; t<globalCycleDuration-1;t++)
			{
				fprintf(fi,  "%d " ,    cyclePartialIBQSTransitionMatrix[l1][t]);
			}

			fprintf(fi,  "\n");
		}

		fclose(fi);
	}

	os<<"../OUTPUT/"<<prefix<<"-"<<cycleName<<"-PartialIBQSPrevious.dat";
		fileName=os.str();
		os.str("");

		fi = fopen( fileName.c_str(),"w");
		if(fi==NULL){
			printf("** error: impossible to open the file %s.\n", fileName.c_str());

		}
		else
		{
			for(int t=0; t<globalCycleDuration-1;t++)
			{
				fprintf(fi,  "%d " ,    t);
			}
			fprintf(fi,  "\n");
			for(int l1=0;l1<nbIPoints;l1++)
			{
				//double IbqsComplet = computeNextIBQSValue(Igrid[l1], 0, globalCycleDuration+100);
				//unsigned long long int index =(unsigned long long int)max(0,(int)floor((IbqsComplet-Igrid[0])/(Igrid[1]-Igrid[0])));
				//fprintf(fi,  "%d " ,    index);
				for(int t=0; t<globalCycleDuration-1;t++)
				{
					fprintf(fi,  "%d " ,    cyclePartialIBQSTransitionMatrixPrevious[l1][t]);
				}

				fprintf(fi,  "\n");
			}

			fclose(fi);
		}

}

CycleData::~CycleData() {
	// TODO Auto-generated destructor stub
}

