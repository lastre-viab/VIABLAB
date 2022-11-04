/*
 * CropData.cpp
 *
 *  Created on: 7 ao�t 2021
 *      Author: adesi
 */

#include "../include/CropData.h"
CropData::CropData() {
	// TODO Auto-generated constructor stub

}

CropData::CropData(string fichierSource) {
	// TODO Auto-generated constructor stub
	FILE * pFile;
	string line;


	ostringstream os;

	string tempStr, tmpStr;
	ptree dataRoot;
	string sedCommand;
	string input_tempfile;

	//elimination des commentaires dans les fichiers d'entrée
	input_tempfile="../EXEC_DATA/"+fichierSource;
	//sedCommand="sed  -e 's/[ ]*#.*//g' -e '/^$/d' ../EXEC_DATA/"+fichierSource+" > "+input_tempfile;
	//unsigned long long int zz=system(sedCommand.c_str()); zz++;

	//- chargement des donnees depuis le fichier json via le parseur boost
	read_json(input_tempfile, dataRoot);


	name = dataRoot.get<string>("NAME", "DefaultCrop");
	practice = dataRoot.get<int>("PRATIQUE", 1);

	duration = dataRoot.get<int>("DUREE_CYCLE", 1);
	cout<< "Name de crop "<< name << " duree "<< duration<< endl;
	admissibleSeason = new int[12];
	if(dataRoot.find("SAISONALITE")!=dataRoot.not_found())
	{
		ptree  dataTree=dataRoot.find("SAISONALITE")->second;
		ptree::const_iterator it = dataTree.begin();

		for(unsigned long long int tabIndice=0;tabIndice<12;tabIndice++)
		{
			admissibleSeason[tabIndice]=(*it).second.get_value<int>();
			it++;
		}
	}

	sensi = dataRoot.get<double>("SENSI_IBQS", 0.5);
	meanYield = dataRoot.get<double>("RENDEMENT_MOYEN", 1.0);
	soilDegradationRate = dataRoot.get<double>("TAUX_DEGRADATION_IBQS", 0.5);
	practiceImpactRate = dataRoot.get<double>("IMPACT_PRATIQUE_IBQS", 0.05);

	absolutePracticeImpact = dataRoot.get<double>("IMPACT_PRATIQUE_ABSOLU_IBQS", 0.05);
	activityGrantPerMonth = dataRoot.get<double>("SUBVENTION_ACTIVITE", 0.0);
	productionGrantPerMonth = dataRoot.get<double>("SUBVENTION_PRATIQUE", 0.0);
	workCostPerMonth = dataRoot.get<double>("COUT_TRAVAIL", 0.0);
	startingConst = dataRoot.get<double>("COUT_INSTALLATION", 1.0);
	otherCosts = dataRoot.get<double>("AUTRES_COUTS", 1.0);
	pricePerTonn = dataRoot.get<double>("PRIX_VENTE_PAR_TONNE", 1.0);
	firstCroppingMonth = dataRoot.get<int>("PREMIER_MOIS_RECONLTE", 1);
	croppingPeriod = dataRoot.get<int>("PERIODE_RECOLTE", 1);
	croppingCost = dataRoot.get<double>("COUT_RECOLTE", 1.0);
	inputCosts = dataRoot.get<double>("COUT_INTRANTS", 1.0);
	climateSensi = dataRoot.get<double>("TAUX_DEGRADATION_CYCLONE", 1.0);
	restoreCost = dataRoot.get<double>("COUT_RESTAURATION_CYCLONE", 1.0);

	cout<< " fini de lire les donnees de spec "<<this->name<<endl;
	//system("pause" );
}

double CropData::calculRendement(  double I)
{
	double Rm = meanYield;
	double s = sensi;
	double R=Rm;
	double Rma=Rm*(1.0+s);
	double Rmi=Rm*(1.0-2.0*s);

	if(s<0.5)
		if(I<0.5)
			R=(1./(1.-2.0*s))*(-(2.0*(1.-2.0*s)*I+2.0*s)+sqrt((2.0*s)*(2.0*s)+8.*(1.-2.0*s)*I))*Rm;
	//R=Rmi+(1./(1.-2.0*s))*(-(2.0*(1.-2.0*s)*I+2.0*s)+sqrt((2.0*s)*(2.0*s)+8.*(1.-2.0*s)*I))*(Rm-Rmi);
		else
			R=Rm+(Rma-Rm)* (1./((1.-2.0*s)))*(-(2.0*(1.-2.0*s)*(I-0.5)-(2.-2.0*s))-sqrt((2.-2.0*s)*(2-2.0*s)-8.*(1.0-2.0*s)*(I-0.5)) );

	else
		if(s>0.5)
			if(I<=0.5)

				R=Rm*(1/((2.0*s-1)))*(-(2.0*(2.0*s-1)*(I)-2.0*s)-sqrt((2.0*s)*(2.0*s)-8.*(2.0*s-1)*(I)) );
			else

				R=Rm+(Rma-Rm)*(1/( 2.0*s-1))*(-(2.0*(2.0*s-1)*(I-0.5)+(2-2.0*s))+sqrt((2-2.0*s)*(2-2.0*s)+8.*(2.0*s-1)*(I-0.5)));

		else
			if(I<=0.5)
				R=2.0*I*Rm;
			else
				R=(2.0*I-1)*(Rma-Rm)+Rm;

	return R;
}

double CropData::calculRendementAssocie(  double I, double Rm)
{
	double s = sensi;
	double R=Rm;
	double Rma=Rm*(1.0+s);
	double Rmi=Rm*(1.0-2.0*s);

	if(s<0.5)
		if(I<0.5)
			//R=(1./(1.-2.0*s))*(-(2.0*(1.-2.0*s)*I+2.0*s)+sqrt((2.0*s)*(2.0*s)+8.*(1.-2.0*s)*I))*Rm;
			R=Rmi+(1./(1.-2.0*s))*(-(2.0*(1.-2.0*s)*I+2.0*s)+sqrt((2.0*s)*(2.0*s)+8.*(1.-2.0*s)*I))*(Rm-Rmi);
		else
			R=Rm+(Rma-Rm)* (1./((1.-2.0*s)))*(-(2.0*(1.-2.0*s)*(I-0.5)-(2.-2.0*s))-sqrt((2.-2.0*s)*(2-2.0*s)-8.*(1.0-2.0*s)*(I-0.5)) );

	else
		if(s>0.5)
			if(I<=0.5)

				R=Rm*(1/((2.0*s-1)))*(-(2.0*(2.0*s-1)*(I)-2.0*s)-sqrt((2.0*s)*(2.0*s)-8.*(2.0*s-1)*(I)) );
			else

				R=Rm+(Rma-Rm)*(1/( 2.0*s-1))*(-(2.0*(2.0*s-1)*(I-0.5)+(2-2.0*s))+sqrt((2-2.0*s)*(2-2.0*s)+8.*(2.0*s-1)*(I-0.5)));

		else
			if(I<=0.5)
				R=2.0*I*Rm;
			else
				R=(2.0*I-1)*(Rma-Rm)+Rm;

	return R;
}

double CropData::calculDeltaIbiomasse(  double I)
{
	if(meanYield>0)
	{
		double R;
		R=calculRendement(I);

		return  (0.5* (R/meanYield)*soilDegradationRate);
	}
	else
	{
		return 0.0;
	}
}

double CropData::calculDeltaIbiomasseAssocie(  double I, double Rm)
{
	if(Rm>0)
	{
		double R;
		R=calculRendementAssocie(I, Rm);

		return  (0.5* (R/Rm)*soilDegradationRate);
	}
	else
	{
		return 0.0;
	}
}

double CropData::calculDeltaIpratique(  double I)
{
	return I*practiceImpactRate + absolutePracticeImpact;
}
double CropData::calculInext(double I, unsigned long long  int t, unsigned long long int T)
{
	if(t+duration <=T)
	{
		return min(1.0,  (I-calculDeltaIbiomasse(I)+calculDeltaIpratique(I)));
	}
	else
	{
		double Inext = min(1.0,  (I-calculDeltaIbiomasse(I)+calculDeltaIpratique(I)));
		double coefficient = (double)(T-t)/((double) duration);
		return I+(Inext- I)*coefficient;
	}
}

unsigned long long int CropData::calculINextIndex(double I, double stepI, double IMin, unsigned long long int nbPointsI)
{
	double dIb, dIp, dIj;
	dIb=calculDeltaIbiomasse(I);
	dIp=calculDeltaIpratique(I);

	int niPlus;
	double Iplus=min(1.0,  (I-dIb+dIp));
	if(Iplus<IMin)
	{
		niPlus=2*nbPointsI;
		//cout<< " dynamique : on sort \n";
	}
	else
		niPlus=min((int)nbPointsI - 1, max(0,(int)round((Iplus-IMin)/(stepI))));

	return (unsigned long long int)niPlus;
}


double CropData::calculInextAssocie(double I, double Rm, unsigned long long  int t, unsigned long long int T)
{
	if(t+duration <=T)
	{
		return min(1.0,  (I-calculDeltaIbiomasseAssocie(I, Rm)+calculDeltaIpratique(I)));
	}
	else
	{
		double Inext = min(1.0,  (I-calculDeltaIbiomasseAssocie(I, Rm)+calculDeltaIpratique(I)));
		double coefficient = (double)(T-t)/((double) duration);
		return I+(Inext- I)*coefficient;
	}
}

unsigned long long int CropData::calculINextIndexAssocie(double I, double Rm, double stepI, double IMin, unsigned long long int nbPointsI)
{
	double dIb, dIp;
	dIb=calculDeltaIbiomasseAssocie(I, Rm);
	dIp=calculDeltaIpratique(I);

	int niPlus;
	double Iplus=min(1.0,  (I-dIb+dIp));
	if(Iplus<IMin)
	{
		niPlus=2*nbPointsI;
		//cout<< " dynamique : on sort \n";
	}
	else
		niPlus=min((int)nbPointsI - 1, max(0, (int)round((Iplus-IMin)/(stepI))));

	return (unsigned long long int)niPlus;
}

// calcul du bilan financier /mois et /par Ha pour la sp�c donn�e
// sans les co�ts fixes, g�r�s au niveau de l'exploitation
double CropData::calculBilanMois(unsigned long long int t, double I, double fixedCost)
{

	if(fixedCost < 1.0)
	{
		cout<< " probleme fixed cost = " << fixedCost<< endl;
		exit(1);
	}
	if(t>duration)
	{
		return 0.0;
	}
	double s= sensi;
	double Rm= meanYield;

	double R=calculRendement(I);
	//cout<< " rendement : "<<R;
	double bl=  -this->workCostPerMonth+activityGrantPerMonth-inputCosts-fixedCost;

	if(Rm>0)

	{
		bl=bl+(R*(pricePerTonn+productionGrantPerMonth)-croppingCost*(R/Rm)*(Rm>0))*(t==firstCroppingMonth);
		//cout<< " recolte cost "<<costs[4]*(R/Rm)*(Rm>0)<< " gain recolte "<<R*(p+Sp)<<endl;
		if(t==1)
		{
			bl=bl-startingConst;
		}
		else
		{
			if(t==duration)
			{
				bl=bl-this->otherCosts;
			}
		}
		if((t>firstCroppingMonth) & (croppingPeriod>0))
		{
			if(t%croppingPeriod==0)//mois de r�colte
			{
				bl=bl+R*(pricePerTonn+productionGrantPerMonth)-croppingCost*(R/Rm)*(Rm>0);
			}
		}
	}
	else
	{

		if(t==1)
		{
			bl=bl-startingConst;
		}
		else
		{
			if(t==duration)
			{
				bl=bl-this->otherCosts;
			}
		}

	}
	// cout<< " bilan mois t="<<t<<" bilan "<<bl<<endl;
	return bl;
}

double CropData::calculDeltaBilanMois(unsigned long long int t, double I, double v)
{


	if(t>duration)
	{
		return 0.0;
	}
	double s= sensi;
	double Rm= meanYield;

	double R=calculRendement(I);
	//cout<< " rendement : "<<R;
	double bl=  0.0;

	if(Rm>0)
	{
		bl=bl+(R*(pricePerTonn*v))*(t==firstCroppingMonth);
		//cout<< " recolte cost "<<costs[4]*(R/Rm)*(Rm>0)<< " gain recolte "<<R*(p+Sp)<<endl;

		if((t>firstCroppingMonth) & (croppingPeriod>0))
		{
			if(t%croppingPeriod==0)//mois de r�colte
			{
				bl=bl+R*(pricePerTonn*v);
			}
		}
	}

	// cout<< " bilan mois t="<<t<<" bilan "<<bl<<endl;
	return bl;
}

double CropData::calculBilanMoisAssocie(unsigned long long int t, double I, double fixedCost, double wSurface, double lerPartiel)
{

	if(fixedCost < 1.0)
	{
		cout<< " probleme fixed cost = " << fixedCost<< endl;
		exit(1);
	}
	if(t>duration)
	{
		return 0.0;
	}
	double Rm= meanYield*lerPartiel;

	double R=calculRendementAssocie(I, Rm);
	//cout<< " rendement : "<<R;
	double bl= wSurface*( -this->workCostPerMonth+activityGrantPerMonth-inputCosts-fixedCost);

	if(Rm>0)

	{
		bl=bl+(R*(pricePerTonn+productionGrantPerMonth)-wSurface*croppingCost*(R/Rm))*(t==firstCroppingMonth);
		//cout<< " recolte cost "<<costs[4]*(R/Rm)*(Rm>0)<< " gain recolte "<<R*(p+Sp)<<endl;
		if(t==1)
		{
			bl=bl-wSurface*startingConst;
		}
		else
		{
			if(t==duration)
			{
				bl=bl-wSurface*this->otherCosts;
			}
		}
		if((t>firstCroppingMonth) & (croppingPeriod>0))
		{
			if(t%croppingPeriod==0)//mois de r�colte
			{
				bl=bl+R*(pricePerTonn+productionGrantPerMonth)-wSurface*croppingCost*(R/Rm)*(Rm>0);
			}
		}
	}
	else
	{

		if(t==1)
		{
			bl=bl-wSurface*startingConst;
		}
		else
		{
			if(t==duration)
			{
				bl=bl-wSurface*this->otherCosts;
			}
		}

	}
	// cout<< " bilan mois t="<<t<<" bilan "<<bl<<endl;
	return bl;
}

double CropData::calculDeltaBilanMoisAssocie(unsigned long long int t,
		double I, double fixedCost, double wSurface, double lerPartiel, double v)
{

	if(t>duration)
	{
		return 0.0;
	}
	double Rm= meanYield*lerPartiel;

	double R=calculRendementAssocie(I, Rm);
	//cout<< " rendement : "<<R;
	double bl= 0.0;

	if(Rm>0)

	{
		bl=bl+(R*(pricePerTonn*v))*(t==firstCroppingMonth);
		//cout<< " recolte cost "<<costs[4]*(R/Rm)*(Rm>0)<< " gain recolte "<<R*(p+Sp)<<endl;


		if((t>firstCroppingMonth) & (croppingPeriod>0))
		{
			if(t%croppingPeriod==0)//mois de r�colte
			{
				bl=bl+R*(pricePerTonn*v);
			}
		}
	}

	// cout<< " bilan mois t="<<t<<" bilan "<<bl<<endl;
	return bl;
}

CropData::~CropData() {
	// TODO Auto-generated destructor stub
}

