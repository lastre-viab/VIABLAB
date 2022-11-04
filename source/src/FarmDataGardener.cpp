/*
 * GaiaModelManager.cpp
 *
 *  Created on: 6 ao�t 2021
 *      Author: adesi
 */

#include "../include/FarmDataGardener.h"

FarmData_Gardener::FarmData_Gardener() {
	// TODO Auto-generated constructor stub

}

void FarmData_Gardener::readFarmFile( string farmSimulationFile)
{
	string line;
	ostringstream os;

	string tempStr, tmpStr;
	ptree dataRoot;
	string sedCommand;
	string input_tempfile;
	string cropName, practName, cycleName;


	input_tempfile="../EXEC_DATA/"+ farmSimulationFile;
	//sedCommand="sed  -e 's/[ ]*#.*//g' -e '/^$/d' ../EXEC_DATA/"+ farmSimulationFile+" > "+input_tempfile;
	//unsigned long long int  zz=system(sedCommand.c_str()); zz++;

	//- chargement des donnees depuis le fichier json via le parseur boost
	read_json(input_tempfile, dataRoot);

	prefix=dataRoot.get<string>("NAME", "test");
	nbParcels = dataRoot.get<int>("NB_PARCELLES", 1);
	surfaces = new double[nbParcels];
	for(unsigned long long int tabIndice=0;tabIndice<nbParcels;tabIndice++)
	{
		surfaces[tabIndice]=1.0;//valeurs par default
	}
	if(dataRoot.find("SURFACES")!=dataRoot.not_found())
	{
		ptree  tabTree=dataRoot.find("SURFACES")->second;
		ptree::const_iterator it = tabTree.begin();
		for(unsigned long long int tabIndice=0;tabIndice<nbParcels;tabIndice++)
		{
			surfaces[tabIndice]=(*it).second.get_value<double>();
			it++;

			cout<< " surface = "<<surfaces[tabIndice]<<endl;
		}
	}
	deficitMax = dataRoot.get<double>("DEFICITE_MAX", 0.0);
	fixedCostPerHa = dataRoot.get<double>("COUTFIXE_PAR_HA", 0.0);
	if(dataRoot.find("CYCLES_LIST")!=dataRoot.not_found())
	{

		ptree  cyclesTree=dataRoot.find("CYCLES_LIST")->second;
		this->nbCycles=cyclesTree.count("CYCLE");
		cout<< " nb cycles is "<<nbCycles<<endl;

		cropsPortfolio = new CycleData *[nbCycles];
		unsigned long long int cycleNum = 0;
		for (ptree::value_type &cy : cyclesTree.get_child(""))
		{
			ptree cycleRoot = cy.second;
			cycleName=cycleRoot.get<string>("NAME", "DefaultCycle");
			string cycleType = cycleRoot.get<string>("TYPE", "single");
			int type =0;
			if(cycleType == "mixte")
			{
				type = 1;
			}
			string pratique = cycleRoot.get<string>("PRATIQUE", "Eco");
			string longCrop = cycleRoot.get<string>("SPEC_LONGUE", "DefaultCrop")+"_"+pratique;
			cout<< " spec longue "<<longCrop<<endl;
			int nbShortCrops = cycleRoot.get<int>("NB_SPEC_COURTES", 0);
			double w = cycleRoot.get<double>("POIDS_SPEC_LONGUE", 1.0);
			double lc = cycleRoot.get<double>("LER_COURTE", 1-w);
			double ll = cycleRoot.get<double>("LER_LONGUE", w);
			string * shortCrops = new string[nbShortCrops];
			if(cycleRoot.find("SPECS_COURTES")!=cycleRoot.not_found())
			{
				ptree  tabTree=cycleRoot.find("SPECS_COURTES")->second;
				ptree::const_iterator it = tabTree.begin();
				for(unsigned long long int tabIndice=0;tabIndice<nbShortCrops;tabIndice++)
				{
					shortCrops[tabIndice]=(*it).second.get_value<string>()+"_"+pratique;;
					it++;

					cout<< " short crop= "<<shortCrops[tabIndice]<<endl;
				}
			}

			cropsPortfolio[cycleNum] = new CycleData(type, cycleName, longCrop, nbShortCrops, shortCrops,w, lc, ll, &allCrops);
			cout<< " Cycle num "<<cycleNum<< " "<< cropsPortfolio[cycleNum]->cycleName<< " duree "<< cropsPortfolio[cycleNum]->GetDuration()<<endl;
			cycleNum++;

		}
	}

	if(dataRoot.find("ANNEX_SURFACE")!=dataRoot.not_found())
	{
		annexCoutFixe = dataRoot.get<double>("COUTFIXE_PAR_HA", 0.0);

		farmWithAnnexSurface = true;
		ptree  annexRoot = dataRoot.find("ANNEX_SURFACE")->second;
		annexSurface = annexRoot.get<double>("SURFACE", 1.0);
		annexInitIbqs = annexRoot.get<double>("INITIAL_IBQS", 1.0);
		ptree annexCycleRoot = annexRoot.find("CYCLE")->second;;
		cycleName=annexCycleRoot.get<string>("NAME", "DefaultCycle");
		string cycleType = annexCycleRoot.get<string>("TYPE", "single");
		int type =0;
		if(cycleType == "mixte")
		{
			type = 1;
		}
		string pratique = annexCycleRoot.get<string>("PRATIQUE", "Eco");
		string longCrop = annexCycleRoot.get<string>("SPEC_LONGUE", "DefaultCrop")+"_"+pratique;
		cout<< "  ANNEXE spec longue "<<longCrop<<endl;
		int nbShortCrops = annexCycleRoot.get<int>("NB_SPEC_COURTES", 0);
		double w = annexCycleRoot.get<double>("POIDS_SPEC_LONGUE", 1.0);
		double lc = annexCycleRoot.get<double>("LER_COURTE", 1-w);
		double ll = annexCycleRoot.get<double>("LER_LONGUE", w);
		string * shortCrops = new string[nbShortCrops];
		if(annexCycleRoot.find("SPECS_COURTES")!=annexCycleRoot.not_found())
		{
			ptree  tabTree=annexCycleRoot.find("SPECS_COURTES")->second;
			ptree::const_iterator it = tabTree.begin();
			for(unsigned long long int tabIndice=0;tabIndice<nbShortCrops;tabIndice++)
			{
				shortCrops[tabIndice]=(*it).second.get_value<string>()+"_"+pratique;;
				it++;

				cout<< " short crop= "<<shortCrops[tabIndice]<<endl;
			}
		}

		annexCycle = new CycleData(type, cycleName, longCrop, nbShortCrops, shortCrops,w, lc, ll, &allCrops);
	}
	else
	{
		farmWithAnnexSurface = false;
		annexSurface = 0.0;
		annexInitIbqs = 0.0;
	}
	//system("pause");
}
void FarmData_Gardener::readAllCropsFile( string allCropsListFile)
{

	string line;
	ostringstream os;

	string tempStr, tmpStr;
	ptree dataRoot;
	string sedCommand;
	string input_tempfile;
	string cropName, practName, cycleName;


	//elimination des commentaires dans les fichiers d'entrée
	input_tempfile="../EXEC_DATA/"+ allCropsListFile;
	//sedCommand="sed  -e 's/[ ]*#.*//g' -e '/^$/d' ../EXEC_DATA/"+ allCropsListFile+" > "+input_tempfile;
	//unsigned long long int zz=system(sedCommand.c_str()); zz++;

	//- chargement des donnees depuis le fichier json via le parseur boost
	read_json(input_tempfile, dataRoot);

	allCrops.clear();

	for (ptree::value_type &crop : dataRoot.get_child("CROPS_LIST"))
	{
		// Get the content of the node
		std::string cropName = crop.second.data();
		allCrops[cropName+"_Eco"]=new CropData(cropName+"_Eco.json");
		allCrops[cropName+"_Conv"]=new CropData(cropName+"_Conv.json");
	}


	this->nbAllSpecs = allCrops.size();
}
void FarmData_Gardener::readParamsFile( string numParamsFile){
	FILE * pFile;
	string line;


	ostringstream os;

	string tempStr, tmpStr;
	ptree dataRoot;
	string sedCommand;
	string input_tempfile;
	string cropName, practName, cycleName;

	input_tempfile="../EXEC_DATA/"+ numParamsFile;
	//sedCommand="sed  -e 's/[ ]*#.*//g' -e '/^$/d' ../EXEC_DATA/"+ numParamsFile+" > "+input_tempfile;
	//unsigned long long int zz=system(sedCommand.c_str()); zz++;

	//- chargement des donnees depuis le fichier json via le parseur boost
	read_json(input_tempfile, dataRoot);

	IBQS_min = dataRoot.get<double>("IBQS_MIN", 0.0);
	IBQS_max = dataRoot.get<double>("IBQS_MAX", 1.0);
	IBQS_target = dataRoot.get<double>("IBQS_TARGET", 0.9);
	this->nbPointsIBQS = dataRoot.get<int>("NB_POINTS_IBQS", 50);
	this->nbPointsValeur = dataRoot.get<int>("NB_POINTS_VALEUR", 50);

	timeHorizon = dataRoot.get<int>("HORIZON", 120);
	IBQS_level_first_parcel = dataRoot.get<double>("IBQS_LEVEL_FIRST_PARCEL", 0.0);

	ptree trajsRoot = dataRoot.find("TRAJECTOIRES")->second;
	nbTrajectoires=trajsRoot.count("INIT_POINT");
	cout<< " nb trajs is "<<nbTrajectoires<<endl;
	initPointsCoords = new unsigned long long int [nbTrajectoires*2];
	initValues = new double[nbTrajectoires];
	int iValue = 0;
	int iCoords = 0;
	double dI=(IBQS_max-IBQS_min)/(nbPointsIBQS-1);
	for (ptree::value_type &p : trajsRoot.get_child(""))
	{
		ptree pointRoot = p.second;
		double ibqsReel  = pointRoot.get<double>("IBQS", 0.1);
		initPointsCoords[iCoords] = min(this->nbPointsIBQS - 1, max((unsigned long long int )0,(unsigned long long int ) round( (ibqsReel - IBQS_min)/dI )));
		//initPointsCoords[iCoords] = min(this->nbPointsIBQS - 1, max((unsigned long long int )0,(unsigned long long int ) floor( (ibqsReel - IBQS_min)/dI )));
		iCoords++;
		initPointsCoords[iCoords] = (unsigned long long int) pointRoot.get<int>("DATE_DEBUT", 0);
		iCoords++;
		initValues[iValue]  = pointRoot.get<double>("BUDGET", 1.0);
		iValue++;

	}

	cout<< "lecture de trajectoires fini init points : "<<endl;
	for(int k = 0; k<nbTrajectoires*2; k++)
	{
		cout<< " "<< initPointsCoords[k];
	}
	cout<<endl;


}

FarmData_Gardener::FarmData_Gardener(string sourceFile)
{
	// TODO Auto-generated constructor stub
	FILE * pFile;
	string line;

	string allCropsListFile = "allCropsList.json";
	string farmDescriptionFile = "exploitation.json";
	string numParamsFile = "exploitation.json";



	ostringstream os;

	string tempStr, tmpStr;
	ptree dataRoot;
	string sedCommand;
	string input_tempfile;
	string cropName, practName, cycleName;

	input_tempfile="../EXEC_DATA/"+ sourceFile;
	//sedCommand="sed  -e 's/[ ]*#.*//g' -e '/^$/d' ../EXEC_DATA/"+ sourceFile+" > "+input_tempfile;
	//unsigned long long int zz=system(sedCommand.c_str()); zz++;

	//- chargement des donnees depuis le fichier json via le parseur boost
	read_json(input_tempfile, dataRoot);

	allCropsListFile = dataRoot.get<string>("SPECLISTE", "LesSpeculations.json");
	this->readAllCropsFile(allCropsListFile);

	farmDescriptionFile = dataRoot.get<string>("EXPLOITATION", "Exploitation.json");
	this->readFarmFile(farmDescriptionFile);

	numParamsFile = dataRoot.get<string>("NUM_PARAMS", "NumParams.json");

	this->readParamsFile(numParamsFile);

	Igrid=new double[nbPointsIBQS];
	double dI=(IBQS_max-IBQS_min)/(nbPointsIBQS-1);
	for(int ii=0;ii<(int)nbPointsIBQS;ii++)
	{
		Igrid[ii]=IBQS_min+ii*dI;
		//	cout<< " I real = "<<IRealValues[ii]<<endl;
	}

	cout<< " computing partial data nb cycles "<<nbCycles<<endl;
	for(int k=0;k<nbCycles;k++)
	{
		cout<< " Cycle "<<cropsPortfolio[k]->cycleName<<endl;
		cropsPortfolio[k]->setFixedCost(this->fixedCostPerHa);
		cropsPortfolio[k]->computeCyclePartialData(prefix, Igrid, nbPointsIBQS, timeHorizon);
	}

	cycleDurations = new unsigned long long int[nbCycles];
	seasonalities =  new int *[nbCycles];
	maxCycleDuration=0;
	for(int c=0;c<nbCycles;c++)
	{
		cycleDurations[c] = cropsPortfolio[c]->GetDuration();
		maxCycleDuration = max(maxCycleDuration, cycleDurations[c]);
		seasonalities[c] = new int[12];
		for(int s=0;s<12;s++)
		{
			seasonalities[c][s] = cropsPortfolio[c]->getCycleSeasonality(s);
		}
	}

	IBQSTransistionMatrix = new unsigned long long int*[nbPointsIBQS];
	IBQSTransistionMatrixPrevious = new unsigned long long int*[nbPointsIBQS];
	cycleCompleteBalance = new double * [nbPointsIBQS];
	cycleCompleteMaxDeficit = new double * [nbPointsIBQS];
	for(int i=0;i<nbPointsIBQS;i++)
	{
		IBQSTransistionMatrixPrevious[i] = new unsigned long long int[nbCycles];
		for(int c=0;c<nbCycles;c++)
		{
			IBQSTransistionMatrixPrevious[i][c] = i;
		}
	}
	unsigned long long int iNextIndex;
	for(int i=0;i<nbPointsIBQS;i++)
	{
		IBQSTransistionMatrix[i] = new unsigned long long int[nbCycles];
		cycleCompleteBalance[i] = new double[nbCycles];
		cycleCompleteMaxDeficit[i] = new double[nbCycles];

		for(int c=0;c<nbCycles;c++)
		{
			double tempI = cropsPortfolio[c]->computeNextIBQSValue(Igrid[i], 0, timeHorizon);
			//cout<< " cycle "<<cropsPortfolio[c]->cycleName<< " ibqs current"<< Igrid[i]<< " ibqs next "<< tempI<<endl;
			if(tempI<Igrid[0])
			{
				IBQSTransistionMatrix[i][c]=2*nbPointsIBQS;
			}
			else
			{
				iNextIndex = (unsigned long long int)min((int)nbPointsIBQS - 1, max(0, (int)round((tempI-Igrid[0])/(Igrid[1]-Igrid[0]))));
				//iNextIndex = (unsigned long long int)min((int)nbPointsIBQS - 1, max(0, (int)floor((tempI-Igrid[0])/(Igrid[1]-Igrid[0]))));
				IBQSTransistionMatrix[i][c] = iNextIndex;
				IBQSTransistionMatrixPrevious[iNextIndex][c] = i;
			}
			cycleCompleteBalance[i][c] = cropsPortfolio[c]->calculBilanCycle(Igrid[i], 0, timeHorizon);
			cycleCompleteMaxDeficit[i][c] = cropsPortfolio[c]->calculMaxDeficitCycle(Igrid[i], 0, timeHorizon);
		}
	}
	bool print =false;
	if(print)
	{
		for(int c=0;c<nbCycles;c++)
		{
			cout<<"====================================================================================================\n";

			cout<< " ====================== Partial data of the cycle "<< cropsPortfolio[c]->cycleName<<endl;
			cout<< " cycleCompleteMaxDeficit "<<endl;
			for(int k=0;k<nbPointsIBQS;k++)
			{

				cout<< " "<<cycleCompleteMaxDeficit[k][c];
				cout<<endl;
			}
			cout<<endl;
			cout<< " cycleCompleteBalance "<<endl;
			for(int k=0;k<nbPointsIBQS;k++)
			{
				cout<< " "<<cycleCompleteBalance[k][c];
				cout<<endl;
			}
			cout<<endl;
			cout<< " IBQSTransistionMatrix "<<endl;
			for(int k=0;k<nbPointsIBQS;k++)
			{
				cout<< " "<<IBQSTransistionMatrix[k][c];
				cout<<endl;
			}
			cout<<endl;
			cout<<"====================================================================================================\n";
		}
	}
	annexCycleBalance = new double [timeHorizon];
	if(farmWithAnnexSurface)
	{
		unsigned long long int duree = annexCycle->GetDuration();
		unsigned long long t = 0;
		double ibqs = annexInitIbqs;
		unsigned long long int  currentMonth = 1;
		while(t<timeHorizon)
		{
			currentMonth = 1;
			while ( currentMonth <= duree && t< timeHorizon)
			{
				annexCycleBalance[t] = annexSurface * annexCycle->calculBilanMoisCycle(currentMonth, ibqs, annexCoutFixe);
				t++;
				currentMonth++;
			}
			ibqs = annexCycle->computeNextIBQSValue(ibqs, t, timeHorizon);

		}
	}
	else
	{
		for(int t = 0;t<timeHorizon;t++)
		{
			annexCycleBalance[t] = 0.0;
		}
	}
	saveFarmDataForMatlab();
	//exit(1);
}

unsigned long long int FarmData_Gardener::getMaxCycleDuration()
{
	return maxCycleDuration;
}
double* FarmData_Gardener::getIgrid()
{
	return Igrid;
}
unsigned long long int FarmData_Gardener::getNbControls()
{
	return nbCycles;
}
unsigned long long int FarmData_Gardener::getNbIPoints()
{
	return nbPointsIBQS;
}
double FarmData_Gardener::getIStar()
{
	return IBQS_target;
}
unsigned long long int FarmData_Gardener::getTimeHorizon()
{
	return timeHorizon;
}
double FarmData_Gardener::getIbqsLevelForSavings()
{
	return IBQS_level_first_parcel;
}

double FarmData_Gardener::getCycleSeasonality( unsigned long long int un, unsigned long long int t)
{
	return seasonalities[un][t];
}

double FarmData_Gardener::getAdmissibleDeficit()
{
	return deficitMax;
}

unsigned long long int  FarmData_Gardener::getItransition ( unsigned long long int In, unsigned long long int t, unsigned long long int un)
{
	if(t+cycleDurations[un]<=timeHorizon)
	{
		//cout<< " ibq = "<<Igrid[In]<< " matrice transi "<<IBQSTransistionMatrix[In][un];
		return IBQSTransistionMatrix[In][un];

	}
	else
	{
		return cropsPortfolio[un]->getIPartialtransition(In, t- (timeHorizon-cycleDurations[un]+1));
	}
}

unsigned long long int  FarmData_Gardener::getItransitionPrevious ( unsigned long long int In, unsigned long long int t, unsigned long long int un)
{

	if(t+cycleDurations[un]<=timeHorizon)
	{
		// cout<< " ibq = "<<Igrid[In]<< " matrice transi previous  "<<IBQSTransistionMatrixPrevious[In][un]<<endl;
		return IBQSTransistionMatrixPrevious[In][un];

	}
	else
	{
		//cout<< "For PARTIAL PREVIOUS t = "<< t<< " timeHorizon " << timeHorizon   << "  cycleDurations[un] "<< cycleDurations[un]<<endl;
		// cout<< "For PARTIAL PREVIOUS t- (timeHorizon-cycleDurations[un]) = "<< t- (timeHorizon-cycleDurations[un])<<endl;
		return cropsPortfolio[un]->getIPartialtransitionPrevious(In, t- (timeHorizon-cycleDurations[un]+1));
	}
}

unsigned long long int FarmData_Gardener::getTimeTransition (  unsigned long long int t, unsigned long long int un)
{
	//cout<< " cycle duration is "<<cycleDurations[un]<<endl;
	return min(cycleDurations[un], timeHorizon - t);
}

unsigned long long int FarmData_Gardener::getCycleDuration (unsigned long long int un)
{
	//cout<< " cycle duration is "<<cycleDurations[un]<<endl;
	return cycleDurations[un];
}

double FarmData_Gardener::getCycleBalance(unsigned long long int In, unsigned long long int t, unsigned long long int un)
{
	if(t+cycleDurations[un]<=timeHorizon)
	{
		return cycleCompleteBalance[In][un];
	}
	else
	{
		return cropsPortfolio[un]->getCyclePartialBalance(In, t- (timeHorizon-cycleDurations[un]+1));
	}
}

double FarmData_Gardener::getCyclePerturbedBalance(unsigned long long int In, unsigned long long int t, unsigned long long int un, double v)
{
	if(t+cycleDurations[un]<=timeHorizon)
	{
		return cycleCompleteBalance[In][un]+cropsPortfolio[un]->calculDeltaBilanCycle(Igrid[In], t, timeHorizon, v);
	}
	else
	{
		return cropsPortfolio[un]->getCyclePartialBalance(In, t- (timeHorizon-cycleDurations[un]+1));
	}
}

double FarmData_Gardener::getCycleBalanceForGivenPeriod(unsigned long long int In, unsigned long long int mStart, unsigned long long int mEnd, unsigned long long int un)
{
	double res = 0.0;

	for(unsigned long long int m = mStart ; m<= mEnd; m++)
	{


		res+= cropsPortfolio[un]->getCycleBalancePerMonth(In, m);

	}
	return res;
}

double FarmData_Gardener::getCycleBalanceWithPricePerturbForGivenPeriod(unsigned long long int In,
		unsigned long long int mStart, unsigned long long int mEnd, unsigned long long int un, double v)
{
	double res = 0.0;

	for(unsigned long long int m = mStart ; m<= mEnd; m++)
	{


		res+= cropsPortfolio[un]->getCycleBalancePerMonth(In, m) + cropsPortfolio[un]->calculDeltaBilanMoisCycle(m, In, v);

	}
	return res;
}

double FarmData_Gardener::getMinCycleBalanceForGivenPeriod(unsigned long long int In, unsigned long long int gamma, unsigned long long int nbMonths)
{
	double res = PLUS_INF;
	for(int k=0;k<nbCycles;k++)
	{
		unsigned long long int delta2 = cropsPortfolio[k]->GetDuration();
		if(delta2>gamma)
		{
			unsigned long long int mStart2 = delta2 - gamma +1;
			unsigned long long int mEnd2 = mStart2+nbMonths-1;
			unsigned long long int previousI = getItransitionPrevious(In, 0, k);
			double bilan2 = getCycleBalanceForGivenPeriod(previousI, mStart2, mEnd2, k);
			res = min(res, bilan2);
		}
	}
	return res;
}

double FarmData_Gardener::getCycleMaxDeficitForGivenPeriod(unsigned long long int In, unsigned long long int mStart, unsigned long long int mEnd, unsigned long long int un)
{
	double df=0.0;
	double dfMax=-PLUS_INF;
	//cout<<"DEFICIT MAX time period is "<<mStart<< " "<<mEnd<<endl;
	for(unsigned long long int m = mStart ; m<= mEnd; m++)
	{
		//cout<< " m = "<<m<< " mEnd = " <<mEnd<<endl;
		df=df-cropsPortfolio[un]->getCycleBalancePerMonth(In, m);
		if(df>dfMax)
		{
			dfMax=df;
		}
	}
	//cout<<" res = "<<dfMax<<endl;
	return dfMax;
}

double FarmData_Gardener::getCycleMaxDeficitForGivenPeriodMultiParcel(unsigned long long int *In, unsigned long long int *mStart, unsigned long long int nbMonths, unsigned long long int *un)
{
	double df=0.0;
	double dfMax=-PLUS_INF;
	//cout<<"DEFICIT MAX time period is "<<mStart<< " "<<mEnd<<endl;
	for(unsigned long long int m = 0 ; m< nbMonths; m++)
	{
		//cout<< " m = "<<m<< " mEnd = " <<mEnd<<endl;

		for(int k=0;k<nbParcels; k++)
		{
			df=df-surfaces[k]*cropsPortfolio[un[k]]->getCycleBalancePerMonth(In[k], mStart[k]+m);
		}

		if(df>dfMax)
		{
			dfMax=df;
		}
	}

	//cout<<" res = "<<dfMax<<endl;
	return dfMax;
}
double FarmData_Gardener::getCycleMaxDeficit(unsigned long long int In, unsigned long long int t, unsigned long long int un)
{

	if(t+cycleDurations[un]<=timeHorizon)
	{
		//cout<< " max def du u = "<<un<< " I "<<In<< " t = "<<t<<endl;
		//cout<< " complete cycle max def "<<cycleCompleteMaxDeficit[In][un]<<endl;
		return cycleCompleteMaxDeficit[In][un];
	}
	else
	{
		//cout<< " max def du u = "<<un<< " I "<<In<< " t = "<<t<<endl;
		// cout<< " partial max def pour tau = "<< t- (timeHorizon-cycleDurations[un])<<" "<<cropsPortfolio[un]->getCyclePartialMaxDeficit(In, t- (timeHorizon-cycleDurations[un]))<<endl;
		return cropsPortfolio[un]->getCyclePartialMaxDeficit(In, t- (timeHorizon-cycleDurations[un]+1));
	}
}

double FarmData_Gardener::getCycleMaxDeficitWithAnnexe(unsigned long long int In, unsigned long long int t, unsigned long long int un)
{

	double df=0.0;
	double dfMax=-PLUS_INF;
	unsigned long long int nbMonths = cycleDurations[un];
	if(t+cycleDurations[un]<=timeHorizon)
	{
		nbMonths = timeHorizon - t;
	}
	//cout<<"DEFICIT MAX time period is "<<mStart<< " "<<mEnd<<endl;
	for(unsigned long long int m = 1 ; m<=nbMonths; m++)
	{
		df=df-(cropsPortfolio[un]->getCycleBalancePerMonth(In, m)*surfaces[0] + annexCycleBalance[t+m-1]);
		if(df>dfMax)
		{
			dfMax=df;
		}
	}
	//cout<<" res = "<<dfMax<<endl;
	return dfMax;
}


int FarmData_Gardener::getNbTrajs()
{
	return nbTrajectoires;
}
unsigned long long int * FarmData_Gardener::getInitPointsCoords()
{
	return initPointsCoords;
}
double * FarmData_Gardener::getInitValues()
{
	return initValues;
}

void FarmData_Gardener::saveFarmDataForMatlab()
{

	ostringstream os;
	string fileName;
	FILE * fi;

	os<<"../OUTPUT/"<<prefix<<"-IBQSGrid.dat";
	fileName=os.str();
	os.str("");

	fi = fopen( fileName.c_str(),"w");
	if(fi==NULL){
		printf("** error: impossible to open the file %s.\n", fileName.c_str());

	}
	else
	{
		for(int l1=0;l1<nbPointsIBQS;l1++)
		{
			fprintf(fi,  "%15.8f\n" ,    Igrid[l1]);
		}

		fclose(fi);
	}

//	os<<"../OUTPUT/"<<prefix<<"-BilansAnnex.dat";
//	fileName=os.str();
//	os.str("");
//
//	fi = fopen( fileName.c_str(),"w");
//	if(fi==NULL){
//		printf("** error: impossible to open the file %s.\n", fileName.c_str());
//
//	}
//	else
//	{
//		for(int l1=0;l1<timeHorizon;l1++)
//		{
//			fprintf(fi,  "%15.8f\n" ,  annexCycleBalance[l1]);
//		}
//
//		fclose(fi);
//	}

	os<<"../OUTPUT/"<<prefix<<"-CycleNames.dat";
	fileName=os.str();
	os.str("");

	fi = fopen( fileName.c_str(),"w");
	if(fi==NULL){
		printf("** error: impossible to open the file %s.\n", fileName.c_str());

	}
	else
	{
		for(int l1=0;l1<nbCycles;l1++)
		{
			fprintf(fi,  "%s\n" ,    cropsPortfolio[l1]->cycleName.c_str());
		}

		fclose(fi);
	}

	os<<"../OUTPUT/"<<prefix<<"-simuData.dat";
	fileName=os.str();
	os.str("");

	fi = fopen( fileName.c_str(),"w");
	if(fi==NULL){
		printf("** error: impossible to open the file %s.\n", fileName.c_str());

	}
	else
	{

		fprintf(fi,  "%d\n" ,    nbCycles);


		fclose(fi);
	}

	os<<"../OUTPUT/simuPrefix.dat";
	fileName=os.str();
	os.str("");

	fi = fopen( fileName.c_str(),"w");
	if(fi==NULL){
		printf("** error: impossible to open the file %s.\n", fileName.c_str());

	}
	else
	{
		fprintf(fi,  "%s\n" , prefix.c_str());

		fclose(fi);
	}
	for(int c=0;c<nbCycles;c++)
	{
		os<<"../OUTPUT/"<<prefix<<"-"<<cropsPortfolio[c]->cycleName<<"-IBQSdyn.dat";
		fileName=os.str();
		os.str("");
		fi = fopen( fileName.c_str(),"w");
		if(fi==NULL){
			printf("** error: impossible to open the file %s.\n", fileName.c_str());

		}
		else
		{
			for(int l1=0;l1<nbPointsIBQS;l1++)
			{
				fprintf(fi,  "%d\n" ,    IBQSTransistionMatrix[l1][c]);
			}

			fclose(fi);
		}
	}
}

double FarmData_Gardener::getAnnexCycleBalance(unsigned long long start, unsigned long long int end)
{
	double result = 0.0;
	for(unsigned long long int k = start; k<=end; k++)
	{
		result += annexCycleBalance[k];
	}
	return result;
}

int  FarmData_Gardener::getNbParcels()
{
	return nbParcels;
}

double * FarmData_Gardener::getParcelsSurfaces()
{
	return surfaces;
}

string FarmData_Gardener::getPrefix()
{
	return prefix;
}



FarmData_Gardener::~FarmData_Gardener() {
	// TODO Auto-generated destructor stub
}

