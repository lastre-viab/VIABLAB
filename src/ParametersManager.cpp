/*
 * ParametersManager.cpp
 *
 *  Created on: 27 juil. 2021
 *      Author: adesi
 */

#include "../include/ParametersManager.h"


ParametersManager::ParametersManager() {
	// TODO Auto-generated constructor stub

}

ParametersManager::ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams *cp,systemParams *sp)
{
	gridParameters =gp;
	algoParameters=avp;
	controlParameters=cp;
	systemParameters=sp;
}

ParametersManager::ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams * cp,  systemParams *sp, int nbThreads, string paramsFile)
{
	gridParameters =gp;
	algoParameters=avp;
	controlParameters=  cp;
	systemParameters=sp;
	nbOmpThreads = nbThreads;
	parametersFileName = paramsFile;
	cout<< " starting reading data from json file \n";
	this->readControlParametersFromJson();
	this->readGridParametersFromJson();
	this->readAlgoParametersFromJson();
	this->readSystemParametersFromJson();
}

gridParams * ParametersManager::getGridParameters(){
	return gridParameters;
}
algoViabiParams * ParametersManager::getAlgoParameters(){
	return algoParameters;
}
controlParams * ParametersManager::getControlParameters(){
	return controlParameters;
}
systemParams * ParametersManager::getSystemParameters(){
	return systemParameters;
}

void ParametersManager::readGridParametersFromJson()
{

	cout<< "  REad grig params \n";
	string line;
	ostringstream os;

	string tempStr, tmpStr;
	ptree allParamsRoot, dataRoot;
	string sedCommand;

	string input_tempfile="../INPUT/"+ parametersFileName;
	//sedCommand="sed  -e 's/[ ]*#.*//g' -e '/^$/d' ../EXEC_DATA/"+ farmSimulationFile+" > "+input_tempfile;
	//unsigned long long int  zz=system(sedCommand.c_str()); zz++;

	//- chargement des donnees depuis le fichier json via le parseur boost
	read_json(input_tempfile, allParamsRoot);
	dataRoot = allParamsRoot.find("GRID_PARAMETERS")->second;


	/*
	 * initialisation des paramètres des contrôles
	 */
	gridParameters->FILE_PREFIX = dataRoot.get<string>("OUTPUT_FILE_PREFIX", "Model-");

	int dim = dataRoot.get<int>("STATE_DIMENSION", 1);
	gridParameters->GRID_METHOD = dataRoot.get<int>("GRID_METHOD", 1);

	gridParameters->DIM=(unsigned long long int) dim;

	gridParameters->GRID_TYPE=0;
	gridParameters->OMP_THREADS = nbOmpThreads;

	gridParameters->LIMINF = new double [dim];
	gridParameters->LIMSUP = new double [dim];

	for(unsigned long long int tabIndice=0;tabIndice<dim;tabIndice++)
	{
		gridParameters->LIMINF[tabIndice]=0.0;//valeurs par default
		gridParameters->LIMSUP[tabIndice]=1.0;//valeurs par default
	}
	this->readTabData(&dataRoot, gridParameters->LIMSUP, "STATE_MAX_VALUES", dim);
	this->readTabData(&dataRoot, gridParameters->LIMINF, "STATE_MIN_VALUES", dim);

	gridParameters->NBPOINTS=new unsigned long long int[dim];
	for(unsigned long long int tabIndice=0;tabIndice<dim;tabIndice++)
	{
		gridParameters->NBPOINTS[tabIndice]=2;//valeurs par default
	}

	this->readTabData(&dataRoot, gridParameters->NBPOINTS, "STATE_GRID_POINTS", dim);



	gridParameters->PERIODIC = new int [dim];

	if(dataRoot.find("STATE_PERIODIC")!=dataRoot.not_found())
	{
		this->readTabData(&dataRoot, gridParameters->PERIODIC, "STATE_PERIODIC", dim);
	}
	else
	{
		for(unsigned long long int tabIndice=0;tabIndice<dim;tabIndice++)
		{
			gridParameters->PERIODIC[tabIndice]=0;//valeurs par default
		}
	}

	if(dataRoot.find("GRID_MAIN_DIR")!=dataRoot.not_found())
	{
		gridParameters->GRID_MAIN_DIR = dataRoot.get<int>("GRID_MAIN_DIR", 0);
	}
	else
	{
		gridParameters->GRID_MAIN_DIR = 0;
	}

	gridParameters->SLICE_DIRS   =  new int [dim];
	if(dataRoot.find("SLICE_DIRECTIONS")!=dataRoot.not_found())
	{
		this->readTabData(&dataRoot, gridParameters->SLICE_DIRS, "SLICE_DIRECTIONS", dim);
	}
	else
	{
		for(unsigned long long int tabIndice=0;tabIndice<dim;tabIndice++)
		{
			gridParameters->SLICE_DIRS[tabIndice]=0;//valeurs par default
		}
	}

	gridParameters->SLICE_VALUES =  new double [dim];
	if(dataRoot.find("SLICE_LEVELS")!=dataRoot.not_found())
	{
		this->readTabData(&dataRoot, gridParameters->SLICE_VALUES, "SLICE_LEVELS", dim);
	}
	else
	{
		for(unsigned long long int tabIndice=0;tabIndice<dim;tabIndice++)
		{
			gridParameters->SLICE_VALUES[tabIndice]=0.0;//valeurs par default
		}
	}

	gridParameters->SORTIE_OK_INF = new int [dim];
	if(dataRoot.find("LOWER_LIMIT_IS_NOT_CONSTRAINT")!=dataRoot.not_found())
	{
		this->readTabData(&dataRoot, gridParameters->SORTIE_OK_INF, "LOWER_LIMIT_IS_NOT_CONSTRAINT", dim);
	}
	else
	{
		for(unsigned long long int tabIndice=0;tabIndice<dim;tabIndice++)
		{
			gridParameters->SORTIE_OK_INF[tabIndice]=0;//valeurs par default
		}
	}

	gridParameters->SORTIE_OK_SUP = new int [dim];
	if(dataRoot.find("UPPER_LIMIT_IS_NOT_CONSTRAINT")!=dataRoot.not_found())
	{
		this->readTabData(&dataRoot, gridParameters->SORTIE_OK_SUP, "UPPER_LIMIT_IS_NOT_CONSTRAINT", dim);
	}
	else
	{
		for(unsigned long long int tabIndice=0;tabIndice<dim;tabIndice++)
		{
			gridParameters->SORTIE_OK_SUP[tabIndice]=0;//valeurs par default
		}
	}
	cout<< " LECTURe SORTIE PARAMS finie test"<< endl;
	printVector(gridParameters->SORTIE_OK_SUP, dim);
	printVector(gridParameters->SORTIE_OK_INF, dim);

	cout<< "  REad grig params : FINISHED \n";
}

void ParametersManager::readSystemParametersFromJson()
{
	cout<< "  REad SYSTEM params \n";
	string line;
	ostringstream os;

	string tempStr, tmpStr;
	ptree allParamsRoot, dataRoot;
	string sedCommand;

	string input_tempfile="../INPUT/"+ parametersFileName;

	//- chargement des donnees depuis le fichier json via le parseur boost
	read_json(input_tempfile, allParamsRoot);
	dataRoot = allParamsRoot.find("SYSTEM_PARAMETERS")->second;

	systemParameters->COMPUTE_LC=dataRoot.get<int>("LIPSCHITZ_CONSTANT_COMPUTE_METHOD", 1);
	systemParameters->COMPUTE_MF=dataRoot.get<int>("DYN_BOUND_COMPUTE_METHOD", 1);

	if(dataRoot.find("LIPSCHITZ_CONSTANT")!=dataRoot.not_found())
	{
		systemParameters->LIP=dataRoot.get<double>("LIPSCHITZ_CONSTANT", 1.0);
	}
	else
	{
		systemParameters->LIP = 1.0;
	}

	if(dataRoot.find("DYN_BOUND")!=dataRoot.not_found())
	{
		systemParameters->MF=dataRoot.get<double>("DYN_BOUND", 1.0);
	}
	else
	{
		systemParameters->MF = 1.0;
	}

	/*
	 * Initialisation des paramètres de systèmes dynamique
	 */
	systemParameters->DYN_TYPE=dataRoot.get<int>("DYNAMICS_TYPE", 1);


	/*
	 * paramètres par défaut. Non utilisés ici
	 */
	systemParameters->globDeltat=dataRoot.get<bool>("IS_TIMESTEP_GLOBAL", false);
	systemParameters->maxTime=dataRoot.get<double>("TIME_HORIZON", 10.0);

	systemParameters->SCHEME=dataRoot.get<int>("TIME_DISCRETIZATION_SCHEME", 1);
	systemParameters->SCALE_PARAM=dataRoot.get<bool>("SCALING", false);
	if(algoParameters->COMPUTE_TMIN)
	{
		systemParameters->L_LIP=1.0;
		systemParameters->L_MF=systemParameters->maxTime;
	}
	else{
		systemParameters->L_LIP=dataRoot.get<double>("COST_LIPSCHITZ_CONSTANT", 1.0);
		systemParameters->L_MF=dataRoot.get<double>("COST_BOUND_CONSTANT", 1.0);
		cout<< "  lu cost constants "<<  systemParameters->L_LIP <<  " and " << systemParameters->L_MF << endl;
	}

	cout<< "  Read SYSTEM params  : FINISHED\n";
}


void ParametersManager::readAlgoParametersFromJson()
{

	cout<< "  Raad algo params \n";
	string line;
	ostringstream os;

	string tempStr, tmpStr;
	ptree allParamsRoot, dataRoot;
	string sedCommand;

	string input_tempfile="../INPUT/"+ parametersFileName;

	//- chargement des donnees depuis le fichier json via le parseur boost
	read_json(input_tempfile, allParamsRoot);
	dataRoot = allParamsRoot.find("ALGORITHM_PARAMETERS")->second;

	algoParameters->NB_OMP_THREADS=nbOmpThreads;

	algoParameters->FILE_PREFIX=gridParameters->FILE_PREFIX;
	algoParameters->TARGET_OR_DEPARTURE=dataRoot.get<int>("TARGET_OR_DEPARTURE", 1);
	algoParameters->COMPUTE_TMIN=dataRoot.get<int>("COMPUTE_MIN_TIME", 0);
	algoParameters->GRID_REFINMENTS_NUMBER = dataRoot.get<int>("GRID_REFINMENTS_NUMBER", 0);
	algoParameters->COMPUTE_SET = dataRoot.get<int>("COMPUTE_VIABLE_SET", 1);
	algoParameters->SET_TYPE = dataRoot.get<int>("SET_TYPE", 1);
	algoParameters->INTERATION_STOP_LEVEL = dataRoot.get<int>("ITERATION_STOP_LEVEL", 0);

	if((algoParameters->COMPUTE_SET==0)  & (algoParameters->GRID_REFINMENTS_NUMBER >0))
	{
		for(int k=0;k<gridParameters->DIM;k++)
		{
			for(int j=0;j<algoParameters->GRID_REFINMENTS_NUMBER;j++)
			{
				gridParameters->NBPOINTS[k]=2*gridParameters->NBPOINTS[k]-1;
			}
		}
	}

	algoParameters->NB_TRAJS = dataRoot.get<int>("NUMBER_OF_TRAJECTORIES", 0);

	algoParameters->TYPE_TRAJ = dataRoot.get<int>("TRAJECTORY_TYPE", 1);
	double * initPoints = new double[gridParameters->DIM * algoParameters->NB_TRAJS];
	if(algoParameters->NB_TRAJS >0)
	{
		if(dataRoot.find("INITIAL_POINTS")!=dataRoot.not_found())
		{
			ptree  pointsTree=dataRoot.find("INITIAL_POINTS")->second;
			int pointCounter = 0;
			for (ptree::value_type &point : pointsTree.get_child(""))
			{
				ptree pointRoot = point.second;
				cout<< " readig point coordinates with dimension "<< gridParameters->DIM<<endl;
				this->readTabData(&pointRoot, initPoints+pointCounter * gridParameters->DIM, "POINT_COORDINATES", gridParameters->DIM);
				pointCounter++;
			}
		}
	}
	algoParameters->INIT_POINTS = initPoints;

	unsigned long long int * initPoints_fd = new unsigned long long int[gridParameters->DIM * algoParameters->NB_TRAJS];
	if(dataRoot.find("INITIAL_POINTS_DISCRETE")!=dataRoot.not_found())
	{
		ptree  pointsTree=dataRoot.find("INITIAL_POINTS_DISCRETE")->second;
		int pointCounter = 0;
		for (ptree::value_type &point : pointsTree.get_child(""))
		{
			ptree pointRoot = point.second;
			this->readTabData(&pointRoot, initPoints_fd+pointCounter * gridParameters->DIM, "POINT_COORDINATES", gridParameters->DIM);
			pointCounter++;
		}
	}

	algoParameters->INIT_POINTS_FD = initPoints_fd;

	double * initValues = new double[algoParameters->NB_TRAJS];
	if(dataRoot.find("INITIAL_VALUES")!=dataRoot.not_found())
	{
		ptree  pointsTree=dataRoot.find("INITIAL_VALUES")->second;
		int pointCounter = 0;
		for (ptree::value_type &point : pointsTree.get_child(""))
		{
			initValues[pointCounter] = point.second.get_value<double>();
			pointCounter++;
		}
	}
	algoParameters->INIT_VALUES = initValues;

	double * initValues_fd = new double[algoParameters->NB_TRAJS];
	if(dataRoot.find("INITIAL_VALUES_DISCRETE")!=dataRoot.not_found())
	{
		ptree  pointsTree=dataRoot.find("INITIAL_VALUES_DISCRETE")->second;
		int pointCounter = 0;
		for (ptree::value_type &point : pointsTree.get_child(""))
		{
			initValues_fd[pointCounter] = point.second.get_value<double>();
			pointCounter++;
		}
	}
	algoParameters->INIT_VALUES_FD = initValues_fd;

	double * initControls = new double[controlParameters->DIMC * algoParameters->NB_TRAJS];
	if(dataRoot.find("INITIAL_CONTROLS")!=dataRoot.not_found())
	{
		ptree  pointsTree=dataRoot.find("INITIAL_CONTROLS")->second;
		int pointCounter = 0;
		for (ptree::value_type &point : pointsTree.get_child(""))
		{
			ptree pointRoot = point.second;
			this->readTabData(&pointRoot, initControls+pointCounter * controlParameters->DIMC, "CONTROL_COORDINATES", controlParameters->DIMC);
			pointCounter++;
		}
	}

	algoParameters->INIT_CONTROLS = initControls;


	algoParameters->INTERMEDIATE_SAVINGS = dataRoot.get<int>("INTERMEDIATE_SAVINGS", 0);
	algoParameters->SAVE_BOUNDARY = dataRoot.get<int>("SAVE_BOUNDARY", 0);
	algoParameters->SAVE_PROJECTION = dataRoot.get<int>("SAVE_PROJECTION", 0);
	algoParameters->SAVE_SLICE =dataRoot.get<int>("SAVE_SLICE", 0);
	algoParameters->SAVE_SLICE_BOUND = dataRoot.get<int>("SAVE_SLICE_BOUND", 0);
	algoParameters->PROJECTION=new unsigned long long int[gridParameters->DIM];
	for(unsigned long long int tabIndice=0;tabIndice<gridParameters->DIM;tabIndice++)
	{
		algoParameters->PROJECTION[tabIndice]=0;//valeurs par default
	}
	this->readTabData(&dataRoot, algoParameters->PROJECTION, "PROJECTION_AXIS", gridParameters->DIM);

	algoParameters->SAVE_SUBLEVEL = dataRoot.get<int>("SAVE_SUBLEVEL", 0);
	algoParameters->LEVEL=dataRoot.get<int>("LEVEL", 0);


	cout<< "  Read algo params  : finished\n";
}
void ParametersManager::readControlParametersFromJson()
{
	cout<< "  Read control parameter \n";
	string line;
	ostringstream os;

	string tempStr, tmpStr;
	ptree allParamsRoot, dataRoot;
	string sedCommand;

	string input_tempfile="../INPUT/"+ parametersFileName;
	cout<< " complete file name "<< input_tempfile << endl;
	//sedCommand="sed  -e 's/[ ]*#.*//g' -e '/^$/d' ../EXEC_DATA/"+ farmSimulationFile+" > "+input_tempfile;
	//unsigned long long int  zz=system(sedCommand.c_str()); zz++;

	//- chargement des donnees depuis le fichier json via le parseur boost
	read_json(input_tempfile, allParamsRoot);
	cout<< " all params root OK ";
	dataRoot = allParamsRoot.find("CONTROL_PARAMETERS")->second;


	/*
	 * initialisation des paramètres des contrôles
	 */
	int dimc = dataRoot.get<int>("CONTROL_DIMENSION", 1);
	controlParameters->DIMC=(unsigned long long int) dimc;
	/*
	 * paramètres par défaut. Non utilisés ici
	 */
	controlParameters->LIMINFC = new double [dimc];
	controlParameters->LIMSUPC = new double [dimc];
	for(unsigned long long int tabIndice=0;tabIndice<dimc;tabIndice++)
	{
		controlParameters->LIMINFC[tabIndice]=0.0;//valeurs par default
		controlParameters->LIMSUPC[tabIndice]=1.0;//valeurs par default
	}
	this->readTabData(&dataRoot, controlParameters->LIMSUPC, "CONTROL_MAX_VALUES", dimc);
	this->readTabData(&dataRoot, controlParameters->LIMINFC, "CONTROL_MIN_VALUES", dimc);

	controlParameters->NBPOINTSC=new unsigned long long int[dimc];
	for(unsigned long long int tabIndice=0;tabIndice<dimc;tabIndice++)
	{
		controlParameters->NBPOINTSC[tabIndice]=1;//valeurs par default
	}
	this->readTabData(&dataRoot, controlParameters->NBPOINTSC, "CONTROL_GRID_POINTS", dimc);

	int dimc_ty;
	if(dataRoot.find("CONTROL_TYCHASTIC_DIMENSION")!=dataRoot.not_found())
	{
		dimc_ty= dataRoot.get<int>("CONTROL_TYCHASTIC_DIMENSION", 1);
		controlParameters->DIM_TY=(unsigned long long int) dimc_ty;
		/*
		 * paramètres par défaut. Non utilisés ici
		 */
		controlParameters->LIMINF_TY = new double [dimc_ty];
		controlParameters->LIMSUP_TY = new double [dimc_ty];
		for(unsigned long long int tabIndice=0;tabIndice<dimc_ty;tabIndice++)
		{
			controlParameters->LIMINF_TY[tabIndice]=0.0;//valeurs par default
			controlParameters->LIMSUP_TY[tabIndice]=1.0;//valeurs par default
		}
		this->readTabData(&dataRoot, controlParameters->LIMSUP_TY, "CONTROL_TY_MAX_VALUES", dimc_ty);
		this->readTabData(&dataRoot, controlParameters->LIMINF_TY, "CONTROL_TY_MIN_VALUES", dimc_ty);

		controlParameters->NBPOINTS_TY=new unsigned long long int[dimc_ty];
		for(unsigned long long int tabIndice=0;tabIndice<dimc_ty;tabIndice++)
		{
			controlParameters->NBPOINTS_TY[tabIndice]=1;//valeurs par default
		}
		this->readTabData(&dataRoot, controlParameters->NBPOINTS_TY, "CONTROL_TY_GRID_POINTS", dimc_ty);
	}
	else
	{
		dimc_ty= 0;
		controlParameters->DIM_TY=0;
		/*
		 * paramètres par défaut. Non utilisés ici
		 */
		controlParameters->LIMINF_TY = new double [dimc_ty];
		controlParameters->LIMSUP_TY = new double [dimc_ty];
		controlParameters->NBPOINTS_TY=new unsigned long long int[dimc_ty];
	}

}

void ParametersManager::readTabData (ptree *dataRoot, unsigned long long int * target, string label, int nbElements )
{
	if(dataRoot->find(label)!=dataRoot->not_found())
	{
		ptree  tabTree=dataRoot->find(label)->second;
		ptree::const_iterator it = tabTree.begin();
		for(unsigned long long int tabIndice=0;tabIndice<nbElements;tabIndice++)
		{
			target[tabIndice]=(unsigned long long int)(*it).second.get_value<int>();
			it++;

			cout<< " lu  depuis "<<  label << " "<<target[tabIndice]<<endl;
		}
	}
}

void ParametersManager::readTabData (ptree *dataRoot, int * target, string label, int nbElements )
{
	if(dataRoot->find(label)!=dataRoot->not_found())
	{
		ptree  tabTree=dataRoot->find(label)->second;
		ptree::const_iterator it = tabTree.begin();
		for(unsigned long long int tabIndice=0;tabIndice<nbElements;tabIndice++)
		{
			target[tabIndice]=(*it).second.get_value<int>();
			it++;

			cout<< " lu  depuis "<<  label << " "<<target[tabIndice]<<endl;
		}
	}
}

void ParametersManager::readTabData (ptree *dataRoot, double * target, string label, int nbElements )
{
	if(dataRoot->find(label)!=dataRoot->not_found())
	{
		ptree  tabTree=dataRoot->find(label)->second;
		ptree::const_iterator it = tabTree.begin();
		for(unsigned long long int tabIndice=0;tabIndice<nbElements;tabIndice++)
		{
			target[tabIndice]=(*it).second.get_value<double>();
			it++;

			cout<< " lu  depuis "<<  label << " "<<target[tabIndice]<<endl;
		}
	}
}
ParametersManager::~ParametersManager() {
	// TODO Auto-generated destructor stub
}

