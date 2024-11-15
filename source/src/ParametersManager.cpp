/*
 * ParametersManager.cpp
 *
 *    VIABLAB : a numerical library for Mathematical Viability Computations
 *    Copyright (C) <2020>  <Anna DESILLES, LASTRE>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Created on: 7 juil. 2021
 *      Author: Anna DESILLES
 */
#include "../include/ParametersManager.h"

ParametersManager::ParametersManager()
    {
    // TODO Auto-generated constructor stub

    }

ParametersManager::ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams *cp, systemParams *sp)
    {
    gridParameters = gp;
    algoParameters = avp;
    controlParameters = cp;
    systemParameters = sp;
    }

ParametersManager::ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams *cp, systemParams *sp, int nbThreads, string paramsFile)
    {
    gridParameters = gp;
    algoParameters = avp;
    controlParameters = cp;
    systemParameters = sp;
    nbOmpThreads = nbThreads;
    parametersFileName = paramsFile;
    ptree allParamsRoot, dataRoot;


    this->readControlParametersFromJson();
    this->readGridParametersFromJson();
    this->readAlgoParametersFromJson();
    this->readSystemParametersFromJson();
    }

gridParams* ParametersManager::getGridParameters()
    {
    return gridParameters;
    }
algoViabiParams* ParametersManager::getAlgoParameters()
    {
    return algoParameters;
    }
controlParams* ParametersManager::getControlParameters()
    {
    return controlParameters;
    }
systemParams* ParametersManager::getSystemParameters()
    {
    return systemParameters;
    }

void ParametersManager::readGridParametersFromJson()
    {
    string line;
    ostringstream os;

    string tempStr, tmpStr;
    ptree allParamsRoot, dataRoot;
    string sedCommand;

    string input_tempfile = "../INPUT/" + parametersFileName;
    spdlog::info("[ParametrManager] : Reading of Grid parameters");
    spdlog::info("[ParametrManager] : Input json file {}", input_tempfile.c_str());

    //- loading data from input file
    read_json(input_tempfile, allParamsRoot);
    dataRoot = allParamsRoot.find("GRID_PARAMETERS")->second;

    gridParameters->FILE_PREFIX = dataRoot.get<string>("OUTPUT_FILE_PREFIX", "Model-");

    int dim = dataRoot.get<int>("STATE_DIMENSION", 1);
    gridParameters->GRID_METHOD = dataRoot.get<int>("GRID_METHOD", 1);

    gridParameters->DIM = (unsigned long long int) dim;

    gridParameters->GRID_TYPE = 0;
    gridParameters->OMP_THREADS = nbOmpThreads;

    gridParameters->LIMINF = new double[dim];
    gridParameters->LIMSUP = new double[dim];

    for (int tabIndice = 0; tabIndice < dim; tabIndice++)
	{
	gridParameters->LIMINF[tabIndice] = 0.0;	//default values
	gridParameters->LIMSUP[tabIndice] = 1.0;	//default values
	}
    readTabData(&dataRoot, gridParameters->LIMSUP, "STATE_MAX_VALUES", dim);
    readTabData(&dataRoot, gridParameters->LIMINF, "STATE_MIN_VALUES", dim);

    gridParameters->NBPOINTS = new unsigned long long int[dim];
    for (int tabIndice = 0; tabIndice < dim; tabIndice++)
	{
	gridParameters->NBPOINTS[tabIndice] = 2;	//default values
	}

    this->readTabData(&dataRoot, gridParameters->NBPOINTS, "STATE_GRID_POINTS", dim);

    gridParameters->PERIODIC = new int[dim];

    if (dataRoot.find("STATE_PERIODIC") != dataRoot.not_found())
	{
	this->readTabData(&dataRoot, gridParameters->PERIODIC, "STATE_PERIODIC", dim);
	}
    else
	{
	for (auto tabIndice = 0; tabIndice < dim; tabIndice++)
	    {
	    gridParameters->PERIODIC[tabIndice] = 0;	//default values
	    }
	}

    if (dataRoot.find("GRID_MAIN_DIR") != dataRoot.not_found())
	{
	gridParameters->GRID_MAIN_DIR = dataRoot.get<int>("GRID_MAIN_DIR", 0);
	}
    else
	{
	gridParameters->GRID_MAIN_DIR = 0;
	}

    gridParameters->SLICE_DIRS = new int[dim];
    if (dataRoot.find("SLICE_DIRECTIONS") != dataRoot.not_found())
	{
	readTabData(&dataRoot, gridParameters->SLICE_DIRS, "SLICE_DIRECTIONS", dim);
	}
    else
	{
	for (auto tabIndice = 0; tabIndice < dim; tabIndice++)
	    {
	    gridParameters->SLICE_DIRS[tabIndice] = 0;	//Default Values
	    }
	}

    gridParameters->SLICE_VALUES = new double[dim];
    if (dataRoot.find("SLICE_LEVELS") != dataRoot.not_found())
	{
	readTabData(&dataRoot, gridParameters->SLICE_VALUES, "SLICE_LEVELS", dim);
	}
    else
	{
	for (auto tabIndice = 0; tabIndice < dim; tabIndice++)
	    {
	    gridParameters->SLICE_VALUES[tabIndice] = 0.0;	//Default Values
	    }
	}

    gridParameters->SLICE_VALUES_FD = new unsigned long long int[dim];
    if (dataRoot.find("SLICE_LEVELS_DISCRETE") != dataRoot.not_found())
	{
	readTabData(&dataRoot, gridParameters->SLICE_VALUES_FD, "SLICE_LEVELS_DISCRETE", dim);
	}
    else
	{
	for (auto tabIndice = 0; tabIndice < dim; tabIndice++)
	    {
	    gridParameters->SLICE_VALUES_FD[tabIndice] = 0;	//Default Values
	    }
	}

    gridParameters->SORTIE_OK_INF = new int[dim];
    if (dataRoot.find("LOWER_LIMIT_IS_NOT_CONSTRAINT") != dataRoot.not_found())
	{
	readTabData(&dataRoot, gridParameters->SORTIE_OK_INF, "LOWER_LIMIT_IS_NOT_CONSTRAINT", dim);
	}
    else
	{
	for (auto tabIndice = 0; tabIndice < dim; tabIndice++)
	    {
	    gridParameters->SORTIE_OK_INF[tabIndice] = 0;	//Default Values
	    }
	}

    gridParameters->SORTIE_OK_SUP = new int[dim];
    if (dataRoot.find("UPPER_LIMIT_IS_NOT_CONSTRAINT") != dataRoot.not_found())
	{
	readTabData(&dataRoot, gridParameters->SORTIE_OK_SUP, "UPPER_LIMIT_IS_NOT_CONSTRAINT", dim);
	}
    else
	{
	for (auto tabIndice = 0; tabIndice < dim; tabIndice++)
	    {
	    gridParameters->SORTIE_OK_SUP[tabIndice] = 0;	//Default Values
	    }
	}
    logVector("[ParametrManager] : Sortie sup: ", gridParameters->SORTIE_OK_SUP, dim);
    logVector("[ParametrManager] : Sortie inf: ", gridParameters->SORTIE_OK_INF, dim);

    spdlog::info("[ParametrManager] : Read grid params : FINISHED");
    }

void ParametersManager::readSystemParametersFromJson()
    {
    string line;
    ostringstream os;

    string tempStr, tmpStr;
    ptree allParamsRoot, dataRoot;
    string sedCommand;

    string input_tempfile = "../INPUT/" + parametersFileName;

    spdlog::info("[ParametrManager] : Reading of System parameters");
    spdlog::info("[ParametrManager] : Input json file {}", input_tempfile.c_str());

    //- loading data from input file
    read_json(input_tempfile, allParamsRoot);
    dataRoot = allParamsRoot.find("SYSTEM_PARAMETERS")->second;

    systemParameters->COMPUTE_LC = dataRoot.get<int>("LIPSCHITZ_CONSTANT_COMPUTE_METHOD", 1);
    systemParameters->COMPUTE_MF = dataRoot.get<int>("DYN_BOUND_COMPUTE_METHOD", 1);

    if (dataRoot.find("LIPSCHITZ_CONSTANT") != dataRoot.not_found())
	{
	systemParameters->LIP = dataRoot.get<double>("LIPSCHITZ_CONSTANT", 1.0);
	}
    else
	{
	systemParameters->LIP = 1.0;
	}

    if (dataRoot.find("DYN_BOUND") != dataRoot.not_found())
	{
	systemParameters->MF = dataRoot.get<double>("DYN_BOUND", 1.0);
	}
    else
	{
	systemParameters->MF = 1.0;
	}

    spdlog::debug("[ParametrManager] : Dyn bounds estimation method {}, Dyn bounds default value  {}", systemParameters->COMPUTE_MF,
	    systemParameters->MF);
    spdlog::debug("[ParametrManager] : Lipschitz constant estimation method {}, Lipschitz constant value {}", systemParameters->COMPUTE_LC,
	    systemParameters->LIP);

    systemParameters->DYN_TYPE = dataRoot.get<int>("DYNAMICS_TYPE", 1);

    systemParameters->globDeltat = dataRoot.get<bool>("IS_TIMESTEP_GLOBAL", false);
    systemParameters->maxTime = dataRoot.get<double>("TIME_HORIZON", 10.0);

    systemParameters->SCHEME = dataRoot.get<int>("TIME_DISCRETIZATION_SCHEME", 1);
    systemParameters->SCALE_PARAM = dataRoot.get<bool>("SCALING", false);

    if (dataRoot.find("COST_LIPSCHITZ_CONSTANT") != dataRoot.not_found())
	{
	systemParameters->L_LIP = dataRoot.get<double>("COST_LIPSCHITZ_CONSTANT", 1.0);
	}
    else
	{
	if (algoParameters->COMPUTE_TMIN)
	    {
	    systemParameters->L_LIP = systemParameters->maxTime;
	    }
	else
	    {
	    systemParameters->L_LIP = 1.0;
	    }
	}

    systemParameters->L_MF = dataRoot.get<double>("COST_BOUND_CONSTANT", 1.0);

    spdlog::debug("[ParametrManager] : Cost function bound{}, Cost function Lipschitz contant {}", systemParameters->L_MF, systemParameters->L_LIP);
    spdlog::info("[ParametrManager] : Read system params : FINISHED");
    }

void ParametersManager::readAlgoParametersFromJson()
    {

    string line;
    ostringstream os;

    string tempStr, tmpStr;
    ptree allParamsRoot, dataRoot;
    string sedCommand;

    string input_tempfile = "../INPUT/" + parametersFileName;
    spdlog::info("[ParametrManager] : Reading of Algorithm parameters");
    spdlog::info("[ParametrManager] : Input json file {}", input_tempfile.c_str());

    //- loading data from input file
    read_json(input_tempfile, allParamsRoot);
    dataRoot = allParamsRoot.find("ALGORITHM_PARAMETERS")->second;

    algoParameters->NB_OMP_THREADS = nbOmpThreads;

    algoParameters->FILE_PREFIX = gridParameters->FILE_PREFIX;
    algoParameters->TARGET_OR_DEPARTURE = dataRoot.get<int>("TARGET_OR_DEPARTURE", 1);
    algoParameters->COMPUTE_TMIN = dataRoot.get<int>("COMPUTE_MIN_TIME", 0);
    algoParameters->GRID_REFINMENTS_NUMBER = dataRoot.get<int>("GRID_REFINMENTS_NUMBER", 0);
    algoParameters->COMPUTE_SET = dataRoot.get<int>("COMPUTE_VIABLE_SET", 1);
    algoParameters->SET_TYPE = dataRoot.get<int>("SET_TYPE", 1);
    algoParameters->INTERATION_STOP_LEVEL = dataRoot.get<int>("ITERATION_STOP_LEVEL", 0);
    algoParameters->SAVE_SUBLEVEL = dataRoot.get<int>("SAVE_SUBLEVEL", 0);
    algoParameters->LEVEL = dataRoot.get<double>("LEVEL", 0.0);
    algoParameters->SAVE_VIAB_LIGHT = dataRoot.get<int>("SAVE_VIABSET_LIGHT", 0);

    if ((algoParameters->COMPUTE_SET == 0) & (algoParameters->GRID_REFINMENTS_NUMBER > 0))
	{
	for (unsigned int k = 0; k < gridParameters->DIM; k++)
	    {
	    for (int j = 0; j < algoParameters->GRID_REFINMENTS_NUMBER; j++)
		{
		gridParameters->NBPOINTS[k] = 2 * gridParameters->NBPOINTS[k] - 1;
		}
	    }
	}

    algoParameters->NB_TRAJS = dataRoot.get<int>("NUMBER_OF_TRAJECTORIES", 0);

    algoParameters->TYPE_TRAJ = dataRoot.get<int>("TRAJECTORY_TYPE", 1);
    double *initPoints = new double[gridParameters->DIM * algoParameters->NB_TRAJS];
    if (algoParameters->NB_TRAJS > 0)
	{
	if (dataRoot.find("INITIAL_POINTS") != dataRoot.not_found())
	    {
	    readDoubleTabData(&dataRoot, initPoints, "INITIAL_POINTS", algoParameters->NB_TRAJS, gridParameters->DIM);
	    }
	}
    algoParameters->INIT_POINTS = initPoints;

    unsigned long long int *initPoints_fd = new unsigned long long int[gridParameters->DIM * algoParameters->NB_TRAJS];
    if (dataRoot.find("INITIAL_POINTS_DISCRETE") != dataRoot.not_found())
	{
	readDoubleTabData(&dataRoot, initPoints_fd, "INITIAL_POINTS_DISCRETE", algoParameters->NB_TRAJS, gridParameters->DIM);
	}

    algoParameters->INIT_POINTS_FD = initPoints_fd;

    double *initValues = new double[algoParameters->NB_TRAJS];
    if (dataRoot.find("INITIAL_VALUES") != dataRoot.not_found())
	{
	ptree pointsTree = dataRoot.find("INITIAL_VALUES")->second;
	int pointCounter = 0;
	for (ptree::value_type &point : pointsTree.get_child(""))
	    {
	    initValues[pointCounter] = point.second.get_value<double>();
	    pointCounter++;
	    }
	}
    algoParameters->INIT_VALUES = initValues;

    double *initValues_fd = new double[algoParameters->NB_TRAJS];
    if (dataRoot.find("INITIAL_VALUES_DISCRETE") != dataRoot.not_found())
	{
	ptree pointsTree = dataRoot.find("INITIAL_VALUES_DISCRETE")->second;
	int pointCounter = 0;
	for (ptree::value_type &point : pointsTree.get_child(""))
	    {
	    initValues_fd[pointCounter] = point.second.get_value<double>();
	    pointCounter++;
	    }
	}
    algoParameters->INIT_VALUES_FD = initValues_fd;

    double *initControls = new double[controlParameters->DIMC * algoParameters->NB_TRAJS];
    if (dataRoot.find("INITIAL_CONTROLS") != dataRoot.not_found())
	{
	readDoubleTabData(&dataRoot, initControls, "INITIAL_CONTROLS", algoParameters->NB_TRAJS, controlParameters->DIMC);
	}

    algoParameters->INIT_CONTROLS = initControls;

    algoParameters->INTERMEDIATE_SAVINGS = dataRoot.get<int>("INTERMEDIATE_SAVINGS", 0);
    algoParameters->SAVE_BOUNDARY = dataRoot.get<int>("SAVE_BOUNDARY", 0);
    algoParameters->SAVE_PROJECTION = dataRoot.get<int>("SAVE_PROJECTION", 0);
    algoParameters->SAVE_SLICE = dataRoot.get<int>("SAVE_SLICE", 0);
    algoParameters->SAVE_SLICE_BOUND = dataRoot.get<int>("SAVE_SLICE_BOUND", 0);
    algoParameters->PROJECTION = new unsigned long long int[gridParameters->DIM];
    for (unsigned long long int tabIndice = 0; tabIndice < gridParameters->DIM; tabIndice++)
	{
	algoParameters->PROJECTION[tabIndice] = 0;	//Default Values
	}
    this->readTabData(&dataRoot, algoParameters->PROJECTION, "PROJECTION_AXIS", gridParameters->DIM);

    spdlog::info("[ParametrManager] : Read algorithm params : FINISHED");
    }
void ParametersManager::readControlParametersFromJson()
    {

    string line;
    ostringstream os;

    string tempStr, tmpStr;
    ptree allParamsRoot, dataRoot;
    string sedCommand;

    string input_tempfile = "../INPUT/" + parametersFileName;
    spdlog::info("[ParametrManager] : Reading of controls parameters");
    spdlog::info("[ParametrManager] : Input json file {}", input_tempfile.c_str());

    //- loading data from input file
    read_json(input_tempfile, allParamsRoot);
    dataRoot = allParamsRoot.find("CONTROL_PARAMETERS")->second;

    /*
     * initialisation des paramètres des contrôles
     */
    int dimc = 0;
    if (dataRoot.find("CONTROL_DIMENSION") != dataRoot.not_found())
	{
	dimc = dataRoot.get<int>("CONTROL_DIMENSION", 1);

	controlParameters->DIMC = (unsigned long long int) dimc;
	/*
	 * paramètres par défaut. Non utilisés ici
	 */
	controlParameters->LIMINFC = new double[dimc];
	controlParameters->LIMSUPC = new double[dimc];
	for (auto tabIndice = 0; tabIndice < dimc; tabIndice++)
	    {
	    controlParameters->LIMINFC[tabIndice] = 0.0;	//Default Values
	    controlParameters->LIMSUPC[tabIndice] = 1.0;	//Default Values
	    }
	this->readTabData(&dataRoot, controlParameters->LIMSUPC, "CONTROL_MAX_VALUES", dimc);
	this->readTabData(&dataRoot, controlParameters->LIMINFC, "CONTROL_MIN_VALUES", dimc);

	controlParameters->NBPOINTSC = new unsigned long long int[dimc];
	for (auto tabIndice = 0; tabIndice < dimc; tabIndice++)
	    {
	    controlParameters->NBPOINTSC[tabIndice] = 1;	//Default Values
	    }
	this->readTabData(&dataRoot, controlParameters->NBPOINTSC, "CONTROL_GRID_POINTS", dimc);
	}
    else
	{
	dimc = 0;
	controlParameters->DIMC = 0;
	controlParameters->LIMINFC = new double[dimc];
	controlParameters->LIMSUPC = new double[dimc];
	controlParameters->NBPOINTSC = new unsigned long long int[dimc];
	}
    int dimc_ty;
    if (dataRoot.find("CONTROL_TYCHASTIC_DIMENSION") != dataRoot.not_found())
	{
	dimc_ty = dataRoot.get<int>("CONTROL_TYCHASTIC_DIMENSION", 1);
	controlParameters->DIM_TY = (unsigned long long int) dimc_ty;
	/*
	 * initializing with empty tabs, default values, not used in this case
	 */
	controlParameters->LIMINF_TY = new double[dimc_ty];
	controlParameters->LIMSUP_TY = new double[dimc_ty];
	for (auto tabIndice = 0; tabIndice < dimc_ty; tabIndice++)
	    {
	    controlParameters->LIMINF_TY[tabIndice] = 0.0;	//Default Values
	    controlParameters->LIMSUP_TY[tabIndice] = 1.0;	//Default Values
	    }
	this->readTabData(&dataRoot, controlParameters->LIMSUP_TY, "CONTROL_TY_MAX_VALUES", dimc_ty);
	this->readTabData(&dataRoot, controlParameters->LIMINF_TY, "CONTROL_TY_MIN_VALUES", dimc_ty);

	controlParameters->NBPOINTS_TY = new unsigned long long int[dimc_ty];
	for (auto tabIndice = 0; tabIndice < dimc_ty; tabIndice++)
	    {
	    controlParameters->NBPOINTS_TY[tabIndice] = 1;	//Default Values
	    }
	this->readTabData(&dataRoot, controlParameters->NBPOINTS_TY, "CONTROL_TY_GRID_POINTS", dimc_ty);
	}
    else
	{
	dimc_ty = 0;
	controlParameters->DIM_TY = 0;
	/*
	 * initializing with empty tabs, default values, not used in this case
	 */
	controlParameters->LIMINF_TY = new double[dimc_ty];
	controlParameters->LIMSUP_TY = new double[dimc_ty];
	controlParameters->NBPOINTS_TY = new unsigned long long int[dimc_ty];
	}
    spdlog::info("[ParametrManager] : Read control : FINISHED");
    }

void ParametersManager::readTabData(ptree *dataRoot, unsigned long long int *target, string label, int nbElements)
    {
    if (dataRoot->find(label) != dataRoot->not_found())
	{
	ptree tabTree = dataRoot->find(label)->second;
	ptree::const_iterator it = tabTree.begin();
	for (auto tabIndice = 0; tabIndice < nbElements; tabIndice++)
	    {
	    target[tabIndice] = (unsigned long long int) (*it).second.get_value<int>();
	    it++;
	    }
	logVector("[ParametrManager] : Read from table " + label, target, nbElements);
	}
    }

void ParametersManager::readDoubleTabData(ptree *dataRoot, unsigned long long int *target, string label, int nbElements, int dim)
    {
    if (dataRoot->find(label) != dataRoot->not_found())
	{
	int cptPoints = 0;
	BOOST_FOREACH(ptree::value_type & rowPair, dataRoot->find(label)->second)
	    {
	    if (cptPoints >= nbElements)
		{
		break;
		}
	    int cptCoords = 0;
	    BOOST_FOREACH(ptree::value_type & itemPair, rowPair.second)
		{
		if(cptCoords >= dim)
		    {
		    break;
		    }
		target[cptPoints * dim + cptCoords] = (unsigned long long int) itemPair.second.get_value<int>();
		cptCoords++;
		}
	    cptPoints++;
	    }
	}
    logVector("[ParametrManager] : Read from table " + label, target, nbElements * dim);
    }

void ParametersManager::readTabData(ptree *dataRoot, int *target, string label, int nbElements)
    {
    if (dataRoot->find(label) != dataRoot->not_found())
	{
	ptree tabTree = dataRoot->find(label)->second;
	ptree::const_iterator it = tabTree.begin();
	for (auto tabIndice = 0; tabIndice < nbElements; tabIndice++)
	    {
	    target[tabIndice] = (*it).second.get_value<int>();
	    it++;
	    }
	logVector("[ParametrManager] : Read from table " + label, target, nbElements);
	}
    }

void ParametersManager::readDoubleTabData(ptree *dataRoot, int *target, string label, int nbElements, int dim)
    {
    if (dataRoot->find(label) != dataRoot->not_found())
	{
	int cptPoints = 0;
	BOOST_FOREACH(ptree::value_type & rowPair, dataRoot->find(label)->second)
	    {
	    if (cptPoints >= nbElements)
		{
		break;
		}
	    int cptCoords = 0;
	    BOOST_FOREACH(ptree::value_type & itemPair, rowPair.second)
		{
		if(cptCoords >= dim)
		    {
		    break;
		    }
		target[cptPoints * dim + cptCoords] = itemPair.second.get_value<int>();
		cptCoords++;
		}
	    cptPoints++;
	    }
	}
    logVector("[ParametrManager] : Read from table " + label, target, nbElements * dim);
    }

void ParametersManager::readDoubleTabData(ptree *dataRoot, double *target, string label, int nbElements, int dim)
    {
    if (dataRoot->find(label) != dataRoot->not_found())
	{
	int cptPoints = 0;
	BOOST_FOREACH(ptree::value_type & rowPair, dataRoot->find(label)->second)
	    {
	    if (cptPoints >= nbElements)
		{
		break;
		}
	    int cptCoords = 0;
	    BOOST_FOREACH(ptree::value_type & itemPair, rowPair.second)
		{
		if(cptCoords >= dim)
		    {
		    break;
		    }
		target[cptPoints * dim + cptCoords] = itemPair.second.get_value<double>();
		cptCoords++;
		}
	    cptPoints++;
	    }
	}
    logVector("[ParametrManager] : Read from table " + label, target, nbElements * dim);
    }

void ParametersManager::readTabData(ptree *dataRoot, double *target, string label, int nbElements)
    {
    if (dataRoot->find(label) != dataRoot->not_found())
	{
	ptree tabTree = dataRoot->find(label)->second;
	ptree::const_iterator it = tabTree.begin();
	for (auto tabIndice = 0; tabIndice < nbElements; tabIndice++)
	    {
	    target[tabIndice] = (*it).second.get_value<double>();
	    it++;
	    }
	logVector("[ParametrManager] : Read from table " + label, target, nbElements);
	}
    }
ParametersManager::~ParametersManager()
    {
    // TODO Auto-generated destructor stub
    }

