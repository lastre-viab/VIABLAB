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
#include "../include/Params.h"
#include "../include/ParametersManager.h"
#include "../include/WeakDeclarations.h"

gridParams gp =
	{
	};
algoViabiParams avp =
	{
	};
controlParams cp =
	{
	};
systemParams sp =
	{
	};

template<typename T>
void ParametersManager::readTabData(ptree *dataRoot, T *target, const string &label, int nbElements, const T &defaultValue)
    {
    ptree::assoc_iterator labelPos = dataRoot->find(label);
    if (labelPos != dataRoot->not_found())
	{
	ptree tabTree = labelPos->second;
	if (tabTree.empty())
	    {
	    throw boost::property_tree::ptree_bad_data("The data read from the JSON file is not an array", tabTree);
	    }
	ptree::const_iterator it = tabTree.begin();
	for (auto tabIndice = 0; tabIndice < nbElements; tabIndice++)
	    {
	    target[tabIndice] = (*it).second.get_value < T > (defaultValue);
	    it++;
	    }
	logVector("[ParametrManager] : Read from table " + label, target, nbElements);
	}
    else
	{
	std::fill(target, target + nbElements, defaultValue);
	}
    }

template<typename T>
void ParametersManager::readTabDataSkipInvalid(ptree *dataRoot, T *target, const string &label, int &nbElements)
    {
    if (dataRoot->find(label) != dataRoot->not_found())
	{
	ptree tabTree = dataRoot->find(label)->second;
	if (tabTree.empty())
	    {
	    throw boost::property_tree::ptree_bad_data("Not an array", tabTree);
	    }

	int tabIndice = 0;
	int ptreeValuesRead = 0;
	for (ptree::const_iterator it = tabTree.begin(); it != tabTree.end(); ++it)
	    {
	    try
	    {
		target[tabIndice] = (*it).second.get_value<T>();
		tabIndice++;
	    }
	    catch (const std::exception &e)
		{
		spdlog::warn("Invalid value read in array named \"{}\" at index {}, skipping to next value", label, ptreeValuesRead);
		nbElements--;
		}
	    ptreeValuesRead++;
	    }
	logVector("[ParametrManager] : Read from table " + label, target, nbElements);
	}
    }

ParametersManager::ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams *cp, systemParams *sp, int nbThreads,
	const string &paramsFile, void *mh)
    {
    gridParameters = gp;
    algoParameters = avp;
    controlParameters = cp;
    systemParameters = sp;
    nbOmpThreads = nbThreads;
    parametersFileName = paramsFile;
    modelHandle = mh;

    ptree allParamsRoot;
    string input_tempfile = "../INPUT/" + parametersFileName;
    spdlog::info("[ParametrManager] : Reading input json file {}", input_tempfile.c_str());
    read_json(input_tempfile, allParamsRoot);

    this->readControlParametersFromJson(allParamsRoot);
    this->readGridParametersFromJson(allParamsRoot);
    this->readAlgoParametersFromJson(allParamsRoot);
    this->readSystemParametersFromJson(allParamsRoot);
    this->readTrajectoryParametersListFromJson(allParamsRoot);
    this->readModelParametersFromJson(allParamsRoot);
    }

const gridParams* ParametersManager::getGridParameters() const
    {
    return gridParameters;
    }
const algoViabiParams* ParametersManager::getAlgoParameters() const
    {
    return algoParameters;
    }
const controlParams* ParametersManager::getControlParameters() const
    {
    return controlParameters;
    }
const systemParams* ParametersManager::getSystemParameters() const
    {
    return systemParameters;
    }
const std::vector<trajectoryParams>& ParametersManager::getTrajectoryParametersList() const
    {
    return trajectoryParametersList;
    }
std::vector<trajectoryParams>& ParametersManager::getTrajectoryParametersList()
    {
    return trajectoryParametersList;
    }
int ParametersManager::getNbTrajectories() const
    {
    return trajectoryParametersList.size();
    }
const modelParams* ParametersManager::getModelParameters() const
    {
    return &modelParameters;
    }
void* ParametersManager::getModelHandle(void)
    {
    return modelHandle;
    }

void ParametersManager::readGridParametersFromJson(ptree &allParamsRoot)
    {
    string line;
    ostringstream os;

    string tempStr, tmpStr;
    ptree dataRoot;
    string sedCommand;

    spdlog::info("[ParametrManager] : Reading of Grid parameters");
    dataRoot = allParamsRoot.find("GRID_PARAMETERS")->second;

    gridParameters->FILE_PREFIX = dataRoot.get < string > ("OUTPUT_FILE_PREFIX", "Model-");

    gridParameters->GRID_METHOD = dataRoot.get < GridMethod > ("GRID_METHOD", BS);
    bool isHybridGrid = (gridParameters->GRID_METHOD == HBS) || (gridParameters->GRID_METHOD == HMM);
    gridParameters->IS_HYBRID = isHybridGrid;
    spdlog::info("[ParametrManager] : is hybrid {}",isHybridGrid);
    int dim;

    if (!isHybridGrid)
	{
	dim = dataRoot.get<int>("STATE_DIMENSION", 1);

	gridParameters->DIM = (unsigned long long int) dim;

	gridParameters->GRID_TYPE = 0;
	gridParameters->OMP_THREADS = nbOmpThreads;

	gridParameters->LIMINF = new double[dim];
	gridParameters->LIMSUP = new double[dim];

	readTabData(&dataRoot, gridParameters->LIMSUP, "STATE_MAX_VALUES", dim, 1.0);
	readTabData(&dataRoot, gridParameters->LIMINF, "STATE_MIN_VALUES", dim, 0.0);

	gridParameters->NBPOINTS = new unsigned long long int[dim];
	readTabData(&dataRoot, gridParameters->NBPOINTS, "STATE_GRID_POINTS", dim, 2ULL);

	gridParameters->PERIODIC = new bool[dim];
	readTabData(&dataRoot, gridParameters->PERIODIC, "STATE_PERIODIC", dim, false);

	gridParameters->GRID_MAIN_DIR = dataRoot.get<int>("GRID_MAIN_DIR", 0);

	gridParameters->SLICE_DIRS = new bool[dim];
	readTabData(&dataRoot, gridParameters->SLICE_DIRS, "SLICE_DIRECTIONS", dim, false);

	gridParameters->SLICE_VALUES = new double[dim];
	readTabData(&dataRoot, gridParameters->SLICE_VALUES, "SLICE_LEVELS", dim, 0.0);

	gridParameters->SLICE_VALUES_FD = new unsigned long long int[dim];
	readTabData(&dataRoot, gridParameters->SLICE_VALUES_FD, "SLICE_LEVELS_DISCRETE", dim, 0ULL);

	gridParameters->SORTIE_OK_INF = new bool[dim];
	readTabData(&dataRoot, gridParameters->SORTIE_OK_INF, "LOWER_LIMIT_IS_NOT_CONSTRAINT", dim, false);

	gridParameters->SORTIE_OK_SUP = new bool[dim];
	readTabData(&dataRoot, gridParameters->SORTIE_OK_SUP, "UPPER_LIMIT_IS_NOT_CONSTRAINT", dim, false);

	logVector("[ParametrManager] : Sortie sup: ", gridParameters->SORTIE_OK_SUP, dim);
	logVector("[ParametrManager] : Sortie inf: ", gridParameters->SORTIE_OK_INF, dim);

	}
    else
	{
	if (dataRoot.find("HYBRID_CONTINUOUS_SUBGRID") == dataRoot.not_found())
	    {
	    spdlog::error("[ParametrManager] : Hybrid continuous subgrid is missing");
	    }

	ptree dataRootContinuous = dataRoot.find("HYBRID_CONTINUOUS_SUBGRID")->second;
	if (dataRoot.find("HYBRID_DISCRETE_SUBGRID") == dataRoot.not_found())
	    {
	    spdlog::error("[ParametrManager] : Hybrid discrete subgrid is missing");
	    }
	ptree dataRootDiscrete = dataRoot.find("HYBRID_DISCRETE_SUBGRID")->second;
	int dim_hc, dim_hd;
	dim_hc = dataRootContinuous.get<int>("STATE_DIMENSION", 1);
	dim_hd = dataRootDiscrete.get<int>("STATE_DIMENSION", 1);

	dim = dim_hc + dim_hd;

	gridParameters->DIM = (unsigned long long int) dim;
	gridParameters->DIM_HC = (unsigned long long int) dim_hc;
	gridParameters->DIM_HD = (unsigned long long int) dim_hd;

	spdlog::info("[ParametrManager] : dim HC {}",dim_hc);
	spdlog::info("[ParametrManager] : dim HD {}",dim_hd);
	spdlog::info("[ParametrManager] : dim  {}",dim);
	//the main_dir parameter is taken from continuous grid
	gridParameters->GRID_MAIN_DIR = dataRootContinuous.get<int>("GRID_MAIN_DIR", 0);

	gridParameters->LIMINF = new double[dim];
	gridParameters->LIMSUP = new double[dim];
	gridParameters->NBPOINTS = new unsigned long long int[dim];
	gridParameters->PERIODIC = new bool[dim];
	gridParameters->SLICE_DIRS = new bool[dim];

	gridParameters->SLICE_VALUES = new double[dim];

	gridParameters->SLICE_VALUES_FD = new unsigned long long int[dim];

	gridParameters->SORTIE_OK_INF = new bool[dim];

	gridParameters->SORTIE_OK_SUP = new bool[dim];

	for (int k = 0; k < dim; k++)
	    {
	    gridParameters->LIMINF[k] = 0.0;
	    gridParameters->LIMSUP[k] = 1.0;
	    gridParameters->NBPOINTS[k] = 1;
	    gridParameters->PERIODIC[k] = false;
	    gridParameters->SLICE_DIRS[k] = 0;
	    gridParameters->SLICE_VALUES[k] = 0.0;
	    gridParameters->SLICE_VALUES_FD[k] = 0;
	    gridParameters->SORTIE_OK_INF[k] = false;
	    gridParameters->SORTIE_OK_SUP[k] = false;
	    }

	//only for continuous subgrid we need to take the lin inf and lim sup values. For discrete sub grid they are not used => can remain be default
	readTabData(&dataRootContinuous, gridParameters->LIMSUP, "STATE_MAX_VALUES", dim_hc, 1.0);
	readTabData(&dataRootContinuous, gridParameters->LIMINF, "STATE_MIN_VALUES", dim_hc, 0.0);

	readTabData(&dataRootContinuous, gridParameters->SORTIE_OK_SUP, "UPPER_LIMIT_IS_NOT_CONSTRAINT", dim_hc, false);
	readTabData(&dataRootContinuous, gridParameters->SORTIE_OK_INF, "LOWER_LIMIT_IS_NOT_CONSTRAINT", dim_hc, false);

	readTabData(&dataRootDiscrete, gridParameters->SLICE_VALUES_FD + dim_hc, "SLICE_LEVELS_DISCRETE", dim_hd, 0ULL);
	readTabData(&dataRootContinuous, gridParameters->SLICE_VALUES, "SLICE_LEVELS", dim_hc, 0.0);

	readTabData(&dataRootContinuous, gridParameters->SLICE_DIRS, "SLICE_DIRECTIONS", dim_hc, false);
	readTabData(&dataRootDiscrete, gridParameters->SLICE_DIRS + dim_hc, "SLICE_DIRECTIONS", dim_hd, false);

	readTabData(&dataRootContinuous, gridParameters->PERIODIC, "STATE_PERIODIC", dim_hc, false);

	readTabData(&dataRootContinuous, gridParameters->NBPOINTS, "STATE_GRID_POINTS", dim_hc, 2ULL);
	readTabData(&dataRootDiscrete, gridParameters->NBPOINTS + dim_hc, "STATE_GRID_POINTS", dim_hd, 2ULL);

	spdlog::info("[ParametrManager] : Read hybrid grid params : FINISHED");

	logVector("[ParametrManager] : Sortie sup hybrid: ", gridParameters->SORTIE_OK_SUP, dim);
	logVector("[ParametrManager] : Sortie inf hybrid: ", gridParameters->SORTIE_OK_INF, dim);
	logVector("[ParametrManager] : Lim sup hybrid: ", gridParameters->LIMSUP, dim);
	logVector("[ParametrManager] : Lim inf hybrid: ", gridParameters->LIMINF, dim);
	logVector("[ParametrManager] : NB Points hybrid: ", gridParameters->NBPOINTS, dim);
	logVector("[ParametrManager] : Periodic hybrid: ", gridParameters->PERIODIC, dim);
	logVector("[ParametrManager] : Slice dirs hybrid: ", gridParameters->SLICE_DIRS, dim);
	logVector("[ParametrManager] : Slice values continuous hybrid: ", gridParameters->SLICE_VALUES, dim);
	logVector("[ParametrManager] : Slice values discrete: ", gridParameters->SLICE_VALUES_FD, dim);
	}

    spdlog::info("[ParametrManager] : Read grid params : FINISHED");
    }

void ParametersManager::readSystemParametersFromJson(ptree &allParamsRoot)
    {
    string line;
    ostringstream os;

    string tempStr, tmpStr;
    ptree dataRoot;
    string sedCommand;

    spdlog::info("[ParametrManager] : Reading of System parameters");

    dataRoot = allParamsRoot.find("SYSTEM_PARAMETERS")->second;

    systemParameters->COMPUTE_LC = dataRoot.get < ComputeMethod > ("LIPSCHITZ_CONSTANT_COMPUTE_METHOD", ANALYTICAL_CALC);
    systemParameters->COMPUTE_MF = dataRoot.get < ComputeMethod > ("DYN_BOUND_COMPUTE_METHOD", ANALYTICAL_CALC);

    systemParameters->LIP = dataRoot.get<double>("LIPSCHITZ_CONSTANT", 1.0);
    systemParameters->MF = dataRoot.get<double>("DYN_BOUND", 1.0);

    spdlog::debug("[ParametrManager] : Dyn bounds estimation method {}, Dyn bounds default value  {}", toString(systemParameters->COMPUTE_MF),
	    systemParameters->MF);
    spdlog::debug("[ParametrManager] : Lipschitz constant estimation method {}, Lipschitz constant value {}", toString(systemParameters->COMPUTE_LC),
	    systemParameters->LIP);

    systemParameters->DYN_TYPE = dataRoot.get < DynType > ("DYNAMICS_TYPE", CC);

    systemParameters->globDeltat = dataRoot.get<bool>("IS_TIMESTEP_GLOBAL", false);

    systemParameters->SCHEME = dataRoot.get < TimeDiscretizationScheme > ("TIME_DISCRETIZATION_SCHEME", EL);
    systemParameters->SCALE_PARAM = dataRoot.get<bool>("SCALING", false);

    systemParameters->L_LIP = 1.0;
    systemParameters->L_MF = 1.0;

    if (gridParameters->GRID_METHOD != BS)
	{
	if (dataRoot.find("COST_LIPSCHITZ_CONSTANT") != dataRoot.not_found())
	    {
	    systemParameters->L_LIP = dataRoot.get<double>("COST_LIPSCHITZ_CONSTANT", 1.0);
	    }
	else if (algoParameters->COMPUTE_TMIN)
	    {
	    systemParameters->L_LIP = dataRoot.get<double>("TIME_HORIZON", 10.0);
	    }
	systemParameters->L_MF = dataRoot.get<double>("COST_BOUND_CONSTANT", 1.0);
	}

    spdlog::debug("[ParametrManager] : Cost function bound {}, Cost function Lipschitz contant {}", systemParameters->L_MF, systemParameters->L_LIP);
    spdlog::info("[ParametrManager] : Read system params : FINISHED");
    }

void ParametersManager::readAlgoParametersFromJson(ptree &allParamsRoot)
    {

    string line;
    ostringstream os;

    string tempStr, tmpStr;
    ptree dataRoot;
    string sedCommand;

    spdlog::info("[ParametrManager] : Reading of Algorithm parameters");

    dataRoot = allParamsRoot.find("ALGORITHM_PARAMETERS")->second;

    algoParameters->NB_OMP_THREADS = nbOmpThreads;

    algoParameters->FILE_PREFIX = gridParameters->FILE_PREFIX;

    algoParameters->COMPUTE_TMIN = dataRoot.get<bool>("COMPUTE_MIN_TIME", 0);
    algoParameters->GRID_REFINMENTS_NUMBER = dataRoot.get<int>("GRID_REFINMENTS_NUMBER", 0);
    algoParameters->COMPUTE_SET = dataRoot.get<bool>("COMPUTE_VIABLE_SET", true);
    algoParameters->SET_TYPE = dataRoot.get < SetType > ("SET_TYPE", VIAB);
    algoParameters->ITERATION_STOP_LEVEL = dataRoot.get<int>("ITERATION_STOP_LEVEL", 0);
    algoParameters->SAVE_SUBLEVEL = dataRoot.get<bool>("SAVE_SUBLEVEL", false);
    algoParameters->LEVEL = dataRoot.get<double>("LEVEL", 0.0);
    algoParameters->SAVE_VIAB_LIGHT = dataRoot.get<bool>("SAVE_VIABSET_LIGHT", false);

    if ((algoParameters->COMPUTE_SET == false) & (algoParameters->GRID_REFINMENTS_NUMBER > 0))
	{
	for (unsigned int k = 0; k < gridParameters->DIM; k++)
	    {
	    for (int j = 0; j < algoParameters->GRID_REFINMENTS_NUMBER; j++)
		{
		gridParameters->NBPOINTS[k] = 2 * gridParameters->NBPOINTS[k] - 1;
		}
	    }
	}

    algoParameters->INTERMEDIATE_SAVINGS = dataRoot.get<bool>("INTERMEDIATE_SAVINGS", 0);
    algoParameters->SAVE_BOUNDARY = dataRoot.get<bool>("SAVE_BOUNDARY", 0);
    algoParameters->SAVE_PROJECTION = dataRoot.get<bool>("SAVE_PROJECTION", 0);
    algoParameters->SAVE_SLICE = dataRoot.get<bool>("SAVE_SLICE", 0);
    algoParameters->SAVE_SLICE_BOUND = dataRoot.get<bool>("SAVE_SLICE_BOUND", 0);
    algoParameters->PROJECTION = new unsigned long long int[gridParameters->DIM];
    this->readTabData(&dataRoot, algoParameters->PROJECTION, "PROJECTION_AXIS", gridParameters->DIM, 0ULL);
    algoParameters->TARGET_OR_DEPARTURE = dataRoot.get < TargetOrDeparture > ("TARGET_OR_DEPARTURE", DEPARTURE);

    spdlog::info("[ParametrManager] : Read algorithm params : FINISHED");
    }

void ParametersManager::readControlParametersFromJson(ptree &allParamsRoot)
    {

    string line;
    ostringstream os;

    string tempStr, tmpStr;
    ptree dataRoot;
    string sedCommand;

    spdlog::info("[ParametrManager] : Reading of controls parameters");

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
	this->readTabData(&dataRoot, controlParameters->LIMSUPC, "CONTROL_MAX_VALUES", dimc, 1.0);
	this->readTabData(&dataRoot, controlParameters->LIMINFC, "CONTROL_MIN_VALUES", dimc, 0.0);

	controlParameters->NBPOINTSC = new unsigned long long int[dimc];
	this->readTabData(&dataRoot, controlParameters->NBPOINTSC, "CONTROL_GRID_POINTS", dimc, 1ULL);
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
	this->readTabData(&dataRoot, controlParameters->LIMSUP_TY, "CONTROL_TY_MAX_VALUES", dimc_ty, 0.0);
	this->readTabData(&dataRoot, controlParameters->LIMINF_TY, "CONTROL_TY_MIN_VALUES", dimc_ty, 1.0);

	controlParameters->NBPOINTS_TY = new unsigned long long int[dimc_ty];
	this->readTabData(&dataRoot, controlParameters->NBPOINTS_TY, "CONTROL_TY_GRID_POINTS", dimc_ty, 1ULL);
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

    int dimc_ht;
    if (dataRoot.find("CONTROL_HYBRID_TRANSITION_DIMENSION") != dataRoot.not_found())
	{
	dimc_ht = dataRoot.get<int>("CONTROL_HYBRID_TRANSITION_DIMENSION", 1);
	controlParameters->DIM_HT = (unsigned long long int) dimc_ht;
	/*
	 * initializing with empty tabs, default values, not used in this case
	 */
	controlParameters->LIMINF_HT = new double[dimc_ht];
	controlParameters->LIMSUP_HT = new double[dimc_ht];
	this->readTabData(&dataRoot, controlParameters->LIMSUP_HT, "CONTROL_HT_MAX_VALUES", dimc_ht, 0.0);
	this->readTabData(&dataRoot, controlParameters->LIMINF_HT, "CONTROL_HT_MIN_VALUES", dimc_ht, 1.0);

	controlParameters->NBPOINTS_HT = new unsigned long long int[dimc_ht];
	this->readTabData(&dataRoot, controlParameters->NBPOINTS_HT, "CONTROL_HT_GRID_POINTS", dimc_ht, 1ULL);
	}
    else
	{
	dimc_ht = 0;
	controlParameters->DIM_TY = 0;
	/*
	 * initializing with empty tabs, default values, not used in this case
	 */
	controlParameters->LIMINF_HT = new double[dimc_ht];
	controlParameters->LIMSUP_HT = new double[dimc_ht];
	controlParameters->NBPOINTS_HT = new unsigned long long int[dimc_ht];
	}

    spdlog::info("[ParametrManager] : Read control : FINISHED");
    }

void ParametersManager::readTrajectoryParametersListFromJson(ptree &allParamsRoot)
    {
    ptree dataRoot;

    spdlog::info("[ParametrManager] : Reading of Trajectory parameters");

    ptree::assoc_iterator trajParamsPos = allParamsRoot.find("TRAJECTORY_PARAMETERS");
    if (trajParamsPos == allParamsRoot.not_found())
	{
	spdlog::info("[ParametrManager] : No trajectory parameters found");
	// Pour s'assurer que le vecteur est vide
	trajectoryParametersList.clear();
	return;
	}

    ptree defaultValues;
    ptree::assoc_iterator defaultValuesPos;
    dataRoot = trajParamsPos->second;
    if (dataRoot.empty())
	{
	spdlog::warn(
		"[ParametrManager] : Trajectory parameters should be an array or an objet with fields DEFAULT_VALUES and TRAJECTORY_SPECIFIC_VALUES. No trajectory calculations");
	}
    else if ((defaultValuesPos = dataRoot.find("DEFAULT_VALUES")) != dataRoot.not_found())
	{
	// On suppose avoir des paramètres de trajectoire avec valeur par défaut
	defaultValues = defaultValuesPos->second;
	ptree::assoc_iterator newDataRootPos = dataRoot.find("TRAJECTORY_SPECIFIC_VALUES");
	if (newDataRootPos == dataRoot.not_found())
	    {
	    spdlog::error("[ParametrManager] : Default values given without trajectory specific values. Ill formed trajectory parameters");
	    spdlog::error("[ParametrManager] : No trajectory calculation");
	    trajectoryParametersList.clear();
	    return;
	    }
	else
	    {
	    dataRoot = newDataRootPos->second;
	    }
	}

    for (ptree::iterator pos = dataRoot.begin(); pos != dataRoot.end(); ++pos)
	{
	trajectoryParams trajParams;
	mergeJSONPtreeInto(defaultValues, pos->second);
	readTrajectoryParameters(pos->second, trajParams);
	trajectoryParametersList.push_back(std::move(trajParams));
	}
    spdlog::info("[ParametrManager] : Read trajectory parameters : FINISHED");
    }

void ParametersManager::readTrajectoryParameters(ptree &trajectoryParamsRoot, trajectoryParams &trajParams)
    {

    ptree &dataRoot = trajectoryParamsRoot;

    trajParams.maxTime = dataRoot.get<double>("TRAJECTORY_TIME_HORIZON", 10.0);

    if (dataRoot.find("TRAJECTORY_TYPE") != dataRoot.not_found())
	{
	trajParams.NB_STRATEGIES = dataRoot.find("TRAJECTORY_TYPE")->second.size();
	}
    else
	{
	trajParams.NB_STRATEGIES = 0;
	}

    if (trajParams.NB_STRATEGIES > 0)
	{
	trajParams.TRAJECTORY_TYPE = STRATEGY_LIST;
	trajParams.STRATEGIES = new ControlPickStrategyName[trajParams.NB_STRATEGIES];
	readTabDataSkipInvalid(&dataRoot, trajParams.STRATEGIES, "TRAJECTORY_TYPE", trajParams.NB_STRATEGIES);
	initBubble(dataRoot, trajParams);
	}
    else
	{
	trajParams.TRAJECTORY_TYPE = dataRoot.get < TypeTraj > ("TRAJECTORY_TYPE", VD);
	if (isCautious(trajParams.TRAJECTORY_TYPE))
	    {
	    initBubble(dataRoot, trajParams);
	    }
	}

    double *initPoint = new double[gridParameters->DIM];
    readTabData(&dataRoot, initPoint, "INITIAL_POINT", gridParameters->DIM, 0.0);
    trajParams.INIT_POINT = initPoint;

    unsigned long long int *initPoint_fd = new unsigned long long int[gridParameters->DIM];
    readTabData(&dataRoot, initPoint_fd, "INITIAL_POINT_DISCRETE", gridParameters->DIM, 0ULL);

    trajParams.INIT_POINT_FD = initPoint_fd;

    trajParams.INIT_VALUE = dataRoot.get<double>("INITIAL_VALUE", 0.0);
    trajParams.INIT_VALUE_FD = dataRoot.get<double>("INITIAL_VALUE_DISCRETE", 0.0);

    double *initControl = new double[controlParameters->DIMC];
    readTabData(&dataRoot, initControl, "INITIAL_CONTROL", controlParameters->DIMC, 0.0);

    for (unsigned long long int d = 0; d < controlParameters->DIMC; ++d)
	{
	if (initControl[d] < controlParameters->LIMINFC[d])
	    {
	    spdlog::error("Initial control n°{} outside of control domain.", d + 1);
	    spdlog::warn("Clamping coordinate of value {} (index n°{}) into [{} ,{}]", trajectoryParametersList.size() + 1, initControl[d], d + 1,
		    controlParameters->LIMINFC[d], controlParameters->LIMSUPC[d]);
	    initControl[d] = controlParameters->LIMINFC[d];
	    }
	else if (initControl[d] > controlParameters->LIMSUPC[d])
	    {
	    spdlog::error("Initial control n°{} outside of control domain.");
	    spdlog::warn("Clamping coordinate of value {} (index n°{}) into [{} ,{}]", trajectoryParametersList.size() + 1, initControl[d], d + 1,
		    controlParameters->LIMINFC[d], controlParameters->LIMSUPC[d]);
	    initControl[d] = controlParameters->LIMSUPC[d];
	    }
	}
    trajParams.INIT_CONTROL = initControl;

    trajParams.REAL_TIME_STEPS_PER_DISCRETE_STEP = dataRoot.get<int>("REAL_TIME_STEPS_PER_DISCRETE_STEP", 1);

    boost::optional<double> optAngle = dataRoot.get_optional<double>("MAX_ANGLE_RADIANS");
    if (optAngle)
	{
	trajParams.MAX_ANGLE_RADIANS = (*optAngle);
	}
    else if ((optAngle = dataRoot.get_optional<double>("MAX_ANGLE_DEGREES")))
	{
	trajParams.MAX_ANGLE_RADIANS = (*optAngle) * pi / 180;
	}
    else
	{
	trajParams.MAX_ANGLE_RADIANS = pi / 2;
	}

    trajParams.SAVE_PICKING_STRATEGY = dataRoot.get<bool>("SAVE_PICKING_STRATEGY", true);
    trajParams.ARE_STRATEGIES_GUARANTEED = dataRoot.get<bool>("ARE_STRATEGIES_GUARANTEED", true);
    initSeed(dataRoot, trajParams);

    if (controlParameters->DIM_TY > 0)
	{
	ptree tycheDataRoot;
	ptree defaultValues;

	ptree::assoc_iterator tycheParamsPos = dataRoot.find("TYCHE_PARAMETERS");
	if (tycheParamsPos == dataRoot.not_found())
	    {
	    spdlog::warn("[ParametrManager] : No tyche parameters found for tychastic problem, all values will be the default");
	    }
	else
	    {
	    trajParams.TYCHE_PARAMS = new tycheParams[controlParameters->DIM_TY];
	    tycheDataRoot = tycheParamsPos->second;
	    ptree::assoc_iterator defaultValuesPos = tycheDataRoot.find("DEFAULT_VALUES");
	    if (defaultValuesPos != tycheDataRoot.not_found())
		{
		defaultValues = defaultValuesPos->second;
		ptree::assoc_iterator newDataRootPos = tycheDataRoot.find("TYCHE_SPECIFIC_VALUES");
		if (newDataRootPos == tycheDataRoot.not_found())
		    {
		    tycheDataRoot.clear();
		    }
		else
		    {
		    tycheDataRoot = newDataRootPos->second;
		    }
		}
	    }
	if (tycheDataRoot.size() < controlParameters->DIM_TY)
	    {
	    spdlog::warn("[ParametrManager] : Tyche parameters list size (which is {}) smaller than CONTROL_TY_DIMENSION (which is {})",
		    tycheDataRoot.size(), controlParameters->DIM_TY);
	    spdlog::warn(
		    "[ParametrManager] : The parameters list will be consumed until no more parameters are left, after which default values will be given");
	    }
	else if (tycheDataRoot.size() > controlParameters->DIM_TY)
	    {
	    spdlog::warn("[ParametrManager] : Tyche parameters list size (which is {}) different to CONTROL_TY_DIMENSION (which is {})",
		    tycheDataRoot.size(), controlParameters->DIM_TY);
	    spdlog::warn("[ParametrManager] : Parameters whose index exceeds the dimension will be ignored");
	    tycheDataRoot.erase(std::next(tycheDataRoot.begin(), controlParameters->DIM_TY), tycheDataRoot.end());
	    }

	ptree::iterator it = tycheDataRoot.begin();
	for (unsigned long long int i = 0; i < tycheDataRoot.size(); ++i)
	    {
	    mergeJSONPtreeInto(defaultValues, it->second);
	    readTycheParameters(it->second, trajParams.TYCHE_PARAMS[i]);
	    ++it;
	    }
	for (unsigned long long int i = tycheDataRoot.size(); i < controlParameters->DIM_TY; ++i)
	    {
	    readTycheParameters(defaultValues, trajParams.TYCHE_PARAMS[i]);
	    }
	}
    }

void ParametersManager::readTycheParameters(ptree &dataRoot, tycheParams &tycheParams)
    {

    tycheParams.TYCHE_DISTRIBUTION = dataRoot.get < TycheDistribution > ("TYCHE_DISTRIBUTION", UNIFORM);
    boost::optional<double> optTycheValue = dataRoot.get_optional<double>("CONSTANT_TYCHE_VALUE");

    if (tycheParams.TYCHE_DISTRIBUTION == CONSTANT)
	{
	if (optTycheValue)
	    {
	    tycheParams.CONSTANT_TYCHE_VALUE = *optTycheValue;
	    }
	else
	    {
	    spdlog::warn("Constant tyche distribution requested but no value was given, defaulting to 0");
	    tycheParams.CONSTANT_TYCHE_VALUE = 0.0;
	    }
	}

    tycheParams.MAX_NB_REROLLS = dataRoot.get<int>("MAX_NB_REROLLS", 10);
    }

void ParametersManager::readModelParametersFromJson(ptree &allParamsRoot)
    {
    ptree dataRoot;

    spdlog::info("[ParametrManager] : Reading of Model parameters");

    ptree::assoc_iterator modelParametersPos = allParamsRoot.find("MODEL_PARAMETERS");
    if (modelParametersPos != allParamsRoot.not_found())
	{
	dataRoot = modelParametersPos->second;

	for (ptree::iterator pos = dataRoot.begin(); pos != dataRoot.end(); ++pos)
	    {
	    if (pos->second.empty())
		{
		modelParameters[pos->first] = pos->second.data();
		}
	    else
		{
		for (ptree::iterator it = pos->second.begin(); it != pos->second.end(); ++it)
		    {
		    modelParameters.addToList(pos->first, it->second.data());
		    }
		}
	    }
	}
    else
	{
	spdlog::info("[ParametrManager] : No model parameters found");
	}
    spdlog::info("[ParametrManager] : Read model parameters : FINISHED");
    }

void ParametersManager::toPixels(double *bubbleRadius)
    {
    for (unsigned long long int i = 0; i < gridParameters->DIM; ++i)
	{
	double unitPerPixels = (gridParameters->LIMSUP[i] - gridParameters->LIMINF[i]) / gridParameters->NBPOINTS[i];
	bubbleRadius[i] = bubbleRadius[i] / unitPerPixels;
	}
    }

void ParametersManager::initBubble(ptree &dataRoot, trajectoryParams &trajParameters)
    {

    trajParameters.BUBBLE_RADIUS = new double[gridParameters->DIM];
    double *const bubbleRadius = trajParameters.BUBBLE_RADIUS;

    try
    {
	this->readTabData(&dataRoot, bubbleRadius, "BUBBLE_RADIUS", gridParameters->DIM, 0.5);
    }
    // Si la taille n'est pas donnée comme un tableau, peut-être qu'elle est donnée comme un nombre ?
    catch (...)
	{
	double radius = dataRoot.get<double>("BUBBLE_RADIUS", 0.5);
	std::fill(bubbleRadius, bubbleRadius + gridParameters->DIM, radius);
	}

    trajParameters.BUBBLE_INTERPRETATION = dataRoot.get < BubbleInterpretation > ("BUBBLE_INTERPRETATION", MOORE);

    bool allEqual = true;

    if (!isInPixelsUnits(trajParameters.BUBBLE_INTERPRETATION))
	{
	toPixels(bubbleRadius);
	}

    bubbleRadius[0] = std::ceil(bubbleRadius[0]);
    for (unsigned long long int i = 1; i < gridParameters->DIM; ++i)
	{
	// Notre code n'exploite qu'une bulle en taille de pixels entiers
	// On veut que la bulle englobe entièrement la zone souhaitée
	// et la bulle doit au moins avoir une taille en pixels de 1, d'où le ceil
	bubbleRadius[i] = std::ceil(bubbleRadius[i]);
	allEqual &= (bubbleRadius[i - 1] == bubbleRadius[i]);
	}

    if (!allEqual)
	{
	if (trajParameters.BUBBLE_INTERPRETATION == EUCLIDEAN_PX)
	    {
	    spdlog::warn(
		    "[ParametrManager] : Euclidean pixel distance requested but different bubble radiuses given for each dimension. Assuming ELLIPTIC_PX");
	    trajParameters.BUBBLE_INTERPRETATION = ELLIPTIC_PX;
	    }
	else if (trajParameters.BUBBLE_INTERPRETATION == EUCLIDEAN)
	    {
	    trajParameters.BUBBLE_INTERPRETATION = ELLIPTIC;
	    }
	}
    else
	{
	if (trajParameters.BUBBLE_INTERPRETATION == ELLIPTIC)
	    {
	    trajParameters.BUBBLE_INTERPRETATION = EUCLIDEAN;
	    }
	else if (trajParameters.BUBBLE_INTERPRETATION == ELLIPTIC_PX)
	    {
	    trajParameters.BUBBLE_INTERPRETATION = EUCLIDEAN_PX;
	    }
	}
    }

void ParametersManager::initSeed(ptree &dataRoot, trajectoryParams &params)
    {
    ptree::assoc_iterator seedRoot = dataRoot.find("SEED");
    if (seedRoot == dataRoot.not_found())
	{
	params.SEED_LENGTH = 0;
	}
    else if (seedRoot->second.empty())
	{
	boost::optional<int> seedValue = dataRoot.get_optional<int>("SEED");
	if (!seedValue)
	    {
	    params.SEED_LENGTH = 0;
	    }
	else
	    {
	    params.SEED_LENGTH = 1;
	    params.SEED = new unsigned long long int[1];
	    params.SEED[0] = *seedValue;
	    }

	}
    else
	{
	params.SEED_LENGTH = seedRoot->second.size();
	params.SEED = new unsigned long long int[params.SEED_LENGTH];
	readTabDataSkipInvalid(&dataRoot, params.SEED, "SEED", params.SEED_LENGTH);
	}
    }

void mergeJSONPtreeInto(const ptree &from, ptree &to)
    {
    int i = 0;
    for (ptree::const_iterator it = from.begin(); it != from.end(); ++it)
	{
	ptree::assoc_iterator child = to.find(it->first);
	if (child == to.not_found())
	    {
	    to.put_child(it->first, it->second);
	    }
	// Si le champ est un champ composé
	else if (!it->second.empty())
	    {
	    // Si on a un tableau
	    if (it->first == "")
		{
		// Les deux tableaux doivent avoir la même taille
		for (auto j = to.size(); j < from.size(); ++j)
		    {
		    ptree newEmpty;
		    to.put_child("", newEmpty);
		    }
		mergeJSONPtreeInto(it->second, std::next(to.begin(), i)->second);
		}
	    // Si on a un objet
	    else
		{
		mergeJSONPtreeInto(it->second, child->second);
		}
	    }
	++i;
	}
    }

ParametersManager::~ParametersManager()
    {
    // TODO Auto-generated destructor stub
    }

TrajectoryParametersManager::TrajectoryParametersManager(ParametersManager *pm, int trajIndex) :
		parametersManager(pm), trajIndex(trajIndex)
    {
    }

const trajectoryParams* TrajectoryParametersManager::getTrajectoryParameters() const
    {
    return &(parametersManager->getTrajectoryParametersList()[trajIndex]);
    }

int TrajectoryParametersManager::getNbTrajectories() const
    {
    return parametersManager->getNbTrajectories();
    }

int TrajectoryParametersManager::getTrajectoryIndex() const
    {
    return trajIndex;
    }

const gridParams* TrajectoryParametersManager::getGridParameters() const
    {
    return parametersManager->getGridParameters();
    }

const algoViabiParams* TrajectoryParametersManager::getAlgoParameters() const
    {
    return parametersManager->getAlgoParameters();
    }

const controlParams* TrajectoryParametersManager::getControlParameters() const
    {
    return parametersManager->getControlParameters();
    }

const systemParams* TrajectoryParametersManager::getSystemParameters() const
    {
    return parametersManager->getSystemParameters();
    }

void* TrajectoryParametersManager::getModelHandle(void)
    {
    return parametersManager->getModelHandle();
    }

const modelParams* TrajectoryParametersManager::getModelParameters() const
    {
    return parametersManager->getModelParameters();
    }
