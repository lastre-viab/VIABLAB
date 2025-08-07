/*
 * main.cpp
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
 *  Created on: reated on: 27 jul. 2013
 *      Author: Anna DESILLES
 */

#include "../include/defs.h"
#include "../include/ViabProblemFactory.h"
#include "../include/DefaultValues.h"

#include "../src/initParams.h"

#include <dlfcn.h>

int getNbOmpThreads(int argc, char **argv) {
    
    int nbOmpThreads = 1;
    if (argc >= 5)
	{
        if (strcmp(argv[3], "-t") == 0)
	    {
            nbOmpThreads = atoi(argv[4]);
	    }
        else
	    {
            printf("Incorrect arguments. Ignored\n");
	    }
	}
    return nbOmpThreads;
}

// À possiblement remplacer par la version de la librairie standard C++20
bool ends_with(const std::string &str, const std::string &suffix) {

    int sufSize = suffix.size();
    int strSize = str.size();
    
    return strSize >= sufSize && (str.compare(strSize - sufSize, sufSize, suffix) == 0);
}

// À possiblement remplacer par la version de la librairie standard C++20
bool starts_with(const std::string &str, const std::string &prefix) {

    int prefSize = prefix.size();
    int strSize = str.size();
    
    return prefSize <= strSize && (str.compare(0, prefSize, prefix) == 0);
}

void *getModelLibrary(int argc, char **argv) {
    
    if (argc < 2) {
        spdlog::error("No model loaded as no library name was given as an argument");
        return nullptr;
    }

    spdlog::info("Loading model information");
    
    string libLocation(argv[1]);
    ostringstream os;

    os << "./" << libLocation;
    
    if (!ends_with(libLocation, ".so")) {
        os << ".so";
    }    

    void *libHandle = dlopen(os.str().c_str(), RTLD_NOW);
    if (libHandle == nullptr) {
        spdlog::error("Could not load model {0}. Have you created a file named data/{0}.cpp containing your model information?", argv[1]);
        spdlog::error("Error message : {}", dlerror());
    }
    return libHandle;
}

void initLoggers(const string modelFilePrefix) {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");
    auto str = oss.str();
    string logFileName = "../LOG/" + modelFilePrefix + str + ".log";
    
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logFileName, true);
    file_sink->set_level(spdlog::level::trace);

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::info);
    console_sink->set_pattern("[ViabLog] [%^%l%$] %v");

    // Liste des sinks réellement utilisés
    spdlog::sinks_init_list sink_list =
        {
            file_sink, console_sink
        };

    //spdlog::logger logger("multi_sink", sink_list.begin(), sink_list.end());

    // or you can even set multi_sink logger as default logger
    spdlog::set_default_logger(std::make_shared<spdlog::logger>("ViabLog", sink_list));

    spdlog::set_level(spdlog::level::trace);

    //spdlog::flush_every(std::chrono::seconds(30));
}

std::string getParamsFileName(int argc, char **argv, void *modelHandle) {
    std::string paramsFile;
    
    if (argc > 2) {
        paramsFile = argv[2];
    }
    else {
        string *paramsFilePtr = (string *) dlsym(modelHandle, "paramsFile");
        if (paramsFilePtr == nullptr) {
            spdlog::error("No parameters file name given in model or as program argument. Initialize a paramsFile variable in your model or give a valid input file name as an argument.");
            dlclose(modelHandle);
            std::exit(1);
        }
        else {
            spdlog::info("Parameters file is {}. Loaded from model file.", *paramsFilePtr);
        }
        paramsFile = *paramsFilePtr;
    }
    
    return paramsFile;
}

int main(int argc, char **argv)
{    
    int nbOmpThreads = getNbOmpThreads(argc, argv);

    ostringstream os;
    string fileName;
    double t1, t2, elapsed_time, t1_glob, t2_glob, elapsed_time_glob;
    timeval tim, tim_glob;

    setvbuf(stdout, NULL, _IONBF, 0);
    char cCurrentPath[FILENAME_MAX];
    cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

    void *modelHandle = getModelLibrary(argc, argv);
    if (modelHandle == nullptr) {
        return 1;
    }
    
    const std::string &paramsFile = getParamsFileName(argc, argv, modelHandle);
    
    ptree allParamsRoot, dataRoot;
    string input_tempfile = "../INPUT/" + paramsFile;
    //- loading data from input file
    read_json(input_tempfile, allParamsRoot);
    dataRoot = allParamsRoot.find("GRID_PARAMETERS")->second;

    string modelFilePrefix = dataRoot.get<string>("OUTPUT_FILE_PREFIX", "Model-");

    initLoggers(modelFilePrefix);
    
    //initialisation des structures à partir des données de l'utilisateur
    spdlog::info("Start of Viablab programm : initialization of model parameters");
    ParametersManager *PM = initParams(modelHandle, paramsFile, gp, avp, cp, sp, nbOmpThreads);
    /**********************************
     * IMPORTANT : appeler la fonction load_data ICI avant TOUT le reste
     */
    spdlog::info("Starting loading model specific data");
    void (*loadModel)(const ParametersManager *) = getUserSymbolOrDefault(modelHandle, loadModelData);
    loadModel(PM);
	

    /*==========================================================*/

    Viabi *viabProblem;

    ViabProblemFactory vpf(PM);

    viabProblem = vpf.constructViabilityProblem(gp.GRID_METHOD);

    if (avp.COMPUTE_SET)
	{
        // l'ensemble doit être  recalculé
        gettimeofday(&tim_glob, NULL);             //mesure le temps d'execution
        t1_glob = tim_glob.tv_sec + (tim_glob.tv_usec / 1000000.0);

        if (avp.SET_TYPE == VIAB)
	    { // option de calcul de noyau de viabilité a été choisie
	      // Initialisation de l'ensemble de contraintes
            spdlog::warn("Initialization of constraints");
            viabProblem->initialiseConstraints();
            viabProblem->ViabilityKernel(avp.ITERATION_STOP_LEVEL);
	    }
        if (avp.SET_TYPE == CAPT)
	    { // option de calcul de noyau de viabilité a été choisie
	      // Initialisation de l'ensemble de contraintes
            spdlog::warn("Initialization of target set");
            viabProblem->initialiseTarget();
            try
            {
                cout << " debut calcul capture \n";
                viabProblem->CaptureBasin();
            }
            catch (...)
            {
                cout << " exception while capture baisin\n" << std::flush;

                exit(1);
            }

	    }
        if (avp.SET_TYPE == VIABG)
	    { // option de calcul de noyau de viabilité a été choisie
	      // Initialisation de l'ensemble de contraintes
            spdlog::warn("Initialization of constraints");
            viabProblem->initialiseConstraints();
            viabProblem->GarantedViabilityKernel(avp.ITERATION_STOP_LEVEL);
	    }

        gettimeofday(&tim_glob, NULL);
        t2_glob = tim_glob.tv_sec + (tim_glob.tv_usec / 1000000.0);
        elapsed_time_glob = (double) ((t2_glob - t1_glob));
        spdlog::warn("Elapsed time  GLOBAL: {} sec", elapsed_time_glob);

	}			//fin de if(computeSet)
    else
	{ // on charge dans la mémoire l'ensemble calculé et enregistré
	  // correspondant au dernier raffinement
        viabProblem->loadViableSets();
	}

    /*
     * Si les trajectoires doivent être calculées
     */
    
    viabProblem->computeTrajectories();
    if (avp.SAVE_PROJECTION | avp.SAVE_SLICE_BOUND)
	{
        viabProblem->saveViableSets();
	}
    delete viabProblem;

    void (*postProcess_)(const ParametersManager *) = getUserSymbolOrDefault(modelHandle, postProcess);
    // À définir par l'utilisateur
    postProcess_(PM);

    dlclose(modelHandle);
}
