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

#include "../data/ModelDataInclusion.h"

#include "../src/initParams.h"

int main(int argc, char **argv)
    {

    int nbOmpThreads = 1;
    if (argc >= 3)
	{
	if (strcmp(argv[1], "-t") == 0)
	    {
	    nbOmpThreads = atoi(argv[2]);
	    }
	else
	    {
	    printf("Incorrect arguments. Ignored\n");
	    }
	}

    ostringstream os;
    string fileName;
    double t1, t2, elapsed_time, t1_glob, t2_glob, elapsed_time_glob;
    timeval tim, tim_glob;

    setvbuf(stdout, NULL, _IONBF, 0);
    char cCurrentPath[FILENAME_MAX];
    cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

    ptree allParamsRoot, dataRoot;
    string input_tempfile = "../INPUT/" + paramsFile;
    //- loading data from input file
    read_json(input_tempfile, allParamsRoot);
    dataRoot = allParamsRoot.find("GRID_PARAMETERS")->second;

    string modelFilePrefix = dataRoot.get<string>("OUTPUT_FILE_PREFIX", "Model-");

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

//	spdlog::sinks_init_list sink_list =
//	{
//			file_sink, console_sink
//	};
    spdlog::sinks_init_list sink_list =
	{
	console_sink
	};

    //spdlog::logger logger("multi_sink", sink_list.begin(), sink_list.end());

    // or you can even set multi_sink logger as default logger
    spdlog::set_default_logger(std::make_shared<spdlog::logger>("ViabLog", sink_list));

    spdlog::set_level(spdlog::level::trace);

    //spdlog::flush_every(std::chrono::seconds(30));
    //initialisation des structures à partir des données de l'utilisateur
    spdlog::info("Start of Viablab programm : initialization of model parameters");
    ParametersManager *PM = initParams(gp, avp, cp, sp, nbOmpThreads);
    /**********************************
     * IMPORTANT : appeler la fonction load_data ICI avant TOUT le reste
     */
    spdlog::info("Starting loading model specific data");
    if (loadModelData)
	{
	loadModelData(PM);
	}

    /*==========================================================*/

    Viabi *viabProblem;

    ViabProblemFactory *vpf = new ViabProblemFactory(PM);

    viabProblem = vpf->constructViabilityProblem(gp.GRID_METHOD);

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
	    viabProblem->ViabilityKernel(true, avp.INTERATION_STOP_LEVEL);
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
	    viabProblem->GarantedViabilityKernel(true, avp.INTERATION_STOP_LEVEL);
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

    postProcess(PM);

    }
