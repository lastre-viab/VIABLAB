/*
 * main.cpp
 *
 *  Created on: 10 sept. 2013
 *      Author: ANYA
 */
/*
 * Inclusions des classes utilisées pour le calcul
 */
//#include "../include/ViabiBitSet.h"
//#include "../include/Viabi.h"

#include "../include/defs.h"
#include "../include/ViabProblemFactory.h"
#include "../data/DefaultValues.h"



#include "../data/ModelDataInclusion.h"


#include "../src/initParams.h"


int main( int argc, char** argv){

	int nbOmpThreads = 1;
	if(argc>=3){
		if(strcmp(argv[1], "-t") == 0){
			nbOmpThreads = atoi(argv[2]);
		}
		else{
			printf("Incorrect arguments. Ignored\n");
		}
	}

	cout<< "main nb threads :"<<nbOmpThreads<<" \n";
	ostringstream os;
	string fileName;
	double t1,t2,elapsed_time, t1_glob, t2_glob, elapsed_time_glob;
	timeval tim, tim_glob;

	setvbuf(stdout, NULL, _IONBF, 0);
	char cCurrentPath[FILENAME_MAX];
	cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

	//initialisation des structures à partir des données de l'utilisateur

	ParametersManager *PM= initParams(gp, avp, cp, sp, nbOmpThreads);
	/**********************************
	 * IMPORTANT : appeler la fonction load_data ICI avant TOUT le reste
	 */

	cout<< " load model data commence \n";
	if(loadModelData)
	{
		loadModelData( PM);
	}

	/*==========================================================*/

	Viabi *viabProblem;

	ViabProblemFactory * vpf=new ViabProblemFactory(PM);

	viabProblem=vpf->constructViabilityProblem(gp.GRID_METHOD);

	if(avp.COMPUTE_SET){
		// l'ensemble doit être  recalculé
		gettimeofday(&tim_glob,NULL);              //mesure le temps d'execution
		t1_glob=tim_glob.tv_sec+(tim_glob.tv_usec/1000000.0);

		if(avp.SET_TYPE==VIAB) { // option de calcul de noyau de viabilité a été choisie
			// Initialisation de l'ensemble de contraintes
			cout<< " initialise constr \n";
			viabProblem->initialiseConstraints();
			viabProblem->ViabilityKernel( true, avp.INTERATION_STOP_LEVEL);
		}
		if(avp.SET_TYPE==CAPT) { // option de calcul de noyau de viabilité a été choisie
			// Initialisation de l'ensemble de contraintes
			cout<< " initialise target  \n";
			viabProblem->initialiseTarget();
			viabProblem->CaptureBasin();
		}
		if(avp.SET_TYPE==VIABG) { // option de calcul de noyau de viabilité a été choisie
			// Initialisation de l'ensemble de contraintes
			cout<< " initialise constr \n";
			viabProblem->initialiseConstraints();
			viabProblem->GarantedViabilityKernel( true, avp.INTERATION_STOP_LEVEL);
		}

		cout<< "=============================================================================="<<endl;
		gettimeofday(&tim_glob,NULL);
		t2_glob=tim_glob.tv_sec+(tim_glob.tv_usec/1000000.0);
		elapsed_time_glob=(double)((t2_glob-t1_glob));
		cout << "Elapsed time  GLOBAL: " << elapsed_time_glob << " sec." << endl << endl;
		cout<< "=============================================================================="<<endl;
	}//fin de if(computeSet)
	else
	{ // on charge dans la mémoire l'ensemble calculé et enregistré
		// correspondant au dernier raffinement
		viabProblem->loadViableSets();
	}

	/*
	 * Si les trajectoires doivent être calculées
	 */

	viabProblem->computeTrajectories();
	if(avp.SAVE_PROJECTION | avp.SAVE_SLICE_BOUND)
	{
		viabProblem->saveViableSets();
	}
	delete viabProblem;
	cout<< " viab detruit retour main\n";
	postProcess(PM);

}
