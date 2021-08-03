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



/*
 * inclusion des déclarations du modèle
 */
//#include  "../data/equilibres4D_data.h"
//#include  "../data/ex1Viabi2D_data.h"

#include  "../data/julia2D_data.h"
//#include  "../data/ex4_Viabi2D_data.h"
//#include  "../data/ex3Dim_Viabi2D_data.h"
//#include  "../data/ex1_Viabi_multiDim_data.h"
//#include  "../data/ex3bis_Viabi2D_data.h"
//#include "../data/ex2_EcoPolut_data.h"

//#include "../data/userModel.h"
//#include "../data/test_capt_discret_data.h"
//#include "../data/corentin_no_control_feed_synth_data_funct_traj_v3.h"

//#include "../data/corentin_no_control_feed_synth_data_funct_traj_v3.h"
/*
 * Code de la fonction d'initialisation de paramètres
 */
#include "../src/initParams.h"


int main( int argc, char** argv){

	if(argc>=3){
		if(strcmp(argv[1], "-t") == 0){
			ompThreads = atoi(argv[2]);
		}
		else{
			printf("Incorrect arguments. Ignored\n");
		}
	}


	ostringstream os;
	string fileName;
	double t1,t2,elapsed_time, t1_glob, t2_glob, elapsed_time_glob;
	timeval tim, tim_glob;

	setvbuf(stdout, NULL, _IONBF, 0);
	char cCurrentPath[FILENAME_MAX];
	cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

	/**********************************
	 * IMPORTANT : appeler la fonction load_data ICI avant TOUT le reste
	 */
	loadModelData();
	/*==========================================================*/



	Viabi *viabProblem;
	// déclaration de structures pour les paramètres
	gridParams gp;
	algoViabiParams avp;
	controlParams cp;
	systemParams sp;
	//initialisation des structures à partir des données de l'utilisateur
	ParametersManager *PM= initParams(gp, avp, cp, sp);
	ViabProblemFactory * vpf=new ViabProblemFactory(PM);

	viabProblem=vpf->constructViabilityProblem(gridMethod);

	if(computeSet){
		// l'ensemble doit être  recalculé
		gettimeofday(&tim_glob,NULL);              //mesure le temps d'execution
		t1_glob=tim_glob.tv_sec+(tim_glob.tv_usec/1000000.0);

		if(setType==VIAB) { // option de calcul de noyau de viabilité a été choisie
			// Initialisation de l'ensemble de contraintes
			cout<< " initialise constr \n";
			viabProblem->initialiseConstraints();
			viabProblem->ViabilityKernel( true, 1);
		}
		if(setType==CAPT) { // option de calcul de noyau de viabilité a été choisie
			// Initialisation de l'ensemble de contraintes
			cout<< " initialise target  \n";
			viabProblem->initialiseTarget();
			viabProblem->CaptureBasin();
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

	delete viabProblem;
	cout<< " viab detruit retour main\n";
	postProcess();

}