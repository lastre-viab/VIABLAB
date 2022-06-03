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

//#include "../data/WeakDeclarations.h"
#include "../data/DefaultValues.h"
//#include  "../data/equilibres4D_data.h"  //-- a jour, non reg OK
//#include  "../data/ex1Viabi2D_data.h"   //-- a jour, non reg OK
//#include  "../data/testPendule_data.h"   //-- a jour, non reg OK
//#include  "../data/resilience_data.h"   //-- a jour, non reg OK
#include "../data/PSP_dataBis.h"
//#include  "../data/julia2D_data.h"  //-- a jour, non reg OK
//#include  "../data/ex1Viabi2D_data.h"  //-- a jour, non reg OK

//#include  "../data/zermelo_tmin_data.h"  //-- a jour, non reg OK
//#include  "../data/zermelo_Lmin_data.h"  //-- a jour, non reg OK
//#include  "../data/testZermelo_bitSet.h"  //-- a jour, non reg OK
//#include  "../data/maupertuis_data.h"

//#include  "../data/ex4_Viabi2D_data.h"
//#include "../data/ex2_EcoPolut_data.h"

/*
 * Code de la fonction d'initialisation de paramètres
 */

//#include "../data/dataAgroEcoDivMultiParcels.h"
//#include "../data/dataAgroEcoDivMultiParcelsD4.h"

//#include "../data/dataAgroEcoDivTest.h"
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
		loadModelData();
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
			viabProblem->ViabilityKernel( true, 50);
		}
		if(avp.SET_TYPE==CAPT) { // option de calcul de noyau de viabilité a été choisie
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
if(avp.SAVE_PROJECTION | avp.SAVE_SLICE_BOUND)
{
	viabProblem->saveViableSets();
}
	delete viabProblem;
	cout<< " viab detruit retour main\n";
	postProcess();

}
