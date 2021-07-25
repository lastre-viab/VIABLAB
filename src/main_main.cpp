/*
 * main.cpp
 *
 *  Created on: 10 sept. 2013
 *      Author: ANYA
 */
/*
 * Inclusions des classes utilisées pour le calcul
 */
#include "../include/ViabiBitSet.h"
#include "../include/Viabi.h"
#include "../include/ViabiMicroMacro.h"
#include "../include/GridMicroMacro.h"
#include  "../include/SysDyn.h"
#include "../include/GridBitSet.h"
#include "../include/defs.h"

/*
 * inclusion des déclarations du modèle
 */
//#include  "../data/equilibres4D_data.h"
//#include  "../data/ex1Viabi2D_data.h"

//#include  "../data/julia2D_data.h"
//#include  "../data/ex4_Viabi2D_data.h"
//#include  "../data/ex3Dim_Viabi2D_data.h"
//#include  "../data/ex1_Viabi_multiDim_data.h"
//#include  "../data/ex3bis_Viabi2D_data.h"
//#include "../data/ex2_EcoPolut_data.h"

//#include "../data/userModel.h"
//#include "../data/test_capt_discret_data.h"
//#include "../data/corentin_no_control_feed_synth_data_funct_traj_v3.h"

#include "../data/corentin_no_control_feed_synth_data_funct_traj_v3.h"
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
  SysDyn * sd;

  Grid_BitSet *grBS;
  ViabiBitSet *viabBS;

  GridMicroMacro *grMM;
  ViabiMicroMacro *viabMM;
  // déclaration de structures pour les paramètres
  gridParams gp;
  algoViabiParams avp;
  controlParams cp;
  systemParams sp;
  //initialisation des structures à partir des données de l'utilisateur
  initParams(gp, avp, cp, sp);

  if(gridMethod==BS)
    {
    /*
     * Calculs réalisés avec une représentation
     * booléenne
     */
    /*
     * instanciation de la grille
     */
    grBS=new Grid_BitSet( gp);
    grBS->printGrid();
    /*
     * instanciation du système dynamique
     */
    sp.L_LIP=0.0;
    sp.L_MF=0;


    sd= new SysDyn(sp, dim, cp, grBS);
    // instanciation de la classe viabilité
    viabBS=new ViabiBitSet(grBS, sd,avp);

    if(computeSet){
      // l'ensemble doit être  recalculé
      gettimeofday(&tim_glob,NULL);              //mesure le temps d'execution
      t1_glob=tim_glob.tv_sec+(tim_glob.tv_usec/1000000.0);

      if(setType==VIAB) { // option de calcul de noyau de viabilité a été choisie
        // Initialisation de l'ensemble de contraintes

        viabBS->setK0();
      }
      //  raffinements successifs ( y compris le premier calcul si refine=0)
      int refIter=-1;
      int seuilArret=1;
      while( refIter<refine){
        gettimeofday(&tim,NULL);              //mesure le temps d'execution
        t1=tim.tv_sec+(tim.tv_usec/1000000.0);
        //  Calcul du noyau de viabilité
        if(setType==VIAB) {
          viabBS->ViabKer(1,seuilArret);
        }
      //  cout<< "=============================================================================="<<endl;
        gettimeofday(&tim,NULL);
        t2=tim.tv_sec+(tim.tv_usec/1000000.0);
        elapsed_time=(double)((t2-t1));
        cout << "Elapsed time : " << elapsed_time << " sec." << endl << endl;
     //   cout<< "=============================================================================="<<endl;
        /*  * sauvegardes
         */
        os<<"OUTPUT/"<<prefix<<"-viab-"<<refIter+1<<".dat";
        fileName=os.str();
        os.str("");
        grBS->saveValOnGrid(fileName);

        if(saveCoupe){
          os<<"OUTPUT/"<<prefix<<"-viab-"<<refIter+1<<"-Slice"<<".dat";
          fileName=os.str();
          os.str("");
          grBS->saveCoupe(fileName);
        }
        if(saveCoupeBound){
          os<<"OUTPUT/"<<prefix<<"-viab-"<<refIter+1<<"-SliceBound"<<".dat";
          fileName=os.str();
          os.str("");
          grBS->saveCoupeBoundary(fileName);
        }
        if(saveBoundary){
          os<<"OUTPUT/"<<prefix<<"-viab-"<<refIter+1<<"-bound.dat";
          fileName=os.str();
          os.str("");
          grBS->saveBoundary(fileName);
        }

        if(saveProjection)
          {
          os<<"OUTPUT/"<<prefix<<"-viab-"<<refIter+1<<"-proj"<<".dat";
          fileName=os.str();
          os.str("");
          /*
           *  calcul et sauvegarde  de la projection du  noyau
           */
          grBS->saveProjetion(fileName, projection);
          }

        refIter++;
        seuilArret*=2;
        /*
         * raffinement du maillage
         */

        if(refIter<refine)
          grBS->refine();
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
      os<<"OUTPUT/"<<prefix<<"-viab-"<<refine<<".dat";
      fileName=os.str();
      os.str("");
      grBS->loadSet(fileName);
      if(saveCoupe){
        os<<"OUTPUT/"<<prefix<<"-viab-"<<refine<<"-Slice"<<".dat";
        fileName=os.str();
        os.str("");
        grBS->saveCoupe(fileName);
      }
      if(saveCoupeBound){
        os<<"OUTPUT/"<<prefix<<"-viab-"<<refine<<"-SliceBound"<<".dat";
        fileName=os.str();
        os.str("");
        grBS->saveCoupeBoundary(fileName);
      }
      if(saveProjection)
        {
        os<<"OUTPUT/"<<prefix<<"-viab-"<<refine<<"-proj"<<".dat";
        fileName=os.str();
        os.str("");
        /*
         *  calcul et sauvegarde  de la projection du  noyau
         */
        grBS->saveProjetion(fileName, projection);
        }
      }

    /*
     * Si les trajectoires doivent être calculées
     */
    if(nbTrajs>0)
      {
      for(int tr=0;tr<nbTrajs;tr++)
        {
        if(typeTraj==VD)
          {
          os<<"OUTPUT/"<<prefix<<"-traj-"<<tr+1<<".dat";
          fileName=os.str();
          os.str("");

          viabBS->computeViableTrajectorySetVal(initPoints+tr*dim, T, fileName);
          }
        if(typeTraj==VL)
          {
          os<<"OUTPUT/"<<prefix<<"-traj-H-"<<tr+1<<".dat";
          fileName=os.str();
          os.str("");

          viabBS->computeViableTrajectoryHeavy(initPoints+tr*dim, initControls+tr*dimc,T, fileName);
          }
        }
      }

    delete viabBS;

    delete grBS;

    delete sd;
    cout<< "Fin du programme VIABLAB\n";
    postProcess();
    }

  if(gridMethod==MM)
    {
    /*
     * Calculs réalisés avec une représentation
     * booléenne
     */
    /*
     * instanciation de la grille
     */
    grMM=new GridMicroMacro( gp);
    grMM->printGrid();
    /*
     * instanciation du système dynamique
     */
    if(compute_tmin)
      {
      sp.L_LIP=1.0;
      sp.L_MF=T;
      }
    else{
      sp.L_LIP=l_Lip;
      sp.L_MF=l_max;
    }
    sd= new SysDyn(sp, dim, cp, grMM);
    // instanciation de la classe viabilité
    viabMM=new ViabiMicroMacro(grMM, sd,avp);

    if(computeSet){
      // l'ensemble doit être  recalculé
      gettimeofday(&tim_glob,NULL);              //mesure le temps d'execution
      t1_glob=tim_glob.tv_sec+(tim_glob.tv_usec/1000000.0);

      if(setType==VIAB) { // option de calcul de noyau de viabilité a été choisie
        // Initialisation de l'ensemble de contraintes
        cout<< " initialise constr \n";
        viabMM->initializeConstrains();
        cout<< "OK. viab \n";
        viabMM->viabKerValFunc();
      }
      if(setType==CAPT) { // option de calcul de noyau de viabilité a été choisie
        // Initialisation de l'ensemble de contraintes
        cout<< " initialise target  \n";
        viabMM->initialiseTargetHJB();
        viabMM->captBasinEpi();
      }

      cout<< "=============================================================================="<<endl;
      gettimeofday(&tim_glob,NULL);
      t2_glob=tim_glob.tv_sec+(tim_glob.tv_usec/1000000.0);
      elapsed_time_glob=(double)((t2_glob-t1_glob));
      cout << "Elapsed time  GLOBAL: " << elapsed_time_glob << " sec." << endl << endl;
      cout<< "=============================================================================="<<endl;

      if(saveSubLevel){
        os<<"OUTPUT/"<<prefix<<"-subLevel.dat";
        fileName=os.str();
        os.str("");
        grMM->saveSubLevelset(level, fileName);
      }
      if(saveBoundary){
        os<<"OUTPUT/"<<prefix<<"-captBound.dat";
        fileName=os.str();
        os.str("");
        grMM->saveValOnGrid(fileName);
      }
      if(saveProjection)
        {
        os<<"OUTPUT/"<<prefix<<"-viab-"<<refine<<"-proj"<<".dat";
        fileName=os.str();
        os.str("");
        /*
         *  calcul et sauvegarde  de la projection du  noyau
         */
        grMM->saveProjetion(fileName, projection);
        }



    }//fin de if(computeSet)
    else
      { // on charge dans la mémoire l'ensemble calculé et enregistré
      // correspondant au dernier raffinement
      /*  os<<"OUTPUT/"<<prefix<<"-viab-"<<refine<<".dat";
      fileName=os.str();
      os.str("");
      grMM->loadSet(fileName);
      if(saveCoupe){
        os<<"OUTPUT/"<<prefix<<"-viab-"<<refine<<"-Slice"<<".dat";
        fileName=os.str();
        os.str("");
        grMM->saveCoupe(fileName);
      }*/

      if(saveSubLevel){
        os<<"OUTPUT/"<<prefix<<"-subLevel.dat";
        fileName=os.str();
        os.str("");
        grMM->saveSubLevelset(level, fileName);
      }
      if(saveBoundary){
        os<<"OUTPUT/"<<prefix<<"-captBound.dat";
        fileName=os.str();
        os.str("");
        grMM->saveValOnGrid(fileName);
      }
      if(saveProjection)
        {
        os<<"OUTPUT/"<<prefix<<"-viab-"<<refine<<"-proj"<<".dat";
        fileName=os.str();
        os.str("");
        /*
         *  calcul et sauvegarde  de la projection du  noyau
         */
        grMM->saveProjetion(fileName, projection);
        }
      }

    /*
     * Si les trajectoires doivent être calculées
     */
    if(nbTrajs>0)
      {
      bool success[nbTrajs];
      for(int tr=0;tr<nbTrajs;tr++)
        {
        os<<"OUTPUT/"<<prefix<<"-traj-"<<tr+1<<".dat";
        fileName=os.str();
        os.str("");

        viabMM->computeOptimalCaptTrajectory(initPoints+tr*dim, fileName, success[tr]);
        }
      }

    delete viabMM;
    cout<< " viab detruit retour main\n";
    delete grMM;
    cout<< "grid detruit retour main\n";
    delete sd;
    cout<< " objets detruits\n";
    postProcess();
    }

}
