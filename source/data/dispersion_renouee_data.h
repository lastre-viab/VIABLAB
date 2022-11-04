/*! \file  testZermelo.h
 *
 *
 *  \author: A. D�silles, LATSRE
 *  \brief  Fichier contant la d�finition d'un mod�le de viabilit� pour les tests
 *
 *
 */


#ifndef TESTDATA_H_
#define TESTDATA_H_


/***********************************************************
 * Some model specific parameters can be defined here
 *************************************************************/

int i = 0;  // entre 0 et N-1



int N = 7;

// propagules
double vect_alpha_R[7] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; 
double vect_E_A[7] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2}; 
double vect_E_R[7] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3}; 
double q_c_A = 1.0; 
double q_c_R = 1.0; 
double q_u_R = 1.0;


// zones
double vect_L[7] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; 
double vect_A[7] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; 
double c_i_star_max = 15.0 ;


// impact crue
double vect_l_R[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
double vect_l_A[7] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; 
double l_A = 0.1; 
double l_R = 0.0; 

// croissance plante
double h = -0.033 ;
double d = 0.097 ; 
double e = 5.0 ;  
double a_ini = 0.5 ; 



// reservoir crue
double r_max = c_i_star_max; 
double A = c_i_star_max/ 2.0 ;  



// dispersion

double mat_dispersion_A[7][7] = { {1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0},
    {0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0},
    {0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0},
    {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
    {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0} }  ;


double mat_dispersion_R[7][7] = { {1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0},
    {0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0},
    {0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0},
    {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
    {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0} }  ;


// bornes noyaux garanti calculés aux étapes précédentes (de taille i), donc pas de perturbation si i=0
/*
double bornes_noyaux_garantis_precedents_a[7] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};
double bornes_noyaux_garantis_precedents_b[7] = {2.0,2.0,2.0,2.0,2.0,2.0,2.0};
double bornes_noyaux_garantis_precedents_u[7] = {2.0,2.0,2.0,2.0,2.0,2.0,2.0};
 */
/*
double *bornes_noyaux_garantis_precedents_u ;
double *bornes_noyaux_garantis_precedents_u ;
double *bornes_noyaux_garantis_precedents_b ;
 */

// calculé dans la fonction loadModelData on met une taille plus grande que nécessaire à l'étape i
float tab [30];
double bornes_noyaux_garantis_precedents_a[7] ;
double bornes_noyaux_garantis_precedents_b[7] ;
double bornes_noyaux_garantis_precedents_u[7] ;






// durée (horison temporel)
double Tfinal = 10;  // un pas de temps = 1 an
double T_chgt_u_U = 2;  // on applique la borne U au controle les T_chgt_u_U premier pas de temps (années)

// contrainte dynamique (à faire varier)
double contrainte_a_i = 0.05; 

// contrainte controle (à faire varier)
double small_u_max = 5; 
double big_U_max = 10; 



////////////////////////////////////


/*! \var dim
 *  \brief State  dimension
 */
const int dim=3;


// controles 
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=1;

/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MIN[dimc]={0};

/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX[dimc]={big_U_max};

/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {11};

/*!
 * \brief Function  defining the mixed  constraints
 *
 * This function defines the set U(x) for admissible controls as function of the state
 * @param x state variable
 * @param u control variable
 * @return  value that caraterise the constraints set
 */
inline double constraintsXU( double * x, double * u )
{
  // doit renvoyer + \infty  si u \not \in U(x).
  double res = 1.0;
  if ( u[0] > small_u_max ) {res = PLUS_INF;
  } else {res = 1.0;};
  return res;
}




// tychastic (voir unused)


/*!
 * \var T maximum time horizon for the study
 */
double T= Tfinal;


/*! \var dicret_type
 *  \brief the  type of discretization scheme used  for the dynamics
 *  EE or 1 = Euler Explicit  scheme
 *  EI or 2 = Euler Implicit
 *  RK2I or 3 = RK2 Implicit (RK2 for -F)
 *  RK2E or 4 = RK2 Explicit (RK2 for F)
 */
const int discret_type=0;

/*!
 * \var dynType  defines the type of dynamics in the model
 *      1 or CC : continuous in time and space
 *      2 or DC : discrete time continuous space
 *      2 or DD : discrete time discrete space
 *      4 or HD : hybrid TODO
 */
//const int dynType=DC;
const int dynType=DD;

/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */
double STATE_MIN[dim]={0, 0, 0};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={r_max, 0.1, 1500 };  // l'aire max (voir vecteur de données)
// TODO mettre aire max des données

/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
//unsigned long long int nbPointsState[dim]    =  {401,401, 11};
unsigned long long int nbPointsState[dim]    =  {20,21,20};

/*
 * Direction le long de laquelle l'ensemble sera représenté par des suites de bits
 */
unsigned long long int dirTramage =1;

/*
 * Sélection de la méthode de représentation de l'ensemble
 * Ceparamètre détermine quelle classe sera utilisée pour les calculs
 *
 *    - BS = BitSet, représentation par fonction caractéristique , pour
 *                   tout type de calculs
 *    - MM = MicroMacro, représentation par valeurs réelles, pour les calculs
 *          d'ensembles épigraphiques , associés aux système micro-macro
 */
int gridMethod=BS;

/*
 * Sélection de l'ensemble à calculer
 */
int setType=VIAB;

/*
 * Nombre de raffinements successifs dans le calcul
 */
int refine= 0; //2;

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0};

/*!
 * \var dbFileName
 * the  name for the  main data base file for the roject
 */

string prefix="F_renouee_dispersion" ;
//string prefix= std::to_string(i) + "_F_renouee_dispersion-" ;


/*
 * Definition of the dynamics  and associated functions and constants
 */

/*
 * Paramètre qui indique si l'ensemble doit être calculé
 * dans le cas contrainte il devra pouvoir être chargé depuis un fichier pour une étude de
 * trajectoires
 */
int computeSet=1;

//int saveBoundary=1;
int saveBoundary=0;

const int nbTrajs=0;
double initPoints[dim*nbTrajs]={};

int sortieOKinf[dim]={0,0,0};
int sortieOKsup[dim]={1,0,1};
// int compute_tmin= 0;   // à quoi  compute_tmin correspond ? laissé dans unused



/////////////////////////////////////
//     fonctions utiles dynamique
////////////////////////////////////

// faire une fonction min / max ?


// définition fonction r
double dynamic_r(double r, double c) { return  min(r - c + A, r_max ) ; };


double p(double c, double c_i_star_max ) { return c/c_i_star_max ; };
double k_l_p( double c, double l, double c_i_star_max  ) { return l + (1 - l) * p(c, c_i_star_max) ; };
double g_j_A(double a, double u,  double c ) { return a * q_c_A * k_l_p(c,l_A,c_i_star_max) ; };


double lg_j(double a, double L, double Amax ) {
  double res = 1.0;
  if ( a <= L/Amax  ) { res = sqrt(a); } else {res = L;};
  return  res ; };

double g_j_R(double a, double u,  double c, double L, double Amax ) { return lg_j(a,L,Amax) * ( q_c_R * k_l_p(c,l_A,c_i_star_max) + q_u_R * u) ; };



double f_i(double a, double u, double a_new, double Amax ) { 
  double res = 1.0;
  if ( u <= 2.5  ) { res = min(a + a_new + e * (2.5 - u), Amax ) ; } else {res =  min( max( max(a + a_new + h*u + d, (a+a_new) / a_ini * (a + a_new + h * u + d) ) , 0.0), Amax ) ;};
  return  res ; };







void dynamics(double * x, double *u, double * image);
/*
void dynamics(double * x, double *u, double * image)
{

  // cas  i == 0 et i > 0 traités en même temps (boucle for)


  // perturbation
  double perturbation_a[7] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  double perturbation_b[7] = {2.0,2.0,2.0,2.0,2.0,2.0,2.0};
  double perturbation_u[7] = {2.0,2.0,2.0,2.0,2.0,2.0,2.0};
  double c = 0.0;

  // reservoir crue r
  // image[0] =
  double temp_R = dynamic_r(x[0], c)  ;

  // apport akène en i
   double temp_apport_akene = 0.0;
   int j;

   for (j = 0 ; j < i ; j++)
    {
        temp_apport_akene = temp_apport_akene + g_j_A (perturbation_a[j], perturbation_u[j] , c) * mat_dispersion_A[j][i];
    }    
   temp_apport_akene = temp_apport_akene + g_j_A (x[1], u[0] , c) * mat_dispersion_A[i][i]  ;


   // apport rhizome en i
   double temp_apport_rhizome = 0.0;

   for (j = 0 ; j < i ; j++)
    {
        temp_apport_rhizome = temp_apport_rhizome + ( vect_alpha_R[j] * perturbation_b[j] + g_j_R (perturbation_a[j], perturbation_u[j] , c, vect_L[j], vect_A[j]) ) * mat_dispersion_R[j][i];
    }    
   temp_apport_rhizome = temp_apport_rhizome + ( vect_alpha_R[i] * x[2] + g_j_R (x[1], u[0] , c, vect_L[i], vect_A[i]) ) * mat_dispersion_A[i][i]  ;


  // aire
  double a_new = (temp_apport_rhizome * vect_E_R[i] + temp_apport_akene * vect_E_A[i] ) * a_ini;

  //image[1] = 
  double temp_aire = f_i( x[1], u[0], a_new, vect_A[i] ) ;


  // accumulateur rhizome
  // image[2] =  
  double temp_prop = x[2] + temp_apport_rhizome * (1-vect_E_R[i]) ;


  //////// projections
  // R
  double pas_R = double( (STATE_MAX[0] - STATE_MIN[0]) / nbPointsState[0] )  ;
  int k_R = int( (temp_R - STATE_MIN[0]) / pas_R ) ;
  image[0] = STATE_MIN[0] + double(k_R * pas_R) ;

  // aire
  double pas_aire = double( (STATE_MAX[1] - STATE_MIN[1]) / nbPointsState[1] )  ;
  int k_aire = int( (temp_aire - STATE_MIN[1]) / pas_aire) ;
  image[1] = STATE_MIN[1] + double(k_aire * pas_aire) ;

  // prop
   double pas_prop = double( (STATE_MAX[2] - STATE_MIN[2]) / nbPointsState[2] )  ;
   int k_prop = int( (temp_prop - STATE_MIN[2]) / pas_prop ) ;
   image[2] = STATE_MIN[2] + double(k_prop * pas_prop) ;

}

 */




int dimc_ty=0;

double *CONTROL_MIN_ty;
/*! \var CONTROL_MAX_ty[dimc]
 *  \brief maximum values  for  tyches vector components
 *   \see controlParams
 */
double *CONTROL_MAX_ty;

unsigned long long *nbPointsControl_ty;




/// test perturbations
void dynamics_ty(double * x, double *u,  double *v, double * image);
void dynamics_ty(double * x, double *u,  double *v, double * image)
{

  // cas  i == 0 et i > 0 traités s&parément/ensemble? 

  // v de taille (au moins) 3*i+1
  // v = (u_0,a_0,b_0,...,u_{i-1},a_{i-1},b_{i-1},c)


  //  if (i > 0 ) {

  // perturbation
  double perturbation_u[7]; //= {2.0,2.0,2.0,2.0,2.0,2.0,2.0};  (7: taille max, au pire on ne remplit pas)
  double perturbation_a[7]; //= {1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  double perturbation_b[7]; //= {2.0,2.0,2.0,2.0,2.0,2.0,2.0};

  double c; // = 0.0;
  c = v[3*i+1];


  // initialize elements of arrays with v
  for ( int j = 0; j <= (i-1); j++ ) {
    perturbation_u[ j ] = v[ 3*j ];
    perturbation_a[ j ] = v[ 3*j +1 ];
    perturbation_b[ j ] = v[ 3*j +2 ];
  }



  // reservoir crue r
  // image[0] =
  double temp_R = dynamic_r(x[0], c)  ;

  // apport akène en i
  double temp_apport_akene = 0.0;
  //int j;

  for ( int j = 0 ; j < i ; j++)
    {
    temp_apport_akene = temp_apport_akene + g_j_A (perturbation_a[j], perturbation_u[j] , c) * mat_dispersion_A[j][i];
    }    
  temp_apport_akene = temp_apport_akene + g_j_A (x[1], u[0] , c) * mat_dispersion_A[i][i]  ;


  // apport rhizome en i
  double temp_apport_rhizome = 0.0;

  for ( int j = 0 ; j < i ; j++)
    {
    temp_apport_rhizome = temp_apport_rhizome + ( vect_alpha_R[j] * perturbation_b[j] + g_j_R (perturbation_a[j], perturbation_u[j] , c, vect_L[j], vect_A[j]) ) * mat_dispersion_R[j][i];
    }    
  temp_apport_rhizome = temp_apport_rhizome + ( vect_alpha_R[i] * x[2] + g_j_R (x[1], u[0] , c, vect_L[i], vect_A[i]) ) * mat_dispersion_A[i][i]  ;


  // aire
  double a_new = (temp_apport_rhizome * vect_E_R[i] + temp_apport_akene * vect_E_A[i] ) * a_ini;

  //image[1] = 

  double temp_aire = f_i( x[1], u[0], a_new, vect_A[i] ) ;





  // accumulateur rhizome

  // image[2] =  

  double temp_prop = x[2] + temp_apport_rhizome * (1-vect_E_R[i]) ;




  //////// projections

  // R

  double pas_R = double( (STATE_MAX[0] - STATE_MIN[0]) / nbPointsState[0] )  ;

  int k_R = int( (temp_R - STATE_MIN[0]) / pas_R ) ;

  image[0] = STATE_MIN[0] + double(k_R * pas_R) ;



  // aire

  double pas_aire = double( (STATE_MAX[1] - STATE_MIN[1]) / nbPointsState[1] )  ;

  int k_aire = int( (temp_aire - STATE_MIN[1]) / pas_aire) ;

  image[1] = STATE_MIN[1] + double(k_aire * pas_aire) ;



  // prop

  double pas_prop = double( (STATE_MAX[2] - STATE_MIN[2]) / nbPointsState[2] )  ;

  int k_prop = int( (temp_prop - STATE_MIN[2]) / pas_prop ) ;

  image[2] = STATE_MIN[2] + double(k_prop * pas_prop) ;



  /*
/////////////
/////////////   si i=0 pas de perturbation
   }  else {    

  double c; // = 0.0;
  c = v[3*i+1];

  // reservoir crue r
  // image[0] =
  double temp_R = dynamic_r(x[0], c)  ;

  // apport akène en i
   double temp_apport_akene = 0.0;
   temp_apport_akene = temp_apport_akene + g_j_A (x[1], u[0] , c) * mat_dispersion_A[i][i]  ;


   // apport rhizome en i
   double temp_apport_rhizome = 0.0;    
   temp_apport_rhizome = temp_apport_rhizome + ( vect_alpha_R[i] * x[2] + g_j_R (x[1], u[0] , c, vect_L[i], vect_A[i]) ) * mat_dispersion_A[i][i]  ;


  // aire
  double a_new = (temp_apport_rhizome * vect_E_R[i] + temp_apport_akene * vect_E_A[i] ) * a_ini;

  //image[1] = 
  double temp_aire = f_i( x[1], u[0], a_new, vect_A[i] ) ;


  // accumulateur rhizome
  // image[2] =  
  double temp_prop = x[2] + temp_apport_rhizome * (1-vect_E_R[i]) ;


  //////// projections
  // R
  double pas_R = double( (STATE_MAX[0] - STATE_MIN[0]) / nbPointsState[0] )  ;
  int k_R = int( (temp_R - STATE_MIN[0]) / pas_R ) ;
  image[0] = STATE_MIN[0] + double(k_R * pas_R) ;

  // aire
  double pas_aire = double( (STATE_MAX[1] - STATE_MIN[1]) / nbPointsState[1] )  ;
  int k_aire = int( (temp_aire - STATE_MIN[1]) / pas_aire) ;
  image[1] = STATE_MIN[1] + double(k_aire * pas_aire) ;

  // prop
   double pas_prop = double( (STATE_MAX[2] - STATE_MIN[2]) / nbPointsState[2] )  ;
   int k_prop = int( (temp_prop - STATE_MIN[2]) / pas_prop ) ;
   image[2] = STATE_MIN[2] + double(k_prop * pas_prop) ;

// fin else i == 0
   };
   */

}








//////////////////////////////////////
//// test ensemble perturbation

inline double constraintsXUV( double * x, double * u, double * v )
{



  double res = 1.0;

  // perturbation (comme dans la dynamique)
  if(i == 0 ) { 
    double c; // = 0.0;
    c = v[3*i+1];
    if ( c <= r_max ) { res = 1; } else { res = PLUS_INF; }  // normalement la condition est vérifiée avec le choix de l'ensemble de contrainte
  }    else {   // i>0

    double perturbation_u[i-1];
    double perturbation_a[i-1];
    double perturbation_b[i-1];

    double c; // = 0.0;
    c = v[3*i+1];

    double res_c;
    if ( c <= r_max ) { res_c = 1; } else {res_c = PLUS_INF;};  // normalement la condition est vérifiée avec le choix de l'ensemble de contrainte


    int j;

    // initialize elements of arrays with v
    for ( j = 0; j <= (i-1); j++ ) {
      perturbation_u[ j ] = v[ 3*j ]; 
      perturbation_a[ j ] = v[ 3*j +1 ]; 
      perturbation_b[ j ] = v[ 3*j +2 ]; 
    }



    ///// contrainte sur u
    double res_u = 1.0;
    int j_u = 0;
    while ( ( j_u <= (i-1) ) & (perturbation_u[j_u] <= bornes_noyaux_garantis_precedents_u[j_u] )  )
      {
      j_u ++;
      }
    if ( j_u == i ) { res_u = 1; } else {res_u = PLUS_INF;};

    ///// contrainte sur a
    double res_a = 1.0;
    int j_a = 0;
    while ( ( j_a <= (i-1) ) & (perturbation_a[j_a] <= bornes_noyaux_garantis_precedents_a[j_a] )  )
      {
      j_a ++;
      }
    if ( j_a == i ) { res_a = 1; } else {res_a = PLUS_INF;};


    ///// contrainte sur b
    double res_b = 1.0;
    int j_b = 0;
    while ( ( j_b <= (i-1) ) & (perturbation_b[j_b] <= bornes_noyaux_garantis_precedents_b[j_b] )  )
      {
      j_b ++;
      }
    if ( j_b == i ) { res_b = 1; } else {res_b = PLUS_INF;};


    // bilan
    if ( res_u == PLUS_INF ) { res = PLUS_INF; 
    } else if ( res_a == PLUS_INF ) { res = PLUS_INF;
    } else if ( res_b == PLUS_INF ) { res = PLUS_INF;
    } else if ( res_c == PLUS_INF ) { res = PLUS_INF; } else { res = 1; };



    // fin de else (i>0)
  }

  return res;

}














/*      *****************************************
 *  Definition of constraints and target
 *************************************************** */

/*!
 * \brief Function  defining the state   constraints, corresonds  to k(x)
 *
 * This function defines the set K for admissible  states
 * @param x state variable
 * @return  value that caraterise the constraints set
 */


inline double constraintsX( double * x )
{
  double r = x[0];
  double a = x[1];
  double b = x[2];

  double res = 1.0;
  if(a <= contrainte_a_i ) {res = 1.0;} else {res = PLUS_INF;};
  return res;
  //return 1.0;
}




////////////
void static updateParams(int local_i, double local_big_U_max, double local_contrainte_a_i){
  i = local_i ;
  //prefix= std::to_string(i) + "_F_renouee_dispersion-" ;
  big_U_max = local_big_U_max ;
  CONTROL_MAX[dimc]={big_U_max};
  contrainte_a_i = local_contrainte_a_i;
}




////////////
void loadModelData();
void loadModelData(){


  /*
   * AD: Initialisation des paramètres du contrôle tychastique
   *
   * en fonction de i
   */




  // open a file in read mode.
  //ifstream infile;
  //infile.open("test_perturbations.csv");

  /*
   double *myArray
   ? value;
   int j = 0;
   while(inFileStr >> value)
   { 
      myArray[j] = value;
      j++;
   }
   */


  /*
    std::ifstream  data("test_perturbations.csv");
    std::string line;
    std::vector<std::vector<std::string> > parsedCsv;
    while(std::getline(data,line))
    {
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<std::string> parsedRow;
        while(std::getline(lineStream,cell,','))
        {
            parsedRow.push_back(cell);
        }

        parsedCsv.push_back(parsedRow);
    }
   */

  // close the opened file.
  //infile.close();


  /*
   float data[2][2];
    std::ifstream file("test_perturbations.csv");

    for(int row = 0; row < 2; ++row)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;

        std::stringstream iss(line);

        for (int col = 0; col < 2; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;

            //std::stringstream convertor(val);
            int convertor = std::stof (val);
            convertor >> data[row][col];
        }
    }
   */

  /*
    std::string line, val;                  // string for line & value //
    std::vector<std::vector<int>> array;    // vector of vector<int>  //

    while (std::getline (f, line)) {        // read each line //
        std::vector<int> v;                 // row vector v //
        std::stringstream s (line);         // stringstream line //
        while (getline (s, val, ','))       // get each value (',' delimited) //
            v.push_back (std::stoi (val));  // add to row vector //
        array.push_back (v);                // add row vector to array //
    }
   */


  ///////////////////////////////////
  ///////////////////////////////////

  //float tab[30];
  //tab=new double [30];
  int TAILLE_MAX = 1000 ;
  FILE* fichier = NULL;
  char chaine[TAILLE_MAX] = "";
  int j;
  fichier = fopen("test_perturbation.csv", "r");
  j = 0;
  if (fichier != NULL)
    {
    while (fgets(chaine, TAILLE_MAX, fichier) != NULL) // On lit le fichier tant qu'on ne reçoit pas d'erreur (NULL)
      {
      printf("%s", chaine); // On affiche la chaîne qu'on vient de lire

      tab[j] = atof(chaine);
      printf("%f\n", tab[j]);

      j =j+1;
      }

    fclose(fichier);
    };



  // initialize elements of arrays with tab (si i=0, pas de valeurs dans le tableau)
  int jj = 0;
  for ( jj = 0; jj <= (i-1); jj++ ) {
    bornes_noyaux_garantis_precedents_u[ jj ] = tab[ 3*jj ];
    bornes_noyaux_garantis_precedents_u[ jj ] = tab[ 3*jj +1 ];
    bornes_noyaux_garantis_precedents_b[ jj ] = tab[ 3*jj +2 ];
  }

  dimc_ty=3*i;

  CONTROL_MIN_ty= new double [dimc_ty];

  for(int k=0;k<dimc_ty; k++)
    {
    CONTROL_MIN_ty[k]=0.0;// remplir dans cette boucle les bornes inf des tyches
    }

  CONTROL_MAX_ty= new double [dimc_ty];

  for(int k=0;k<dimc_ty; k++)
    {
    CONTROL_MAX_ty[k]=1.0;// remplir dans cette boucle les bornes sup des tyches
    }

  nbPointsControl_ty = new unsigned long long [dimc_ty];
  for(int k=0;k<dimc_ty; k++)
    {
    nbPointsControl_ty[k]=2.0;// remplir dans cette boucle les nombres de valeurs possibles des tyches
    }



  ///// test
  //printf("yolo  \n");
  //printf("test_perturb[0] ",data[0]);
  //printf("test_perturb[0] %d",data[0]);  
  //printf("sort[0] %d",sortieOKsup[0]);
  //string test_affiche= std::to_string(i) ;
  //printf(test_affiche);
  //printf("\n");
  //printf("j %d",j);

};




#include "dispersion_renouee_unused.h"
#endif /* TESTDATA_H_ */
