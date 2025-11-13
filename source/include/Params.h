#ifndef PARAMS_H
#define PARAMS_H

#include <string>
using std::string;

#include "Enums.h"
#include "ControlPickStrategyName.h"
#include "TrajectoryHelpers.h"

/*!
 *  \struct gridParams
 *  \brief Grid parameters container
 */
struct gridParams
    {
    unsigned long long int DIM;       //!< dimension de la  grille
    unsigned long long int *NBPOINTS; //!<   tableau indicant le nombre de points par axe
    double *LIMINF;    //!<   tableau  de limites inf pour chaque variable
    double *LIMSUP;    //!<   tableau  de limites sup pour chaque variable
    bool *PERIODIC; //!<   tableau  de boolens pour indiquer   les variables  qui sont p�riodiques (1=p�riodique, 0=non)
    int GRID_TYPE; //!<  0=les points de calcul sont les noeuds de la grille, 1= les points de calcul sont les milieux des mailles
    string FILE_PREFIX;
    int GRID_MAIN_DIR; //!< for bitSet grids only :  la grille sera stock�e comme tableau de segments paralleles  a cette dimension
    bool *SLICE_DIRS;
    double *SLICE_VALUES;
    unsigned long long int *SLICE_VALUES_FD;
    bool *SORTIE_OK_INF;
    bool *SORTIE_OK_SUP;
    int OMP_THREADS;
    GridMethod GRID_METHOD;
    bool IS_HYBRID;
    unsigned long long int DIM_HD;
    unsigned long long int DIM_HC;
    unsigned long long int *NBPOINTS_HD; //!<   tableau indicant le nombre de points par axe

    };

/*!
 *  \struct controlParams
 *  \brief controle parameters container
 *  \see SysDyn
 */

struct controlParams
    {
    unsigned long long int DIMC;        //!<  dimension  de controles
    unsigned long long int *NBPOINTSC; //!<  nombre de points par axe en cas de discr�tisation
    double *LIMINFC;      //!<  limites inf de chaque variable de controle
    double *LIMSUPC;       //!<  limites sup de chaque variable de controle

    unsigned long long int DIM_TY;        //!<  dimension  de controles
    unsigned long long int *NBPOINTS_TY; //!<  nombre de points par axe en cas de discr�tisation
    double *LIMINF_TY;      //!<  limites inf de chaque variable de controle
    double *LIMSUP_TY;       //!<  limites sup de chaque variable de controle

    unsigned long long int DIM_HT;        //!<  dimension  de controles
    unsigned long long int *NBPOINTS_HT; //!<  nombre de points par axe en cas de discr�tisation
    double *LIMINF_HT;      //!<  limites inf de chaque variable de controle
    double *LIMSUP_HT;       //!<  limites sup de chaque variable de controle
    };

/*!
 * \struct systemParams
 * \brief  dynamical system functions cntainer
 * \see SysDyn
 */

struct systemParams
    {
    void (*DYNAMICS)(const double*, const double*, double*);
    void (*DYNAMICS_TYCH)(const double*, const double*, const double*, double*);
    void (*DYNAMICS_FD)(const unsigned long long int*, const unsigned long long int*, unsigned long long int*);
    void (*DYNAMICS_TYCH_FD)(const unsigned long long int*, const unsigned long long int*, const unsigned long long int*, unsigned long long int*);

    void (*DYNAMICS_HYBRID_D)(const double*, const unsigned long long int*, const unsigned long long int*, unsigned long long int*);
    void (*DYNAMICS_HYBRID_C)(const double*, const unsigned long long int*, const double*, double*);
    double (*CONSTR_XU_HYBRID)(const double*, const unsigned long long int*, const double*, const unsigned long long int*);
    void (*RESET_MAP_HYBRID)(const double *, const unsigned long long int*, const unsigned long long int*, double *, const unsigned long long int*);

    double (*CONSTR_XU)(const double*, const double *u);
    double (*CONSTR_XV_TYCH)(const double*, const double *v);
    double (*CONSTR_XU_fd)(const unsigned long long int*, const unsigned long long int*);
    double (*CONTROL_ELIGIBILITY_FOR_TRAJ_fd)(const unsigned long long int*, const unsigned long long int*, const unsigned long long int*);

    double (*CONSTR_X)(const double*);
    double (*CONSTR_X_HYBRID)(const double*, const unsigned long long int*);
    double (*CONSTR_X_fd)(const unsigned long long int*);
    double (*DYN_CONSTR_FOR_TRAJ)(const double*, double*);
    double (*TARGET)(const double*);
    double (*TARGET_FD)(const unsigned long long int*);
    double (*L_FUNC_FD)(const unsigned long long int*, const unsigned long long int*);
    double (*L_FUNC_TYCH_FD)(const unsigned long long int*, const unsigned long long int*, const unsigned long long int*);

    double (*L_FUNC)(const double*, const double*);
    double (*L_FUNC_TYCH)(const double*, const double*, const double*);
    double (*MU_FUNC_FD)(const unsigned long long int*, const unsigned long long int*);
    double (*M_FUNC)(const double*, const double*);
    double (*M_FUNC_TYCH)(const double*, const double*, const double*);
    void (*JACOBIAN)(const double *x, const double *u, double **jacob);
    void (*JACOBIAN_HYBRID)(const double *x, const unsigned long long int*, const double *u, double **jacob);
    void (*JACOBIAN_TYCH)(const double *x, const double *u, const double *v, double **jacob);
    void (*LOCAL_DYN_BOUNDS)(const double *x, double *res);
    void (*LOCAL_DYN_BOUNDS_HYBRID)(const double *x, const unsigned long long int*, double *res);

    ComputeMethod COMPUTE_LC;
    double LIP;
    double L_LIP;
    bool globDeltat;
    ComputeMethod COMPUTE_MF;
    double MF;
    double L_MF;    
    TimeDiscretizationScheme SCHEME;
    bool SCALE_PARAM;
    DynType DYN_TYPE;     // type of dynamics : FD, CD,
    FdDynType FD_DYN_TYPE; // pour indiquer le type de repr�sentation de la dynamique : fonction ou fihciers r�tro
    std::string RETRO_FILE_NAME;
    };

struct algoViabiParams
    {

    SetType SET_TYPE;
    bool COMPUTE_SET;
    int NB_OMP_THREADS;
    bool COMPUTE_TMIN;

    string FILE_PREFIX;
    int GRID_REFINMENTS_NUMBER;
    bool INTERMEDIATE_SAVINGS;
    bool SAVE_SLICE;
    bool SAVE_SLICE_BOUND;
    bool SAVE_BOUNDARY;
    bool SAVE_PROJECTION;
    bool SAVE_SUBLEVEL;
    bool SAVE_VIAB_LIGHT;
    double LEVEL;
    unsigned long long int *PROJECTION;
    int ITERATION_STOP_LEVEL;
    TargetOrDeparture TARGET_OR_DEPARTURE;
    };

struct algoViabiBoolParams
    {
    unsigned long long int NB_GO_RAM;
    };

struct tycheParams {
    // Nombre maximal de tirage de ce tyché en cas de non correspondance de la contrainte état-tyché
    int MAX_NB_REROLLS;

    // Fonctions utilisateur
    userTyche_t USER_TYCHE;
    probabilityDensity_t PROBABILITY_DENSITY;
    cumulativeDistribution_t CUMULATIVE_DISTRIBUTION;

    // Valeur du tyché pour une distribution constante
    double CONSTANT_TYCHE_VALUE;
    TycheDistribution TYCHE_DISTRIBUTION;
};

struct trajectoryParams {

    double *INIT_POINT;
    double INIT_VALUE;
    double INIT_VALUE_FD;
    unsigned long long int *INIT_POINT_FD;
    double *INIT_CONTROL;

    double maxTime;    
    // Rayon de la bulle sur chaque dimenseion
    double *BUBBLE_RADIUS;    

    // Noms des stratégies des trajectoires par stratégies
    ControlPickStrategyName *STRATEGIES;
    // Tableau d'entier utilisé pour la seed de TychePicker
    unsigned long long int *SEED;
    // Tableau de taille DIM_TY des paramètres de trajectoires tychastiques
    tycheParams *TYCHE_PARAMS;

    // Fonctions utilisateur
    controlWeight_t CONTROL_WEIGHT;
    indexSorter_t SORT_INDEXES;
    neighborValidator_t IS_VALID_NEIGHBOR;
    temporalControl_t TEMPORAL_CONTROL;

    // Angle maximal en radians d'une stratégie SMOOTH
    double MAX_ANGLE_RADIANS;

    int SEED_LENGTH;
    int NB_STRATEGIES;
    // Nombre de pas de temps de la trajectoire réelle calculés pour chaque choix de contrôle
    int REAL_TIME_STEPS_PER_DISCRETE_STEP;

    BubbleInterpretation BUBBLE_INTERPRETATION;
    TypeTraj TRAJECTORY_TYPE;

    // Doit-on enregistrer dans le fichier de sortie les trajectoires ayant contribué au choix de contrôle ?
    bool SAVE_PICKING_STRATEGY;
    // Les stratégies doivent-elle uniquement choisir des contrôles viables garantis (viables pour tout tyché en position courante) ou alors seulement viable pour le tyché courant ?
    bool ARE_STRATEGIES_GUARANTEED;
};

#endif /* PARAMS_H */
