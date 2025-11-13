#include <exception>

#include "../include/WeakDeclarations.h"

/*
 Fichier contenant les implémentations par défaut des fonctions redéfinissable par l'utilisateur dans son modèle
 */

// Si le fichier est vide (valeur par défaut) le code doit renvoyer une erreur

class ParametersManager;

void loadModelData(const ParametersManager*)
    {

    }

double l(const double *x, const double *u)
    {
    return 1.0;
    }

double m(const double *x, const double *u)
    {
    return 0.0;
    }

double l_tych(const double *x, const double *u, const double *v)
    {
    return 1.0;
    }

double l_tych_fd(const unsigned long long int *x, const unsigned long long int *u, const unsigned long long int *v)
    {
    return 1.0;
    }

double m_tych(const double *x, const double *u, const double *v)
    {
    return 0.0;
    }

void dynamics_hybrid(const double *x, const double *u, double *image)
    {

    }

void postProcess(const ParametersManager *PM)
    {

    }

double target(const double *x)
    {
    return PLUS_INF;
    }

double constraintsX(const double *x)
    {
    return 0.0;
    }

double constraintsXU(const double *x, const double *u)
    {
    return 0.0;
    }

void localDynBounds(const double *x, double *res)
    {

    }

double constraintsXU_fd(const unsigned long long int *x, const unsigned long long int *u)
    {
    return 0.0;
    }

void jacobian(const double *x, const double *u, double **jacob)
    {

    }

void jacobian_tych(const double *x, const double *u, const double *v, double **jacob)
    {

    }

void dynamics(const double *x, const double *u, double *image)
    {

    }

void dynamics_tych(const double *x, const double *u, const double *v, double *image)
    {

    }

void dynamics_fd(const unsigned long long int*, const unsigned long long int*, unsigned long long int*)
    {

    }

void dynamics_tych_fd(const unsigned long long int*, const unsigned long long int*, const unsigned long long int*, unsigned long long int*)
    {

    }

double dynConstraintsForTraj(const double *x, double *image)
    {
    return 1.0;
    }

double constraintsX_fd(const unsigned long long int*)
    {
    return 0.0;
    }

double target_fd(const unsigned long long int *x)
    {
    return 0.0;
    }

double l_fd(const unsigned long long int *x, const unsigned long long int *u)
    {
    return 0.0;
    }

double l_fd_tych(const unsigned long long int *x, const unsigned long long int *u, const unsigned long long int *v)
    {
    return 0.0;
    }

double controlEligibilityForTraj_fd(const unsigned long long int *x, const unsigned long long int *u, const unsigned long long int *previousU)
    {
    return 1.0;
    }

double constraintsXUY_fd(const unsigned long long int *x, const unsigned long long int *u)
    {
    return -PLUS_INF;
    }

double constraintsXV_tych(const double *x, const double *v)
    {
    return 0.0;
    }

void dynamics_hybrid_c(const double *x, const unsigned long long int *xd, const double *u, double *image)
    {

    }
void dynamics_hybrid_d(const double *x, const unsigned long long int *xd, const unsigned long long int *u, unsigned long long int *image)
    {

    }
double constraintsXU_hybrid(const double *x, const unsigned long long int *xd, const double *u, const unsigned long long int *ud)
    {
    return 0.0;
    }

void resetMap_hybrid(const double * xc, const unsigned long long int* xd, const unsigned long long int* resetControl, double * imagec, const unsigned long long int* imaged)
    {

    }

double constraintsX_hybrid(const double *x, const unsigned long long int *xd)
    {
    return 0.0;
    }
void jacobian_hybrid(const double *x, const unsigned long long int *xd, const double *u, double **jacob)
    {

    }
void localDynBounds_hybrid(const double *x, const unsigned long long int *xd, double *res)
    {

    }
