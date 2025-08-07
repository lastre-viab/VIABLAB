#ifndef TRAJECTORYHELPERS_H
#define TRAJECTORYHELPERS_H


// À (possiblement) transformer en fonctions utilisant ControlPickCriteria
//
// Je ne suis pas certain de ce qui serait le mieux pour l'utilisateur
// entre un nombre important de paramètres ou une seule classe en paramètre
// ayant des méthodes pour récupérer les paramètres utiles à l'utilisateur
using controlWeight_t = double (*)(const double*, const double*, double, int, int);
using indexSorter_t = void (*)(int *indexes,
                               const double *x, double *const *const us, long long unsigned int nbCTotal, double t,
                               controlWeight_t controlWeight, int trajIndex, int stratIndex);
using neighborValidator_t = bool (*)(const unsigned long long int *, const unsigned long long int *,
                                     const double *, int, int);

using temporalControl_t = void (*)(double t, int trajIndex, int stratIndex, double *userControl);

using userTyche_t = double (*)(const double *x, double t, int coordIndex, int trajIndex);

using cumulativeDistribution_t = double (*)(double tycheValue, int coordIndex, const double *x, double t, int trajIndex);
using probabilityDensity_t = double (*)(double tycheValue, int coordIndex, const double *x, double t, int trajIndex);

/*
 * Fonctions servant à trier les contrôles.
 * La fonction utilisée est décidée selon le type de trajectoire demandée par
 * l'utilisateur.
 * Ce sont des fonctions dont le type est indexSorter_t
 */

void noSort(int *, const double *, double *const *const, long long unsigned int, double, controlWeight_t, int, int);

void shuffleControlIndexes(int *indexes, const double *, double *const *const us, long long unsigned int nbCTotal, double, controlWeight_t, int, int);

void sortControlIndexesByWeight(int *indexes, const double *x, double *const *const us, long long unsigned int nbCTotal, double normalizedTime, controlWeight_t controlWeight, int trajIndex, int stratIndex);

/*
 * Fonctions de filtrage des points d'une bulle.
 * La fonction utilisée est décidée selon la bulle demandée par l'utilisateur
 * lors d'une trajectoire "prudente"
 * Ce sont des fonctions dont le type est neighborValidator_t
 */
bool belowInfiniteDistance(const unsigned long long int *coord1, const unsigned long long int *coord2, const double *radius, int, int);

bool belowEllipticDistance(const unsigned long long int *coord1, const unsigned long long int *coord2, const double *radius, int, int);

bool belowEuclideanDistance(const unsigned long long int *coord1, const unsigned long long int *coord2, const double *radius, int, int);

#endif /* TRAJECTORYHELPERS_H */
