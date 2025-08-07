#include "../include/TrajectoryHelpers.h"
#include "../include/WeakDeclarations.h"
#include "../include/DefaultValues.h"

#include <random>
#include <algorithm>

void sortControlIndexesByWeight(int *indexes, const double *x, double *const * const us, long long unsigned int nbCTotal, double normalizedTime, controlWeight_t controlWeight, int trajIndex, int stratIndex) {
    
    std::iota(indexes, indexes+nbCTotal, 0);
    std::sort(indexes, indexes+nbCTotal, [=](int index1, int index2) {
        return controlWeight(x, us[index1], normalizedTime, trajIndex, stratIndex) > controlWeight(x, us[index2], normalizedTime, trajIndex, stratIndex);
    });
}

void shuffleControlIndexes(int *indexes, const double *, double *const *const, long long unsigned int nbCTotal, double, controlWeight_t, int, int) {
    
    // Source d'aléatoire à possiblement modifier
    static std::mt19937_64 mt{};
    
    std::shuffle(indexes, indexes+nbCTotal, mt);
}

void noSort(int *, const double *, double *const *const, long long unsigned int, double, controlWeight_t, int, int) {}

bool belowInfiniteDistance(const unsigned long long int *coord1, const unsigned long long int *coord2, const double *radius, int, int) {
    // coord1 et coord2 sont déjà pioché par une distance infinie, donc c'est toujours vrai
    return true;
}

bool belowEuclideanDistance(const unsigned long long int *coord1, const unsigned long long int *coord2, const double *radius, int, int) {

    const int dim = gp.DIM;
    
    // suppose que radius[i] == radius[j] pour tous les i et j
    
    long long int squaredDist = 0.0;
    for (int i = 0; i < dim; ++i) {

        long long int square = coord1[i] - coord2[i];
        
        squaredDist += square * square;
    }
    return squaredDist <= radius[0]*radius[0];
}

bool belowEllipticDistance(const unsigned long long int *coord1, const unsigned long long int *coord2, const double *radius, int, int) {

    const int dim = gp.DIM;
    double cartesianEquation = 0.0;    
    
    for (int i = 0; i < dim; ++i) {

        long long int diff = coord1[i] - coord2[i];
        double coefficient = ((double) diff) / radius[i];
        
        cartesianEquation += coefficient * coefficient;
    }
    return cartesianEquation <= 1.0;
}
