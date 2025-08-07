#ifndef BUBBLE_H
#define BUBBLE_H

#include "SysDyn.h"
#include "Grid.h"
#include "ParametersManager.h"

class Bubble {
private:
    typedef unsigned long long int ull;
public:
    Bubble(const TrajectoryParametersManager *params, const Grid *grid, int strategyIndex = 0);
   ~Bubble();
    void setCenter(const ull *newCenter);
    void setCenter(int index, ull value);

    /*!
     * Renvoie si un point de grille repéré par ses coordonnées de grille
     * est compris dans une bulle définie par un centre et des distances
     * en nombre de points à ce centre sur chaque dimension
     */
    bool isInBubble(const ull *point) const;
    
    bool isOnBorder(const ull *point);
    /*!
     * Renvoie si un point de grille non contenu dans l'ensemble de viabilité
     * se trouve dans une bulle définie par un centre et des distances en
     * nombre de points à ce centre sur chaque dimension
     */
    bool isTouchingBoundaryLexicographic();
    bool isTouchingBoundaryInwards();
    
    void getClosestPointOnBorder(ull *point);
    
    int getStrategyIndex();
private:
    const Grid *grid;
    ull *center;
    
    double *radius;
    int *radiusCopy;
    int maxRadius;

    ull *offsetCoords;
    int *centerOffsets;
    int dim;

    int trajIndex, strategyIndex;
    neighborValidator_t isValidNeighbor;
    
    bool isOutsideViabilityKernel(const ull *point);

    // Fonctions itérant sur l'ensemble des points de la bulle

    // Dans l'ordre lexicographique des coordonnées
    ull *getNextNeighborLexicographic(ull *offsetCoords);
    ull *getNextNeighborLexicographicAux(ull *offsetCoords, const ull *center, int *centerOffsets, const int dim, const int *radius);

    // En partant des points les plus extérieurs de la bulle et se rapprochant du centre
    ull *getNextNeighborInwards(ull *offsetCoords);
    ull *getNextNeighborInwardsAux(ull *offsetCoords, const ull *center, int *centerOffsets, const int dim, const int *radius);

    typedef ull * (Bubble::*getNextNeighborMemberPtr)(ull *offsetCoords);
    typedef bool (Bubble::*satisfiesMemberPtr)(const ull *);    
    ull *getFirstPointSatisfying(ull *offsetCoords,
                                 getNextNeighborMemberPtr getNextNeighbor,
                                 satisfiesMemberPtr satisfies);

    double squaredDistanceBetween(const double *p1, const double *p2, int dim) const;
};

#endif /* BUBBLE_H */
