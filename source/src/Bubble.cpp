#include "../include/Bubble.h"

Bubble::Bubble(const TrajectoryParametersManager *params, const Grid *grid, int stratIndex) :
    grid(grid)
{
    const trajectoryParams *tp = params->getTrajectoryParameters();
    const gridParams *gp = params->getGridParameters();

    dim = gp->DIM;
    // centerOffsets peut contenir une information supplémentaire
    // d'où le +1
    centerOffsets = new int[dim+1];
    center = new ull[dim];
    offsetCoords = new ull[dim];
    radius = tp->BUBBLE_RADIUS;

    // La copie évite de devoir modifier le tableau radius lors d'algorithmes récursifs
    // (et, accessoirement, de faire une conversion en entier à chaque fois)
    radiusCopy = new int[dim];
    for (int i = 0; i < dim; ++i) {
        radiusCopy[i] = (int) radius[i];
    }
    
    maxRadius = (int) *std::max_element(radius, radius+dim);
    
    isValidNeighbor = tp->IS_VALID_NEIGHBOR;
    trajIndex = params->getTrajectoryIndex();
    strategyIndex = stratIndex;
}

void Bubble::setCenter(const ull *newCenter) {
    std::copy(newCenter, newCenter+dim, center);
}

void Bubble::setCenter(int i, ull v) {
    center[i] = v;
}

int Bubble::getStrategyIndex() {
    return strategyIndex;
}

Bubble::ull *Bubble::getFirstPointSatisfying(
    ull *offsetCoords,
    getNextNeighborMemberPtr getNextNeighbor,
    satisfiesMemberPtr isSatisfying) {
    
    for (int i = 0; i < dim; ++i) {
        centerOffsets[i] = -radiusCopy[i];
        offsetCoords[i] = center[i] - radiusCopy[i];
    }

    // Compteur utile pour certaines versions de getNextNeighbor
    centerOffsets[dim] = 0;

    ull *foundNeighbor = offsetCoords;
    // Il se peut que le point de départ où toutes les coordonnées
    // sont à leur minimum ne soit pas valide
    if (!isValidNeighbor(center, foundNeighbor, radius, trajIndex, strategyIndex)) {
        foundNeighbor = std::invoke(getNextNeighbor, this, foundNeighbor);
    }
    
    while (foundNeighbor && !std::invoke(isSatisfying, this, foundNeighbor)) {
        foundNeighbor = std::invoke(getNextNeighbor, this, foundNeighbor);
    }
    
    return foundNeighbor;
}

bool Bubble::isInBubble(const ull *currentPos) const {
    int i = 0;
    bool inBubble = true;
    long long int dist;
    while (inBubble && i < dim) {
        dist = currentPos[i] - center[i];
        inBubble &= ( -radius[i] <= dist && dist <= radius[i] );
        i++;
    }
    // isValidNeighbor suppose que l'on est au moins sûr d'être
    // dans le voisinage de Moore, d'où la vérification inBubble
    return inBubble && isValidNeighbor(currentPos, center, radius, trajIndex, strategyIndex);
}

Bubble::ull *Bubble::getNextNeighborLexicographic(ull *offsetCoords) {
    ull *res;
    do {
        res = getNextNeighborLexicographicAux(offsetCoords, center, centerOffsets, dim, radiusCopy);
    } while (res != nullptr && !isValidNeighbor(offsetCoords, center, radius, trajIndex, strategyIndex));
    return res;
}

Bubble::ull *Bubble::getNextNeighborLexicographicAux(ull *offsetCoords, const ull *center, int *centerOffsets, const int dim, const int *radius) {

    // centerOffsets est préparé pour l'itération suivante et offsetCoords vaut le résultat de l'itération précédente
    // On veut pouvoir utiliser cette méthode pour itérer sur tous les points
    // dont le premier, donc on met seulement à jour les compteurs pour l'itération suivante
    
    int i = dim - 1;
    centerOffsets[i]++;
    offsetCoords[i] = center[i] + centerOffsets[i];

    // Gestion du dépassement dans une coordonnée, il faut faire une "retenue"
    while (i > 0 && centerOffsets[i] > radiusCopy[i]) {
        // Si on a dépassé, on renvient sur le minimum de la coordonnée cette itération
        centerOffsets[i] = -radiusCopy[i];
        offsetCoords[i] = center[i] + centerOffsets[i];
        // Et on applique la retenue sur la dimension suivante,
        // qui pourra aussi entrainer une retenue, d'où le while
        i--;
        centerOffsets[i]++;
        offsetCoords[i] = center[i] + centerOffsets[i];
    }   

    // Si la retenue s'est propagée jusqu'à la dernière dimension
    // et que celle-ci est à son maximum, alors on a itéré sur tous les voisins
    if (centerOffsets[0] > radius[0]) {
        return nullptr;
    }
    else {
        return offsetCoords;
    }
}

Bubble::ull *Bubble::getNextNeighborInwards(ull *offsetCoords) {
    ull *res;
    int *iteration = centerOffsets + dim;
    do {
        res = getNextNeighborInwardsAux(offsetCoords, center, centerOffsets, dim, radiusCopy);
        if (!res) {
            if (*iteration < maxRadius) {
                // Tous les points pour le rayon donné on été trouvé
                // On diminue le rayon de 1 pour obtenir les points intérieurs (inwards)
                (*iteration)++;
                for (int i = 0; i < dim; ++i) {
                    radiusCopy[i] = max(0, radiusCopy[i] - 1);
                    centerOffsets[i] = -radiusCopy[i];
                    offsetCoords[i] = center[i] + centerOffsets[i];
                }
                res = offsetCoords;
            }
            // Il n'y a plus de voisin suivant, radiusCopy doit valoir radius pour l'itération suivante
            else {
                for (int i = 0; i < dim; ++i) {
                    radiusCopy[i] = (int) radius[i];
                }
            }
        }
    } while (res != nullptr && !isValidNeighbor(offsetCoords, center, radius, trajIndex, strategyIndex));
    return res;
}

Bubble::ull *Bubble::getNextNeighborInwardsAux(
    ull *offsetCoords, const ull *center, int *centerOffsets, const int dim, const int *radius) {
    
    ull *res = nullptr;

    // Cas de base : on ne veut que les extrema de cette dimension
    if (dim == 1) {
        // On est au maximum, il n'y a plus d'extrema disponibles
        if (centerOffsets[0] == radius[0]) {
            centerOffsets[0] = -radius[0];
            res = nullptr;
        }
        else {            
            centerOffsets[0] = radius[0];
            res = offsetCoords;
        }
        offsetCoords[0] = center[0] + centerOffsets[0];
    }
    // Si on est au min ou au max de notre dimension, il faut itérer sur toutes les autres dimensions
    else if (std::abs(centerOffsets[0]) == radius[0]) {
        res = getNextNeighborLexicographicAux(offsetCoords, center, centerOffsets, dim, radius);
    }
    // Si on n'est pas au min ou au max, on itère récursivement sur les autres coordonnées
    else {
        res = getNextNeighborInwardsAux(offsetCoords + 1, center + 1, centerOffsets + 1, dim - 1, radius + 1);
        // On a fini l'itération sur les coordonnées en dimension - 1
        if (res == nullptr) {
            centerOffsets[0]++;
            offsetCoords[0] = center[0] + centerOffsets[0];            
        }
        res = offsetCoords;
    }
    return res;
}

bool Bubble::isOutsideViabilityKernel(const ull *point) {
    // Si un point n'est pas dans la grille, il n'est pas dans l'ensemble de viabilité non plus
    return !grid->isPointInGrid_fd(point) || !grid->isInSet(point);
}

bool Bubble::isOnBorder(const ull *currentPos) {

    // La position courante n'est pas nécessairement dans le noyau de viabilité
    // La position ne peut être sur un bord du noyau que si elle même est dans le noyau
    if (!(grid->isPointInGrid_fd(currentPos) && grid->isInSet(currentPos))) {
        return false;
    }

    // On veut que la méthode ait un paramètre const
    // Bien que currentPos soit modifié, la valeur finale
    // est identique mais cela nous empêche de marquer la méthode comme const
    // On fait donc une copie de tableau, en espérant que le compilateur comprenne
    // que cette copie peut être évitée
    //
    // Tant que cette implémentation ne cause pas de problème de lenteur
    // particulier, j'opterai pour la solution actuelle
    ull *currentPosCopy = new ull[dim];
    std::copy(currentPos, currentPos+dim, currentPosCopy);
    
    bool onBorder = false;
    int i = 0;
    while (!onBorder && i < dim) {
        currentPosCopy[i]++;
        onBorder |= !(grid->isPointInGrid_fd(currentPosCopy) && grid->isInSet(currentPosCopy));
        currentPosCopy[i] -= 2;
        onBorder |= !(grid->isPointInGrid_fd(currentPosCopy) && grid->isInSet(currentPosCopy));
        currentPosCopy[i]++;
        ++i;
    }

    delete[] currentPosCopy;
    
    return onBorder;
}

void Bubble::getClosestPointOnBorder(ull *point) {

    const int dim = grid->getDim();
    
    double *neighborCoordsAsDouble = new double[dim];
    double *pointAsDouble = new double[dim];
    grid->intCoordsToDoubleCoords(point, pointAsDouble);
    
    for (int i = 0; i < dim; ++i) {
        centerOffsets[i] = -radiusCopy[i];
        offsetCoords[i] = center[i] - radiusCopy[i];
    }
    centerOffsets[dim] = 0;

    ull *foundNeighbor = offsetCoords;
    
    if (!isValidNeighbor(center, foundNeighbor, radius, trajIndex, strategyIndex)) {
        foundNeighbor = getNextNeighborLexicographic(foundNeighbor);
    }
    double minDist = PLUS_INF;
    
    while (foundNeighbor) {
        grid->intCoordsToDoubleCoords(foundNeighbor, neighborCoordsAsDouble);
        double tmpDist = squaredDistanceBetween(pointAsDouble, neighborCoordsAsDouble, dim);
        if (tmpDist < minDist && isOnBorder(foundNeighbor)) {
            minDist = tmpDist;
            std::copy(foundNeighbor, foundNeighbor+dim, point);
        }
        foundNeighbor = getNextNeighborLexicographic(foundNeighbor);
    }   

    delete[] neighborCoordsAsDouble;
    delete[] pointAsDouble;    
}

double Bubble::squaredDistanceBetween(const double *point1, const double *point2, int dim) const {
    double sum = 0;
    for (int i = 0; i < dim; ++i) {
        double diff = point1[i] - point2[i];
        sum += diff*diff;

        // On a que
        // (x - y)^2 = x^2 + y^2 -2xy,
        // cette formule est plus stable numériquement pour des points p1 et p2 proches
        // Mais elle est (probablement) moins efficace à calculer
        //sum += point1[i]*point1[i] + point2[i]*point2[i] - 2*point1[i]*point2[i];
    }
    return sum;
}
                                   
bool Bubble::isTouchingBoundaryLexicographic() {
    return getFirstPointSatisfying(offsetCoords,
                                   &Bubble::getNextNeighborLexicographic,
                                   &Bubble::isOutsideViabilityKernel) != nullptr;
}

bool Bubble::isTouchingBoundaryInwards() {
    return getFirstPointSatisfying(offsetCoords,
                                   &Bubble::getNextNeighborInwards,
                                   &Bubble::isOutsideViabilityKernel) != nullptr;    
}


Bubble::~Bubble() {
    delete [] center;
    delete [] centerOffsets;
    delete [] radiusCopy;
    delete [] offsetCoords;
}
