#include <numeric> //Pour std::accumulate

#include "../include/TychePicker.h"

std::mt19937_64 TychePicker::mt{};

TychePicker::TychePicker(const TrajectoryParametersManager *tpm, SysDyn *sysDyn) :
    // Seed par défaut
    methods(nullptr),
    sysDyn(sysDyn),
    tycheValueIndexes(nullptr),
    nonDeterministicAxes(nullptr),
    maxNbRerolls(nullptr),
    cumulativeDistribution(nullptr),
    cumulativeDistributionValues(nullptr),
    density(nullptr),
    densityValues(nullptr),
    trajIndex(tpm->getTrajectoryIndex()),
    dimTy(sysDyn->getDimTy()),
    nbNonDeterministicAxes(0),
    maxTotalNbRerolls(0)
{
    methods = new tychePickMethod[dimTy];
    nonDeterministicAxes = new int[dimTy];
    tycheValueIndexes = new unsigned long long int[dimTy];
    const trajectoryParams * tp = tpm->getTrajectoryParameters();
    const tycheParams *tycheParamsArray = tp->TYCHE_PARAMS;

    unsigned long long int *nbPointsTy = sysDyn->getNbPointsTy();

    densityValues = new double*[dimTy]();
    cumulativeDistributionValues = new double*[dimTy]();
    maxNbRerolls = new int[dimTy];
        
    // Valeur identique pour chaque tycheParams
    userTyche = tycheParamsArray[0].USER_TYCHE;
    density = tycheParamsArray[0].PROBABILITY_DENSITY;
    cumulativeDistribution = tycheParamsArray[0].CUMULATIVE_DISTRIBUTION;
    
    for (int i = 0; i < dimTy; ++i) {
        maxNbRerolls[i] = tycheParamsArray[i].MAX_NB_REROLLS;
            
        switch (tycheParamsArray[i].TYCHE_DISTRIBUTION) {
        case UNIFORM:
            methods[i] = &TychePicker::pickUniformTyche;
            
            nonDeterministicAxes[nbNonDeterministicAxes] = i;
            nbNonDeterministicAxes++;
            maxTotalNbRerolls += maxNbRerolls[i];
            break;
        case CONSTANT:
            methods[i] = &TychePicker::returnConstantTycheValue;
            // On cherche l'indice de la coordonnée dans les CONTROL_TY_GRID_POINTS[i]
            // telle que la valeur de l'utilisateur est la plus proche possible
            tycheValueIndexes[i] = closestTychePointIndex(tycheParamsArray[i].CONSTANT_TYCHE_VALUE, i);
            break;
        case CUSTOM_DETERMINED:
            methods[i] = &TychePicker::pickUserTyche;
            break;
        case CONSTANT_CUMULATIVE_DISTRIBUTION:
            methods[i] = &TychePicker::returnCustomlyDistributedTyche;            
            cumulativeDistributionValues[i] = new double[nbPointsTy[i]]{};
            getCumulativeDistribution(cumulativeDistributionValues[i], nbPointsTy[i], nullptr, 0.0, i);
            
            nonDeterministicAxes[nbNonDeterministicAxes] = i;
            nbNonDeterministicAxes++;
            maxTotalNbRerolls += maxNbRerolls[i];
            break;
        case CUMULATIVE_DISTRIBUTION:
            methods[i] = &TychePicker::pickCumulativeDistributionTyche;            
            cumulativeDistributionValues[i] = new double[nbPointsTy[i]]{};
            
            nonDeterministicAxes[nbNonDeterministicAxes] = i;
            nbNonDeterministicAxes++;
            maxTotalNbRerolls += maxNbRerolls[i];
            break;
        case CONSTANT_PROBABILITY_DENSITY:
            methods[i] = &TychePicker::returnCustomlyDistributedTyche;            
            densityValues[i] = new double[nbPointsTy[i]]{};
            cumulativeDistributionValues[i] = new double[nbPointsTy[i]]{};
            getDensity(densityValues[i], nbPointsTy[i], nullptr, 0.0, i);
            getCumulativeDistributionFromDensity(densityValues[i], cumulativeDistributionValues[i], nbPointsTy[i]);
            
            nonDeterministicAxes[nbNonDeterministicAxes] = i;
            nbNonDeterministicAxes++;
            maxTotalNbRerolls += maxNbRerolls[i];
            break;        
        case PROBABILITY_DENSITY:
            methods[i] = &TychePicker::pickDensityTyche;            
            densityValues[i] = new double[nbPointsTy[i]]{};
            cumulativeDistributionValues[i] = new double[nbPointsTy[i]]{};
            
            nonDeterministicAxes[nbNonDeterministicAxes] = i;
            nbNonDeterministicAxes++;
            maxTotalNbRerolls += maxNbRerolls[i];
            break;
        }
    }

    if (tp->SEED_LENGTH > 0) {
        std::seed_seq seq(tp->SEED, tp->SEED + tp->SEED_LENGTH);
        mt.seed(seq);
    }
    
    mt.discard(100);
}

unsigned long long int TychePicker::pickTyche(const double *xCoordsDouble, double t) {

    unsigned long long int *nbPoints = sysDyn->getNbPointsTy();
    unsigned long long int index = 0;
    
    for (int i = 0; i < dimTy; ++i) {
        const unsigned long long int pointIndex = ((this->*(methods[i])))(xCoordsDouble, t, i);
        index = index*nbPoints[i] + pointIndex;
    }

    // Vérifier si la contrainte état-tyché est validée
    double **tycheCoords = sysDyn->getTychCoords();    
    // Si le nombre d'axes non-déterministe est 0 et que la contrainte n'est pas vérifiée
    // alors on ne peut pas continuer
    if (nbNonDeterministicAxes == 0) {
        if (sysDyn->constraintsXV_tych(xCoordsDouble, tycheCoords[index]) >= PLUS_INF) {
            spdlog::error("Deterministic tyche doesn't meet state-tyche requirements");
            index = sysDyn->getTotalNbPointsTy();
        }
    }
    else {
        int rerollIndex = 0;
        int totalNbRerolls = 0;
        int *nbRerolls = new int[dimTy]{0};

        int i = nonDeterministicAxes[rerollIndex];

        bool isValid = false;
        
        do {
            isValid = (sysDyn->constraintsXV_tych(xCoordsDouble, tycheCoords[index]) < PLUS_INF);
            while (!isValid && nbRerolls[i] < maxNbRerolls[i]) {
        
                int rerolledPointIndex = ((this->*(methods[i])))(xCoordsDouble, t, i);
                nbRerolls[i]++;
                totalNbRerolls++;
                
                /* = $\prod_{k=i+1}^{dimTy}{nbPoints[i]}$
                   J'ai supposé que dans la plupart des cas, on aurait pas à faire ce calcul.
                   Il se fait uniquement pour une fonction constraintsXV_tych de l'utilisateur qui,
                   si elle existe, devrait dans une majorité des cas ne pas contraindre le tyché.
                   Je n'ai donc pas cherché à stocker une information inutile dans une majorité des cas */
                const unsigned long long int nbPointsProd = std::accumulate(nbPoints+i+1, nbPoints+dimTy, 1,
                                                                            [](unsigned long long a, unsigned long long b){return a*b;});

                // On calcule la coordonnée d'indice i
                const unsigned long long int oldPointIndex = (index / nbPointsProd) % nbPoints[i];
                // Et on remplace l'ancien indice par le nouveau
                index += (rerolledPointIndex - oldPointIndex)*nbPointsProd;

                isValid = (sysDyn->constraintsXV_tych(xCoordsDouble, tycheCoords[index]) < PLUS_INF);
                
                rerollIndex = (rerollIndex + 1) % nbNonDeterministicAxes;
                i = nonDeterministicAxes[rerollIndex];
            }
            rerollIndex = (rerollIndex + 1) % nbNonDeterministicAxes;            
            i = nonDeterministicAxes[rerollIndex];
        } while (!isValid && totalNbRerolls < maxTotalNbRerolls);

        if (totalNbRerolls == maxTotalNbRerolls && !isValid) {
            spdlog::error("Unable to find tyche validating constraintsXV_tych after {} rerolls.", totalNbRerolls);
            index = sysDyn->getTotalNbPointsTy();
        }
    }
    
    return index;
}

unsigned long long int TychePicker::returnConstantTycheValue(const double *, double, int coordIndex) {
    return tycheValueIndexes[coordIndex];
}

unsigned long long int TychePicker::pickUserTyche(const double *xCoordsDouble, double t, int coordIndex) {
    double userValue = userTyche(xCoordsDouble, t, coordIndex, trajIndex);
    return closestTychePointIndex(userValue, coordIndex);
}

unsigned long long int TychePicker::pickUniformTyche(const double *, double, int coordIndex) {
    unsigned long long int nbPoints = sysDyn->getNbPointsTy()[coordIndex];
    std::uniform_int_distribution<unsigned long long int> distrib(0, nbPoints - 1);
    // Si constraintsXV, alors faire une permutation aléatoire d'un vecteur de 0 à nbPointsTy[dim] - 1,
    // et parcourir jusqu'à ce que ça marche...
    return distrib(mt);
}

unsigned long long int TychePicker::pickCumulativeDistributionTyche(const double *xCoordsDouble, double t, int coordIndex) {
    unsigned long long int nbPoints = sysDyn->getNbPointsTy()[coordIndex];
    getCumulativeDistribution(cumulativeDistributionValues[coordIndex], nbPoints, xCoordsDouble, t, coordIndex);
    return returnCustomlyDistributedTyche(xCoordsDouble, t, coordIndex);
}

unsigned long long int TychePicker::pickDensityTyche(const double *xCoordsDouble, double t, int coordIndex) {
    unsigned long long int nbPoints = sysDyn->getNbPointsTy()[coordIndex];
    getDensity(densityValues[coordIndex], nbPoints, xCoordsDouble, t, coordIndex);
    getCumulativeDistributionFromDensity(densityValues[coordIndex], cumulativeDistributionValues[coordIndex], nbPoints);
    return returnCustomlyDistributedTyche(xCoordsDouble, t, coordIndex);
}

unsigned long long int TychePicker::returnCustomlyDistributedTyche(const double *xCoordsDouble, double t, int coordIndex) {
    const unsigned long long int nbPoints = sysDyn->getNbPointsTy()[coordIndex];
    std::uniform_real_distribution<double> distrib(0.0, 1.0);

    double uniformReal = distrib(mt);

    unsigned long long int i = 0;
    while (i < nbPoints - 1 && cumulativeDistributionValues[coordIndex][i] <= uniformReal) {
        // Si constraintsXV, on suppose que l'on normalise parmi les choix restants en tirant un nouveau nombre ?
        ++i;
    }
    return i;
}

void TychePicker::getDensity(double *densityValues, unsigned long long int nbPoints, const double *xCoords, double t, int coordIndex) {
    for (unsigned long long int i = 0; i < nbPoints; ++i) {
        double val = getValueOfIndex(coordIndex, i);
        densityValues[i] = density(val, coordIndex, xCoords, t, trajIndex);
    }
}

void TychePicker::getCumulativeDistributionFromDensity(double *density, double *cumulative, unsigned long long int nbPoints) {

    double sum = 0;
    for (unsigned long long int i = 0; i < nbPoints; ++i) {
        sum += density[i];
        cumulative[i] = sum;
    }
    if (sum > 1.0) {
        spdlog::warn("Cumulative distribution gotten through density is greater than 1");
        spdlog::warn("The greatest value given by the distribution will have probabilty of {}", 1.0 - cumulative[nbPoints - 2]);
    }
}

void TychePicker::getCumulativeDistribution(double *cumulative, unsigned long long int nbPoints,
                                            const double *xCoords, double t, int coordIndex) {
    for (unsigned long long int i = 0; i < nbPoints; ++i) {
        double val = getValueOfIndex(coordIndex, i);
        cumulative[i] = cumulativeDistribution(val, coordIndex, xCoords, t, trajIndex);
    }
    if (cumulative[nbPoints-1] > 1.0) {
        spdlog::warn("Cumulative distribution gotten through user funciton greater than 1");
        spdlog::warn("The greatest value given by the distribution will have probabilty of {}", 1.0 - cumulative[nbPoints - 2]);
    }
}

double TychePicker::getValueOfIndex(int coordIndex, int pointIndex) {
    unsigned long long int nbPoints = sysDyn->getNbPointsTy()[coordIndex];
    double limInf = sysDyn->getLimInfTy()[coordIndex];
    double step = sysDyn->getStepTy()[coordIndex];
    return limInf + pointIndex*step;
}

unsigned long long int TychePicker::closestTychePointIndex(double value, int coordIndex) {

    double limInfTy = sysDyn->getLimInfTy()[coordIndex];
    double stepTy = sysDyn->getStepTy()[coordIndex];
    unsigned long long int nbPointsTy = sysDyn->getNbPointsTy()[coordIndex];
    
    double unroundedIndex = (value - limInfTy)/stepTy;
    if (unroundedIndex < 0.0 || unroundedIndex >= nbPointsTy) {
        double limSupTy = sysDyn->getLimSupTy()[coordIndex];
        spdlog::warn("Tyche value n°{} outside of control domain", coordIndex+1);
        spdlog::warn("Clamping coordinate of value {} (index n°{}) into [{} ,{}]",
                     value, coordIndex+1, limInfTy, limSupTy);
        unroundedIndex = (unroundedIndex < 0.0) ? 0.0 : (nbPointsTy - 1);
    }    
    return lround(unroundedIndex);
}

TychePicker::~TychePicker() noexcept {
    delete [] methods;
    delete [] tycheValueIndexes;
    delete [] nonDeterministicAxes;
    
    for (int i = 0; i < dimTy; ++i) {
        delete [] cumulativeDistributionValues[i];
    }
    delete [] cumulativeDistributionValues;
    for (int i = 0; i < dimTy; ++i) {
        delete [] densityValues[i];
    }
    delete [] densityValues;
}
