#ifndef TYCHEPICKER_H
#define TYCHEPICKER_H

#include <random>

#include "ParametersManager.h"
#include "SysDyn.h"

class TychePicker
{
public:
    TychePicker(const TrajectoryParametersManager *tpm, SysDyn *sysDyn);
    //! Copy constructor
    TychePicker(const TychePicker &other) = default;
    //! Move constructor
    TychePicker(TychePicker &&other) noexcept = default;
    //! Destructor
    ~TychePicker() noexcept;
    //! Copy assignment operator
    TychePicker& operator=(const TychePicker &other) = delete;
    //! Move assignment operator
    TychePicker& operator=(TychePicker &&other) noexcept = delete;

    unsigned long long int pickTyche(const double *xCoordsDouble, double t);
private:    
    
    using tychePickMethod = unsigned long long int (TychePicker::*)(const double *xCoordsDouble, double t, int i);

    // Statique pour ne pas réinitialiser la seed à chaque création de TychePicker (avoir de l'aléatoire par défaut entre deux trajectoires)
    // Son état peut être modifié dans le constructeur si une seed est donnée par les paramètres de trajectoire
    static std::mt19937_64 mt;
    
    unsigned long long int pickUserTyche(const double *xCoordsDouble, double t, int coordIndex);
    unsigned long long int pickUniformTyche(const double *xCoordsDouble, double t, int coordIndex);
    unsigned long long int returnConstantTycheValue(const double *xCoordsDouble, double t, int coordIndex);
    
    unsigned long long int returnCustomlyDistributedTyche(const double *xCoordsDouble, double t, int coordIndex);
    unsigned long long int pickCumulativeDistributionTyche(const double *xCoordsDouble, double t, int coordIndex);
    unsigned long long int pickDensityTyche(const double *xCoordsDouble, double t, int coordIndex);
    
    unsigned long long int closestTychePointIndex(double value, int coordIndex);

    double getValueOfIndex(int coordIndex, int pointIndex);
    
    void getDensity(double *density, unsigned long long int nbPoints,
                    const double *xCoordsDouble, double t, int coordIndex);
    void getCumulativeDistributionFromDensity(double *density, double *cumulative, unsigned long long int nbPoints);
    void getCumulativeDistribution(double *cumulative, unsigned long long int nbPoints,
                                   const double *xCoordsDouble, double t, int coordIndex);

    
    tychePickMethod *methods;
    SysDyn *sysDyn;
    unsigned long long int *tycheValueIndexes;
    
    int *nonDeterministicAxes;
    int *maxNbRerolls;
    
    userTyche_t userTyche;
    cumulativeDistribution_t cumulativeDistribution;
    double **cumulativeDistributionValues;
    probabilityDensity_t density;
    double **densityValues;
    
    int trajIndex;
    int dimTy;
    int nbNonDeterministicAxes;
    int maxTotalNbRerolls;
};

#endif /* TYCHEPICKER_H */
