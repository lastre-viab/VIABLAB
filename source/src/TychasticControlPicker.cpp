#include <dlfcn.h>
#include <numeric> // pour std::iota

#include "../include/TychasticControlPicker.h"
#include "../include/TychasticControlPickCriteria.h"
#include "../include/TychasticTrajectorySimulation.h"
#include "../include/TrajectoryPointsSimulation.h"
#include "../include/TychasticControlPickStrategy.h"

TychasticControlPicker::TychasticControlPicker(SysDyn *sysDyn, TrajectoryParametersManager *tpm, indexSorter_t *sorterPtr) :
    strategies{},
    sysDyn(sysDyn),
    tp(*tpm->getTrajectoryParameters())
{
    trajIndex = tpm->getTrajectoryIndex();    

    unsigned long long int nbTotalC = sysDyn->getTotalNbPointsC();    
    preferedControlIndexes = new int[nbTotalC];
    std::iota(preferedControlIndexes, preferedControlIndexes+nbTotalC, 0);
    controlWeight = tp.CONTROL_WEIGHT;
    sortIndexes = sorterPtr;
}    

OptionalCu TychasticControlPicker::pickControl(TychasticTrajectory &traj, TrajectoryPoints &trajDiscrete,
                                               double rho, StrategyIndexBitFlag &flag) {
    return pickControlFromSubPickerList(0, traj, trajDiscrete, rho, flag);
}

OptionalCu TychasticControlPicker::pickControlFromSubPickerList(
    int subListStartIndex, TychasticTrajectory &traj, TrajectoryPoints &trajDiscrete,
    double rho, StrategyIndexBitFlag &flag) {
    
    OptionalCu cuEither(UNSATISFIED_STRATEGY);
    OptionalCu prevCuEither = cuEither;
    
    const int nbStrats = strategies.size();

    // Les pickers définissent eux même l'ordre des contrôles
    // Donc, on ne veut aucun tri de l'ordre des préférences
    // par défaut
    *sortIndexes = &noSort;

    // asSimulation pour que les stratégies ne modifient pas les trajectoires
    TychasticTrajectorySimulation trajSimu = traj.asSimulation();
    TrajectoryPointsSimulation trajDiscreteSimu = trajDiscrete.asSimulation();
    
    for (int i = subListStartIndex; i < nbStrats; ++i) {
        TychasticControlPickCriteria criteria(this, tp, sysDyn, &trajSimu, &trajDiscreteSimu, i, flag, tycheIndex, rho);
        prevCuEither = cuEither;
        cuEither = strategies[i]->pickControl(cuEither, criteria);

        flag |= ((cuEither != prevCuEither) ? 1<<i : 0);
    }

    return cuEither;
}

TychasticControlPicker *TychasticControlPicker::addPicker(TychasticControlPickStrategy *strategy) {

    constexpr int maxNbStrategies = 8*sizeof(StrategyIndexBitFlag) - 1;
    
    if (strategies.size() >= maxNbStrategies) {
        spdlog::error("Possibilty to have strategy list longer than {} strategies not supported. Sorry.", maxNbStrategies);
        spdlog::error("All strategies past index {} will be ignored. Try creating your own personalized strategy.", maxNbStrategies - 1);
    }
    else {
        strategies.push_back(strategy);
    }
    return this;
}

TychasticControlPicker *TychasticControlPicker::addUserPicker(TrajectoryParametersManager *params, const std::string &strategyName) {

    void *modelHandle = params->getModelHandle();

    std::string factoryFunctionName("new" + strategyName);
    int strategyIndex = strategies.size();
    
    TychasticUserPickStrategy *(*createStrategy)(int strategyIndex, const TrajectoryParametersManager *) =
        (TychasticUserPickStrategy *(*)(int, const TrajectoryParametersManager *)) dlsym(modelHandle, factoryFunctionName.c_str());
    if (createStrategy == nullptr) {
        spdlog::error("Could not load {} control picker, no {} function defined | {}",
                      strategyName, factoryFunctionName, dlerror()
            );
        return this;
    }
    TychasticUserPickStrategy *userStrategy = createStrategy(strategyIndex, params);
    if (userStrategy == nullptr) {
        spdlog::error("Could not create {} control picker, function {} returned nullptr",
                      strategyName, factoryFunctionName
            );
        return this;
    }

    userStrategy->setName(strategyName);
    addPicker(userStrategy);        
    return this;
}

std::vector<std::string> TychasticControlPicker::getStrategyNames() const {
    std::vector<std::string> names;
    for (TychasticControlPickStrategy *strategy : strategies) {
        names.push_back(strategy->getName());
    }
    return names;
}

int *TychasticControlPicker::getPreferedControlIndexes() {
    return preferedControlIndexes;
}

int *TychasticControlPicker::sortPreferedControlIndexes(const double *x, double t, int strategyIndex) {
    unsigned long long int nbCTotal = sysDyn->getTotalNbPointsC();
    double **controlCoords = sysDyn->getControlCoords();
    (*sortIndexes)(preferedControlIndexes, x, controlCoords, nbCTotal, t, controlWeight, trajIndex, strategyIndex);
    return preferedControlIndexes;
}

void TychasticControlPicker::setIndexSorter(indexSorter_t sorter) {    
    *sortIndexes = sorter;
}

int TychasticControlPicker::getTrajectoryIndex() const {
    return trajIndex;
}

int TychasticControlPicker::getNbStrategies() const {
    return strategies.size();
}

void TychasticControlPicker::setTycheIndex(int index) {
    tycheIndex = index;
}

TychasticControlPicker::~TychasticControlPicker() {
    std::for_each(strategies.begin(), strategies.end(), [](TychasticControlPickStrategy *s) {
        delete s;
    });
    delete[] preferedControlIndexes;
}
