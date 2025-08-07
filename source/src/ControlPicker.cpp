#include <dlfcn.h>
#include <numeric> // pour std::iota

#include "../include/ControlPicker.h"
#include "../include/ControlPickCriteria.h"
#include "../include/TrajectorySimulation.h"
#include "../include/TrajectoryPointsSimulation.h"
#include "../include/ControlPickStrategy.h"

ControlPicker::ControlPicker(SysDyn *sysDyn, TrajectoryParametersManager *tpm, indexSorter_t *sorterPtr) :
    sysDyn(sysDyn),
    strategies{}
{
    const trajectoryParams *tp = tpm->getTrajectoryParameters();
    trajIndex = tpm->getTrajectoryIndex();    

    unsigned long long int nbTotalC = sysDyn->getTotalNbPointsC();    
    preferedControlIndexes = new int[nbTotalC];
    std::iota(preferedControlIndexes, preferedControlIndexes+nbTotalC, 0);
    controlWeight = tp->CONTROL_WEIGHT;
    sortIndexes = sorterPtr;
}    

OptionalCu ControlPicker::pickControl(Trajectory &traj, TrajectoryPoints &trajDiscrete,
                                      double rho, StrategyIndexBitFlag &flag) {
    return pickControlFromSubPickerList(0, traj, trajDiscrete, rho, flag);
}

OptionalCu ControlPicker::pickControlFromSubPickerList(
    int subListStartIndex, Trajectory &traj, TrajectoryPoints &trajDiscrete,
    double rho, StrategyIndexBitFlag &flag) {
    
    OptionalCu cuEither(UNSATISFIED_STRATEGY);
    OptionalCu prevCuEither = cuEither;
    
    const int nbStrats = strategies.size();

    // Les pickers définissent eux même l'ordre des contrôles
    // Donc, on ne veut aucun tri de l'ordre des préférences
    // par défaut
    *sortIndexes = &noSort;

    // asSimulation pour que les stratégies ne modifient pas les trajectoires
    TrajectorySimulation trajSimu = traj.asSimulation();
    TrajectoryPointsSimulation trajDiscreteSimu = trajDiscrete.asSimulation();
    
    for (int i = subListStartIndex; i < nbStrats; ++i) {
        ControlPickCriteria criteria(this, sysDyn, &trajSimu, &trajDiscreteSimu, i, flag, rho);
        prevCuEither = cuEither;
        cuEither = strategies[i]->pickControl(cuEither, criteria);

        flag |= ((cuEither != prevCuEither) ? 1<<i : 0);
    }

    return cuEither;
}

ControlPicker *ControlPicker::addPicker(ControlPickStrategy *strategy) {

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

ControlPicker *ControlPicker::addUserPicker(TrajectoryParametersManager *params, const std::string &strategyName) {

    void *modelHandle = params->getModelHandle();

    std::string factoryFunctionName("new" + strategyName);
    int strategyIndex = strategies.size();
    
    UserPickStrategy *(*createStrategy)(int strategyIndex, const TrajectoryParametersManager *) =
        (UserPickStrategy *(*)(int, const TrajectoryParametersManager *)) dlsym(modelHandle, factoryFunctionName.c_str());
    if (createStrategy == nullptr) {
        spdlog::error("Could not load {} control picker, no {} function defined | {}",
                      strategyName, factoryFunctionName, dlerror()
            );
        return this;
    }
    UserPickStrategy *userStrategy = createStrategy(strategyIndex, params);
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

std::vector<std::string> ControlPicker::getStrategyNames() const {
    std::vector<std::string> names;
    for (ControlPickStrategy *strategy : strategies) {
        names.push_back(strategy->getName());
    }
    return names;
}

int *ControlPicker::getPreferedControlIndexes() {
    return preferedControlIndexes;
}

int *ControlPicker::sortPreferedControlIndexes(const double *x, double t, int strategyIndex) {
    unsigned long long int nbCTotal = sysDyn->getTotalNbPointsC();
    double **controlCoords = sysDyn->getControlCoords();
    (*sortIndexes)(preferedControlIndexes, x, controlCoords, nbCTotal, t, controlWeight, trajIndex, strategyIndex);
    return preferedControlIndexes;
}

void ControlPicker::setIndexSorter(indexSorter_t sorter) {    
    *sortIndexes = sorter;
}

int ControlPicker::getTrajectoryIndex() const {
    return trajIndex;
}

int ControlPicker::getNbStrategies() const {
    return strategies.size();
}

ControlPicker::~ControlPicker() {
    std::for_each(strategies.begin(), strategies.end(), [](ControlPickStrategy *s) {
        delete s;
    });
    delete[] preferedControlIndexes;
}
