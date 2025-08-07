#ifndef CONTROLPICKER_H
#define CONTROLPICKER_H

#include "ViabiBitSet.h"
#include "SysDyn.h"
#include "TrajectoryHelpers.h"
#include "Trajectory.h"
#include "TrajectoryPoints.h"
#include "Either.h"

class ControlPickStrategy;

enum ControlPickFailureReason {
    /* La raison pour laquelle aucun contrôle n'a été renvoyé
       est que la stratégie n'a pas trouvé de contrôle la satisfaisant */
    UNSATISFIED_STRATEGY,
    // D'autres valeurs pourront être ajoutées plus tard
};

struct pickedControl {
    double timeStep;
    int controlIndex;

    bool operator==(const pickedControl &other) const {
        return this->timeStep == other.timeStep
            && this->controlIndex == other.controlIndex;
    }
};

using OptionalCu = Either<pickedControl, ControlPickFailureReason>;

class ControlPicker {
public:
    using StrategyIndexBitFlag = Trajectory::BitFlag;
    
    OptionalCu pickControl(
        Trajectory &traj, TrajectoryPoints &trajDiscrete,
        double rho, StrategyIndexBitFlag &flag);

    OptionalCu pickControlFromSubPickerList(int startIndex,
        Trajectory &traj, TrajectoryPoints &trajDiscrete,
        double rho, StrategyIndexBitFlag &flag);

    ControlPicker operator=(const ControlPicker &) = delete;
    ControlPicker(const ControlPicker &) = delete;
    virtual ~ControlPicker();

    void setIndexSorter(indexSorter_t sortIndexes);
    int *getPreferedControlIndexes();
    int *sortPreferedControlIndexes(const double *x, double t, int strategyIndex);

    int getNbStrategies() const;
    int getTrajectoryIndex() const;

    std::vector<std::string> getStrategyNames() const;
    
protected:    
    ControlPicker(SysDyn *sysDyn, TrajectoryParametersManager *tpm, indexSorter_t *sorterPtr);
    ControlPicker *addUserPicker(TrajectoryParametersManager *params, const std::string &name);
    ControlPicker *addPicker(ControlPickStrategy *);
    
    SysDyn *sysDyn;
private:
    std::vector<ControlPickStrategy *> strategies;
    controlWeight_t controlWeight;
    indexSorter_t *sortIndexes;
    int *preferedControlIndexes;
    int trajIndex;
};

#endif /* CONTROLPICKER_H */
