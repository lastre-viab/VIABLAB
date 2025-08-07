#ifndef TYCHASTICCONTROLPICKER_H
#define TYCHASTICCONTROLPICKER_H

#include "ViabiBitSet.h"
#include "SysDyn.h"
#include "TrajectoryHelpers.h"
#include "TychasticTrajectory.h"
#include "TrajectoryPoints.h"
#include "ControlPicker.h"
#include "TychePicker.h"

class TychasticControlPickStrategy;

class TychasticControlPicker {
public:
    using StrategyIndexBitFlag = Trajectory::BitFlag;
    
    OptionalCu pickControl(
        TychasticTrajectory &traj, TrajectoryPoints &trajDiscrete,
        double rho, StrategyIndexBitFlag &flag);

    OptionalCu pickControlFromSubPickerList(int startIndex,
        TychasticTrajectory &traj, TrajectoryPoints &trajDiscrete,
        double rho, StrategyIndexBitFlag &flag);

    TychasticControlPicker operator=(const TychasticControlPicker &) = delete;
    TychasticControlPicker(const ControlPicker &) = delete;
    virtual ~TychasticControlPicker();

    void setIndexSorter(indexSorter_t sortIndexes);
    int *getPreferedControlIndexes();
    int *sortPreferedControlIndexes(const double *x, double t, int strategyIndex);

    void setTycheIndex(int index);

    int getNbStrategies() const;
    int getTrajectoryIndex() const;

    std::vector<std::string> getStrategyNames() const;
    
protected:    
    TychasticControlPicker(SysDyn *sysDyn, TrajectoryParametersManager *tpm, indexSorter_t *sorterPtr);
    TychasticControlPicker *addUserPicker(TrajectoryParametersManager *params, const std::string &name);
    TychasticControlPicker *addPicker(TychasticControlPickStrategy *);
    
private:
    std::vector<TychasticControlPickStrategy *> strategies;
    SysDyn *sysDyn;
    const trajectoryParams &tp;
    controlWeight_t controlWeight;
    indexSorter_t *sortIndexes;
    int *preferedControlIndexes;
    int trajIndex;
    int tycheIndex;
};

#endif /* TYCHASTICCONTROLPICKER_H */
