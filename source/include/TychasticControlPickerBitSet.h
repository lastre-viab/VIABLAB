#ifndef TYCHASTICCONTROLPICKERBITSET_H
#define TYCHASTICCONTROLPICKERBITSET_H

class ViabiBitSetTrajectoryHelper;
#include "../include/TychasticControlPicker.h"

class TychasticControlPickerBitSet final : public TychasticControlPicker {
public:
    TychasticControlPickerBitSet(ViabiBitSetTrajectoryHelper *viabi, SysDyn *sysDyn, TrajectoryParametersManager *params);
    ~TychasticControlPickerBitSet() = default;
private:    
    void addGuaranteedStrategies(ViabiBitSetTrajectoryHelper *viabi, SysDyn *sysDyn, TrajectoryParametersManager *params);
    void addNonGuaranteedStrategies(ViabiBitSetTrajectoryHelper *viabi, SysDyn *sysDyn, TrajectoryParametersManager *params);
    
    TychasticControlPickerBitSet *addFirstGuaranteedPicker();

    TychasticControlPickerBitSet *addFirstNonGuaranteedPicker();
    TychasticControlPickerBitSet *addHeavyPicker(int initCu);
    TychasticControlPickerBitSet *addFirstOnlyPicker();
    TychasticControlPickerBitSet *addShufflePicker();
    TychasticControlPickerBitSet *addSortPicker();
    TychasticControlPickerBitSet *addResetOrderPicker();
    TychasticControlPickerBitSet *addClosestPicker();
    TychasticControlPickerBitSet *addPreferedPicker();
    TychasticControlPickerBitSet *addSmoothPicker(double maxAngleRadians);
    TychasticControlPickerBitSet *addTemporalControlPicker(const TrajectoryParametersManager *params);
    TychasticControlPickerBitSet *addBubbleBorderPicker(const TrajectoryParametersManager *params, const Grid *grid);

    ViabiBitSetTrajectoryHelper *viabiHelper;
};

#endif /* CONTROLPICKERBITSET_H */
