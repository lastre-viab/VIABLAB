#ifndef CONTROLPICKERBITSET_H
#define CONTROLPICKERBITSET_H

class ViabiBitSetTrajectoryHelper;
#include "../include/ControlPicker.h"

class ControlPickerBitSet final : public ControlPicker {
public:
    ControlPickerBitSet(ViabiBitSetTrajectoryHelper *viabi, SysDyn *sysDyn, TrajectoryParametersManager *params);
    ~ControlPickerBitSet() = default;
private:
    ControlPickerBitSet *addHeavyPicker(int initCu);
    ControlPickerBitSet *addFirstPicker();
    ControlPickerBitSet *addFirstOnlyPicker();
    ControlPickerBitSet *addShufflePicker();
    ControlPickerBitSet *addSortPicker();
    ControlPickerBitSet *addResetOrderPicker();
    ControlPickerBitSet *addClosestPicker();
    ControlPickerBitSet *addPreferedPicker();
    ControlPickerBitSet *addBubblePicker(const TrajectoryParametersManager *params, const Grid *grid);
    ControlPickerBitSet *addSmoothPicker(double maxAngleRadians);
    ControlPickerBitSet *addTemporalControlPicker(const TrajectoryParametersManager *params);
    ControlPickerBitSet *addBubbleBorderPicker(const TrajectoryParametersManager *params, const Grid *grid);
    
    ViabiBitSetTrajectoryHelper *viabiHelper;
};

#endif /* CONTROLPICKERBITSET_H */
