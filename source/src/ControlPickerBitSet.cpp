#include "../include/ControlPickerBitSet.h"
#include "../include/ViabiBitSetTrajectoryHelper.h"
#include "../include/ControlPickStrategies.h"

ControlPickerBitSet::ControlPickerBitSet(ViabiTrajectoryHelper *viabiHelper, TrajectoryParametersManager *params) :
    ControlPicker(viabiHelper, params, &viabiHelper->sortIndexes)
{
    const trajectoryParams *tp = params->getTrajectoryParameters();   
    const ControlPickStrategyName *strategies = tp->STRATEGIES;
    
    for (int i = 0; i < tp->NB_STRATEGIES; ++i) {
        switch (strategies[i].getPredefinedStrategyName()) {
        case FIRST:
            addFirstPicker();
            break;
        case HEAVY: {
            int initCu = viabiHelper->getClosestControlTo(tp->INIT_CONTROL);
            addHeavyPicker(initCu);
            break;
        }
        case SHUFFLE:
            addShufflePicker();
            break;
        case SORT:
            if (tp->CONTROL_WEIGHT == nullptr) {
                spdlog::warn("SORT strategy requested but no controlWeight function defined ? Ignoring sort strategy");
            }
            else {
                addSortPicker();
            }
            break;
        case RESET_ORDER:
            addResetOrderPicker();
            break;
        case FIRST_ONLY:
            addFirstOnlyPicker();
            break;
        case CLOSEST:
            addClosestPicker();
            break;
        case PREFERED:
            addPreferedPicker();
            break;
        case BUBBLE:
            addBubblePicker(params, sysDyn->getGrid());
            break;
        case BUBBLE_BORDER:
            addBubbleBorderPicker(params, sysDyn->getGrid());
            break;
        case SMOOTH:
            addSmoothPicker(tp->MAX_ANGLE_RADIANS);
            break;
        case TEMPORAL_CONTROL: {
            if (tp->TEMPORAL_CONTROL == nullptr) {
                spdlog::warn("TEMPORAL_CONTROL strategy requested but no temporalControl function defined ? Ignoring temporalControl strategy");
            }
            else {
                addTemporalControlPicker(params);
            }
            break;
        }
        case USER_STRATEGY:
            addUserPicker(params, strategies[i].getUserStrategyName());
            break;            
        }
    }
}

ControlPickerBitSet *ControlPickerBitSet::addHeavyPicker(int initCu) {       
    addPicker(new HeavyPickStrategy(initCu));
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addFirstPicker() {     
    addPicker(new FirstPickStrategyBitSet(dynamic_cast<ViabiBitSetTrajectoryHelper*>(trajectoryHelper)));
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addFirstOnlyPicker() {     
    addPicker(new FirstOnlyPickStrategy());
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addShufflePicker() {      
    addPicker(new ShuffleStrategy());
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addSortPicker() {
    addPicker(new SortStrategy());
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addResetOrderPicker() {
    addPicker(new ResetOrderStrategy());
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addClosestPicker() {
    addPicker(new ClosestPickStrategy());
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addPreferedPicker() {
    addPicker(new PreferedPickStrategy());
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addBubblePicker(const TrajectoryParametersManager *params, const Grid *grid) {
    int strategyIndex = getNbStrategies();
    addPicker(new BubblePickStrategy(strategyIndex, params, grid));
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addBubbleBorderPicker(const TrajectoryParametersManager *params, const Grid *grid) {
    int strategyIndex = getNbStrategies();
    addPicker(new BubbleBorderPickStrategy(strategyIndex, params, grid));
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addSmoothPicker(double maxAngleRadians) {
    addPicker(new SmoothPickStrategy(maxAngleRadians));
    return this;
}

ControlPickerBitSet *ControlPickerBitSet::addTemporalControlPicker(const TrajectoryParametersManager *params) {
    addPicker(new TemporalControlPickStrategy(params));
    return this;
}
