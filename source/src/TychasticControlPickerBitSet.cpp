#include "../include/TychasticControlPickerBitSet.h"
#include "../include/ViabiBitSetTrajectoryHelper.h"
#include "../include/TychasticControlPickStrategies.h"

TychasticControlPickerBitSet::TychasticControlPickerBitSet(ViabiBitSetTrajectoryHelper *viabiHelper, SysDyn *sysDyn, TrajectoryParametersManager *params) :
    TychasticControlPicker(sysDyn, params, &viabiHelper->sortIndexes),
    viabiHelper(viabiHelper)
{
    const trajectoryParams *tp = params->getTrajectoryParameters();   
    
    if (tp->ARE_STRATEGIES_GUARANTEED) {
        addGuaranteedStrategies(viabiHelper, sysDyn, params);        
    }
    else {
        addNonGuaranteedStrategies(viabiHelper, sysDyn, params);
    }
}

void TychasticControlPickerBitSet::addNonGuaranteedStrategies(ViabiBitSetTrajectoryHelper *viabi, SysDyn *sysDyn, TrajectoryParametersManager *params) {
    
    const trajectoryParams *tp = params->getTrajectoryParameters();   
    const ControlPickStrategyName *strategies = tp->STRATEGIES;
    
    for (int i = 0; i < tp->NB_STRATEGIES; ++i) {
        switch (strategies[i].getPredefinedStrategyName()) {
        case FIRST:
            addFirstNonGuaranteedPicker();
            break;
        case HEAVY: {
            int initCu = viabiHelper->getClosestControlTo(tp->INIT_CONTROL);
            addHeavyPicker(initCu);        
            break;
        }
        case FIRST_ONLY:
            addFirstOnlyPicker();
            break;
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
        case CLOSEST:
            addClosestPicker();
            break;
        case PREFERED:
            addPreferedPicker();
            break;
        case BUBBLE:
            spdlog::error("Tychastic bubble strategy cannot exist, skipping");
            break;
        case SMOOTH:
            addSmoothPicker(tp->MAX_ANGLE_RADIANS);
            break;
        case TEMPORAL_CONTROL:
            if (tp->TEMPORAL_CONTROL == nullptr) {
                spdlog::warn("TEMPORAL_CONTROL strategy requested but no temporalControl function defined ? Ignoring temporalControl strategy");
            }
            else {
                addTemporalControlPicker(params);
            }
            break;
        case BUBBLE_BORDER:
            addBubbleBorderPicker(params, sysDyn->getGrid());
            break;
        case USER_STRATEGY:
            addUserPicker(params, strategies[i].getUserStrategyName());
            break;
        }
    }
}

void TychasticControlPickerBitSet::addGuaranteedStrategies(ViabiBitSetTrajectoryHelper *viabi, SysDyn *sysDyn, TrajectoryParametersManager *params) {
    
    const trajectoryParams *tp = params->getTrajectoryParameters();   
    const ControlPickStrategyName *strategies = tp->STRATEGIES;
    
    for (int i = 0; i < tp->NB_STRATEGIES; ++i) {
        switch (strategies[i].getPredefinedStrategyName()) {
        case HEAVY: {
            int initCu = viabiHelper->getClosestControlTo(tp->INIT_CONTROL);
            addHeavyPicker(initCu);
            break;
        }
        case FIRST_ONLY:
            addFirstOnlyPicker();
            break;
        case CLOSEST:
            addClosestPicker();
            break;
        case PREFERED:
            addPreferedPicker();
            break;
        case TEMPORAL_CONTROL:
            addTemporalControlPicker(params);
            break;
        case BUBBLE_BORDER:
            addBubbleBorderPicker(params, sysDyn->getGrid());
            break;
        case FIRST:
            addFirstGuaranteedPicker();
            break;
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
        case USER_STRATEGY:
            addUserPicker(params, strategies[i].getUserStrategyName());
            break;
        case SMOOTH:
            spdlog::error("Tychastic smooth guaranteed strategy not implemented");
            break;
        case BUBBLE:
            spdlog::error("Tychastic bubble strategy cannot exist, skipping");
            break;
        }
    }
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addFirstGuaranteedPicker() {     
    addPicker(new TychasticFirstGuaranteedPickStrategy());
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addFirstNonGuaranteedPicker() {
    addPicker(new TychasticFirstNonGuaranteedPickStrategyBitSet(viabiHelper));
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addHeavyPicker(int initCu) {
    addPicker(new TychasticHeavyPickStrategy(initCu));
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addFirstOnlyPicker() {
    addPicker(new TychasticFirstOnlyPickStrategy());
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addShufflePicker() {
    addPicker(new TychasticShuffleStrategy());
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addSortPicker() {
    addPicker(new TychasticSortStrategy());
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addResetOrderPicker() {
    addPicker(new TychasticResetOrderStrategy());
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addClosestPicker() {
    addPicker(new TychasticClosestPickStrategy());
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addPreferedPicker() {
    addPicker(new TychasticPreferedPickStrategy());
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addSmoothPicker(double maxAngleRadians) {
    addPicker(new TychasticSmoothPickStrategy(maxAngleRadians));
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addTemporalControlPicker(const TrajectoryParametersManager *params) {
    addPicker(new TychasticTemporalControlPickStrategy(params));
    return this;
}

TychasticControlPickerBitSet *TychasticControlPickerBitSet::addBubbleBorderPicker(const TrajectoryParametersManager *params, const Grid *grid) {
    int strategyIndex = getNbStrategies();
    addPicker(new TychasticBubbleBorderPickStrategy(strategyIndex, params, grid));
    return this;
}
