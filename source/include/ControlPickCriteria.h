#ifndef CONTROLPICKCRITERIA_H
#define CONTROLPICKCRITERIA_H

#include <vector>
#include "viablab_export.h"

class SysDyn;
class Trajectory;
class TrajectoryPoints;
class Grid;
struct pickedControl;

#include "ControlPicker.h"
#include "TrajectoryHelpers.h"

class VIABLAB_LIBRARY_EXPORT ControlPickCriteria
{
public:
    using StrategyIndexBitFlag = ControlPicker::StrategyIndexBitFlag;
    
    ControlPickCriteria(ControlPicker *picker,
                        Trajectory *traj, TrajectoryPoints *trajDiscrete,
                        int strategyIndex, StrategyIndexBitFlag &flag);
    ControlPickCriteria(ControlPicker *picker,
                        Trajectory *traj, TrajectoryPoints *trajDiscrete,
                        int strategyIndex, StrategyIndexBitFlag &flag, double timeStep);
    ControlPicker *getPicker();
    const ControlPicker *getPicker() const;
    void addContributionFromStrategy(int i);
    void addContributionFromStrategies(StrategyIndexBitFlag flag);
    
    SysDyn *getSysDyn();
    Trajectory *getTrajectory();
    TrajectoryPoints *getDiscreteTrajectory();
    const Grid *getGrid() const;

    const double *getCurrentPos() const;
    const double *getCurrentDiscretePos() const;
    unsigned long long int getCurrentPosGridIndex() const;
    const std::vector<std::vector<double>> &getPoints() const;
    int getDim() const;
    
    double **getControlCoords() const;
    const int *sortPreferedControlIndexes() const;
    void resetPreferedControlIndexes();
    void setIndexSorter(indexSorter_t indexSorter);
    unsigned long long int getNbControlCoords() const;
    int getDimC() const;

    double getTime() const;
    double getGridTimeStep() const;
    double getTimeStep() const;

    int getStrategyIndex() const;
    int getTrajectoryIndex() const;

    /*! Renvoie si la paire contrôle+pas de temps est viable en position courante
      La position courante étant la dernière position dans la trajectoire traj */
    bool isViableControl(pickedControl &pickedControl) const;
    bool isViableControl(pickedControl &pickedControl, double *res) const;

    /*! Renvoie si le contrôle mène à la même position dans la trajectoire discrete */
    bool isSameGridPosition(pickedControl &pickedControl) const;
    bool isSameGridPosition(pickedControl &pickedControl, double *res) const;

    /*! Modifie le pickedControl de telle sorte à ce que le pas de temps soit un
      pas de temps viable en position courante
      Un pas de temps est toujours disponible : 0. */
    void findViableTimeStep(pickedControl &controlIndex) const;
    
    ~ControlPickCriteria() = default;
private:
    ViabiTrajectoryHelper * trajectoryHelper;
    ControlPicker *picker;
    SysDyn *sysDyn;
    Trajectory *traj;
    TrajectoryPoints *trajDiscrete;
    double timeStep;
    StrategyIndexBitFlag &flag;
    int strategyIndex;
};


#endif /* CONTROLPICKCRITERIA_H */
