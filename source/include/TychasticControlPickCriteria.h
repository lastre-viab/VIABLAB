#ifndef TYCHASTICCONTROLPICKCRITERIA_H
#define TYCHASTICCONTROLPICKCRITERIA_H

#include <vector>

class SysDyn;
class TychasticTrajectory;
class TrajectoryPoints;
class Grid;
struct pickedControl;

#include "TychasticControlPicker.h"
#include "TrajectoryHelpers.h"

class TychasticControlPickCriteria
{
public:
    using StrategyIndexBitFlag = TychasticControlPicker::StrategyIndexBitFlag;
    
    TychasticControlPickCriteria(TychasticControlPicker *picker,
                                 const trajectoryParams &tp,
                                 TychasticTrajectory *traj, TrajectoryPoints *trajDiscrete,
                                 int strategyIndex, StrategyIndexBitFlag &flag, int tycheIndex);
    TychasticControlPickCriteria(TychasticControlPicker *picker,
                                 const trajectoryParams &tp,
                                 TychasticTrajectory *traj, TrajectoryPoints *trajDiscrete,
                                 int strategyIndex, StrategyIndexBitFlag &flag, int tycheIndex, double timeStep);
    
    TychasticControlPicker *getPicker();
    const TychasticControlPicker *getPicker() const;
    void addContributionFromStrategy(int i);
    void addContributionFromStrategies(StrategyIndexBitFlag flag);
    
    SysDyn *getSysDyn();
    TychasticTrajectory *getTrajectory();
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

    // Soit isViableGuaranteedControl soit isViableControlForTyche
    // Ces deux méthodes sont publiques pour les stratégies étant nécessairement garanties ou non
    bool (TychasticControlPickCriteria::*isViableControl)(pickedControl &pickedControl, int) const;
    bool isViableControlForTyche(pickedControl &pickedControl, int tycheIndex) const;
    /*! Renvoie si la paire contrôle+pas de temps est viable garantie (pour tout tyché) en position courante
      La position courante étant la dernière position dans la trajectoire traj.
    L'indice de tyché est ici ignoré, la fonction cherche à correspondre à la signature de isViableControl*/
    bool isViableGuaranteedControl(pickedControl &pickedControl, int) const;
    
    void getTychasticImage(pickedControl &pickedControl, int tycheIndex, double *imageVect) const;

    /*! Modifie le pickedControl de telle sorte à ce que le pas de temps soit un
      pas de temps viable pourle tyché donné en position courante
      Un pas de temps est toujours disponible : 0. */
    void (TychasticControlPickCriteria::*findViableTimeStep)(pickedControl &cu, int tyIndex) const;

    /*! Modifie le pickedControl de telle sorte à ce que le pas de temps soit un
      pas de temps viable pourle tyché donné en position courante
      Un pas de temps est toujours disponible : 0. */
    void findViableTimeStepForTyche(pickedControl &cu, int tyIndex) const;
    
    /*! Modifie le pickedControl de telle sorte à ce que le pas de temps soit un
      pas de temps viable garanti (pour tout tyché) en position courante
      Un pas de temps est toujours disponible : 0.
      L'indice de tyché est ignoré. On cherche à correspondre à la signature de findViableTimeStep*/
    void findViableGuaranteedTimeStep(pickedControl &cu, int) const;

    bool (TychasticControlPickCriteria::*isSameGridPosition)(pickedControl &cu, int) const;
    /*! Renvoie vrai si pour tout tyché, cu mène à la même position.
      L'indice de tyché est ignoré, on cherche à correspondre à la signature de isSameGridPosition*/
    bool isGuaranteedSameGridPosition(pickedControl &cu, int) const;
    bool isSameGridPositionForTyche(pickedControl &cu, int tyIndex) const;
    
    int getCurrentTycheIndex() const;
    double **getTychCoords();
    
    ~TychasticControlPickCriteria() = default;
private:
    ViabiTrajectoryHelper *trajectoryHelper;
    TychasticControlPicker *controlPicker;
    SysDyn *sysDyn;
    TychasticTrajectory *traj;
    TrajectoryPoints *trajDiscrete;
    double timeStep;
    StrategyIndexBitFlag &flag;
    int strategyIndex;
    int tycheIndex;
};

#endif /* TYCHASTICCONTROLPICKCRITERIA_H */
