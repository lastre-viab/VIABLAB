#ifndef TYCHASTICTRAJECTORYSIMULATION_H
#define TYCHASTICTRAJECTORYSIMULATION_H

#include "TychasticTrajectoryStorage.h"

class TychasticTrajectorySimulation final : public TychasticTrajectory {
public:
    TychasticTrajectorySimulation(TychasticTrajectoryStorage *);
    ~TychasticTrajectorySimulation();
    
    TychasticTrajectorySimulation(TychasticTrajectorySimulation &other);
    
    void addPoint(int cu, int tycheIndex, const double *newTrajPoint, double time, BitFlag flag = 0) override;
    void addPoints(int cu, int tycheIndex, int nbPoints, double endTime, double *res, BitFlag flag = 0) override;
    
    const point &getLastPoint() const override;
    double getLastTimeStamp() const override;
    
    /* Comme la trajectoire peut ne pas avoir de contrôles
       (on est encore au point initial), on ajoute une valeur "par défaut"
       si aucun contrôle n'est encore inscrit dans la liste */
    int getLastControlIndex(int defaultIndex) const override;
    double getDuration() const override;

    const std::vector<point> &getPoints() const override;      
    const std::vector<double> &getTimeStamps() const override; 
    const std::vector<int> &getControlIndexes() const override;
    const std::vector<int> &getTycheIndexes() const override;
    
    TychasticTrajectorySimulation asSimulation() override;
    
private:
    TychasticTrajectoryStorage *storage;
    std::vector<point>::size_type simulationStartIndex;
};

#endif /* TYCHASTICTRAJECTORYSIMULATION_H */
