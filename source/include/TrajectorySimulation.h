#ifndef TRAJECTORYSIMULATION_H
#define TRAJECTORYSIMULATION_H

#include "TrajectoryStorage.h"

class TrajectorySimulation final : public Trajectory {
public:
    TrajectorySimulation(TrajectoryStorage *);
    ~TrajectorySimulation();
    
    TrajectorySimulation(TrajectorySimulation &other);
    
    void addPoint(int cu, const double *newTrajPoint, double time, BitFlag flag = 0) override;
    void addPoints(int cu, int nbPoints, double endTime, double *res, BitFlag flag = 0) override;
    
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
    
    TrajectorySimulation asSimulation() override;
    
private:
    TrajectoryStorage *storage;
    std::vector<point>::size_type simulationStartIndex;
};

#endif /* TRAJECTORYSIMULATION_H */
