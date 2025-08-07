#ifndef TRAJECTORYPOINTSSIMULATION_H
#define TRAJECTORYPOINTSSIMULATION_H

#include "TrajectoryPointsStorage.h"

class TrajectoryPointsSimulation final : public TrajectoryPoints {
public:
    TrajectoryPointsSimulation(TrajectoryPointsStorage *storage);

    TrajectoryPointsSimulation(TrajectoryPointsSimulation &other);

    ~TrajectoryPointsSimulation();
    
    void addPoint(const double *newTrajPoint, double time) override;
    
    void writeToFile(const std::string &filename);
    
    const point &getLastPoint() const override;
    double getLastTimeStamp() const override;

    const std::vector<point> &getPoints() const override;
    const std::vector<double> &getTimeStamps() const override;
    double getDuration() const override;

    TrajectoryPointsSimulation asSimulation() override;
private:
    TrajectoryPointsStorage *storage;
    std::vector<point>::size_type simulationStartIndex;
};

#endif /* TRAJECTORYPOINTSSIMULATION_H */
