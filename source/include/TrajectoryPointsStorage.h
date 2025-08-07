#ifndef TRAJECTORYPOINTSSTORAGE_H
#define TRAJECTORYPOINTSSTORAGE_H

#include "SysDyn.h"
#include "TrajectoryPoints.h"

class TrajectorySimulation;
class TrajectoryPointsSimulation;
class TychasticTrajectorySimulation;

class TrajectoryPointsStorage final : public TrajectoryPoints {
public:
    TrajectoryPointsStorage(const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime = 0.0);
    void addPoint(const double *newTrajPoint, double time) override;
    
    void writeToFile(const std::string &filename);
    
    const point &getLastPoint() const override;
    double getLastTimeStamp() const override;

    const std::vector<point> &getPoints() const override;
    const std::vector<double> &getTimeStamps() const override;
    double getDuration() const override;

    TrajectoryPointsSimulation asSimulation() override;
    
    int getDim() const;
    const SysDyn *getSysDyn() const;
private:
    std::vector<point> trajPoints;
    std::vector<double> timeStamps;

    const SysDyn *sysDyn;
    double duration;
    int dim;
    
    friend TrajectoryPointsSimulation;
    friend TrajectorySimulation;
    friend TychasticTrajectorySimulation;
};

class EmptyTrajectoryPoints final : public TrajectoryPoints {
public:
    EmptyTrajectoryPoints() = default;
    void addPoint(const double *newTrajPoint, double time) override {};
    const point &getLastPoint() const override;
    double getLastTimeStamp() const override;
    
    const std::vector<point> &getPoints() const override;
    const std::vector<double> &getTimeStamps() const override;
    
    double getDuration() const;

    TrajectoryPointsSimulation asSimulation() override;
};

#endif /* TRAJECTORYPOINTSSTORAGE_H */
