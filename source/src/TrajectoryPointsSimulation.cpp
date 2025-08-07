#include "../include/TrajectoryPointsSimulation.h"

TrajectoryPointsSimulation::TrajectoryPointsSimulation(TrajectoryPointsStorage *storage) :
    storage(storage),
    simulationStartIndex(storage->getPoints().size())
{}

TrajectoryPointsSimulation::~TrajectoryPointsSimulation()
{
    std::vector<point> &p = storage->trajPoints;
    std::vector<double> &t = storage->timeStamps;

    p.erase(p.begin() + simulationStartIndex, p.end());
    t.erase(t.begin() + simulationStartIndex, t.end());
}

TrajectoryPointsSimulation::TrajectoryPointsSimulation(TrajectoryPointsSimulation &other)
{
    this->storage = other.storage;
    this->simulationStartIndex = storage->getPoints().size();
}


void TrajectoryPointsSimulation::addPoint(const double *newTrajPoint, double time) {
    storage->addPoint(newTrajPoint, time);
}    
    
const TrajectoryPointsSimulation::point &TrajectoryPointsSimulation::getLastPoint() const {
    return storage->getLastPoint();
}

double TrajectoryPointsSimulation::getLastTimeStamp() const {
    return storage->getLastTimeStamp();
}

const std::vector<TrajectoryPointsSimulation::point> &TrajectoryPointsSimulation::getPoints() const {
    return storage->getPoints();
}

const std::vector<double> &TrajectoryPointsSimulation::getTimeStamps() const {
    return storage->getTimeStamps();
}

double TrajectoryPointsSimulation::getDuration() const {
    return storage->getDuration();
}

TrajectoryPointsSimulation TrajectoryPointsSimulation::asSimulation() {
    return *this;
}
