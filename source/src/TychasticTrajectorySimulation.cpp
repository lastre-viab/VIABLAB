#include "../include/TychasticTrajectorySimulation.h"

TychasticTrajectorySimulation::TychasticTrajectorySimulation(TychasticTrajectoryStorage *storage) :
    storage(storage),
    simulationStartIndex(storage->getPoints().size())
{}

TychasticTrajectorySimulation::TychasticTrajectorySimulation(TychasticTrajectorySimulation &other) {
    this->storage = other.storage;
    this->simulationStartIndex = storage->getPoints().size();
}

TychasticTrajectorySimulation::~TychasticTrajectorySimulation()
{
    std::vector<point> &points = storage->points.trajPoints;
    std::vector<double> &timeStamps = storage->points.timeStamps;
    std::vector<int> &controlIndexes = storage->controlIndexes;
    std::vector<int> &tycheIndexes = storage->tycheIndexes;
    std::vector<BitFlag> &flags = storage->flags;    

    points.erase(points.begin() + simulationStartIndex, points.end());
    timeStamps.erase(timeStamps.begin() + simulationStartIndex, timeStamps.end());
    // -1 car il y aura toujours un contrôle, tyche, bitflag en moins que de points
    controlIndexes.erase(controlIndexes.begin() + simulationStartIndex - 1, controlIndexes.end());
    if (flags.size() > 0) {
        flags.erase(flags.begin() + simulationStartIndex - 1, flags.end());
    }
    tycheIndexes.erase(tycheIndexes.begin() + simulationStartIndex - 1, tycheIndexes.end());
}

void TychasticTrajectorySimulation::addPoint(int cu, int tycheIndex, const double *newTrajPoint, double time, BitFlag flag) {
    storage->addPoint(cu, tycheIndex, newTrajPoint, time, flag);
}    
    
const TychasticTrajectorySimulation::point &TychasticTrajectorySimulation::getLastPoint() const {
    return storage->getLastPoint();
}

void TychasticTrajectorySimulation::addPoints(int cu, int tycheIndex, int nbPoints, double endTime, double *res, BitFlag flag) {
    storage->addPoints(cu, tycheIndex, nbPoints, endTime, res, flag);
}

double TychasticTrajectorySimulation::getLastTimeStamp() const {
    return storage->getLastTimeStamp();
}

int TychasticTrajectorySimulation::getLastControlIndex(int defaultIndex) const {
    return storage->getLastControlIndex(defaultIndex);
}

double TychasticTrajectorySimulation::getDuration() const {
    return storage->getDuration();
}

const std::vector<TychasticTrajectorySimulation::point> &TychasticTrajectorySimulation::getPoints() const {
    return storage->getPoints();
}

const std::vector<double> &TychasticTrajectorySimulation::getTimeStamps() const {
    return storage->getTimeStamps();
}

const std::vector<int> &TychasticTrajectorySimulation::getControlIndexes() const {
    return storage->getControlIndexes();
}

const std::vector<int> &TychasticTrajectorySimulation::getTycheIndexes() const {
    return storage->getTycheIndexes();
}

TychasticTrajectorySimulation TychasticTrajectorySimulation::asSimulation() {
    return *this;
}
