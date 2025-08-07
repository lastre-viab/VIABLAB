#include "../include/TrajectorySimulation.h"

TrajectorySimulation::TrajectorySimulation(TrajectoryStorage *storage) :
    storage(storage),
    simulationStartIndex(storage->getPoints().size())
{}

TrajectorySimulation::TrajectorySimulation(TrajectorySimulation &other) {
    this->storage = other.storage;
    this->simulationStartIndex = storage->getPoints().size();
}

TrajectorySimulation::~TrajectorySimulation()
{
    std::vector<point> &p = storage->points.trajPoints;
    std::vector<double> &t = storage->points.timeStamps;
    std::vector<int> &c = storage->controlIndexes;
    std::vector<BitFlag> &f = storage->flags;

    p.erase(p.begin() + simulationStartIndex, p.end());
    t.erase(t.begin() + simulationStartIndex, t.end());
    // -1 car il y aura toujours un contrôle en moins que de points
    c.erase(c.begin() + simulationStartIndex - 1, c.end());
    if (f.size() > 0) {
        f.erase(f.begin() + simulationStartIndex - 1, f.end());
    }
}

void TrajectorySimulation::addPoint(int cu, const double *newTrajPoint, double time, BitFlag flag) {
    storage->addPoint(cu, newTrajPoint, time, flag);
}    
    
const TrajectorySimulation::point &TrajectorySimulation::getLastPoint() const {
    return storage->getLastPoint();
}

void TrajectorySimulation::addPoints(int cu, int nbPoints, double endTime, double *res, BitFlag flag) {
    storage->addPoints(cu, nbPoints, endTime, res, flag);
}

double TrajectorySimulation::getLastTimeStamp() const {
    return storage->getLastTimeStamp();
}

int TrajectorySimulation::getLastControlIndex(int defaultIndex) const {
    return storage->getLastControlIndex(defaultIndex);
}

double TrajectorySimulation::getDuration() const {
    return storage->getDuration();
}

const std::vector<TrajectorySimulation::point> &TrajectorySimulation::getPoints() const {
    return storage->getPoints();
}

const std::vector<double> &TrajectorySimulation::getTimeStamps() const {
    return storage->getTimeStamps();
}

const std::vector<int> &TrajectorySimulation::getControlIndexes() const {
    return storage->getControlIndexes();
}

TrajectorySimulation TrajectorySimulation::asSimulation() {
    return *this;
}
