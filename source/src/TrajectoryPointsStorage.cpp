#include <stdexcept>

#include "../include/TrajectoryPointsStorage.h"
#include "../include/TrajectoryPointsSimulation.h"

TrajectoryPointsStorage::TrajectoryPointsStorage(
    const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime) :
    sysDyn(sysDyn),
    duration(duration),
    dim(sysDyn->getDim())
{
    addPoint(initPoint, startingTime);
}

void TrajectoryPointsStorage::addPoint(const double *newTrajPoint, double time) {
    trajPoints.push_back(toPoint(newTrajPoint, dim));
    timeStamps.push_back(time);
}

void TrajectoryPointsStorage::writeToFile(const string &filename) {
    
    FILE *fi = fopen(filename.c_str(), "w");
    
    if (fi == nullptr) {
        spdlog::error("Impossible to open the file {}", filename);
    }
    else {

        for (std::size_t i = 0; i < trajPoints.size(); ++i) {
            for (int l1 = 0; l1 < dim; l1++)
            {
                fprintf(fi, "%15.8f ", trajPoints[i][l1]);
            }
            fprintf(fi, "%15.8f ", timeStamps[i]);
            fprintf(fi, "\n");
        }
        fclose(fi);
    }
}

const TrajectoryPointsStorage::point &TrajectoryPointsStorage::getLastPoint() const {
    return trajPoints.back();
}

double TrajectoryPointsStorage::getLastTimeStamp() const {
    return timeStamps.back();
}

const std::vector<TrajectoryPointsStorage::point> &TrajectoryPointsStorage::getPoints() const {
    return trajPoints;
}

const std::vector<double> &TrajectoryPointsStorage::getTimeStamps() const {
    return timeStamps;
}

double TrajectoryPointsStorage::getDuration() const {
    return duration;
}

TrajectoryPointsSimulation TrajectoryPointsStorage::asSimulation() {
    return TrajectoryPointsSimulation(this);
}

int TrajectoryPointsStorage::getDim() const {
    return dim;
}

const SysDyn *TrajectoryPointsStorage::getSysDyn() const {
    return sysDyn;
}



const EmptyTrajectoryPoints::point &EmptyTrajectoryPoints::getLastPoint() const {
    throw std::logic_error("getLastPoint is invalid for empty trajectory points");
}

double EmptyTrajectoryPoints::getLastTimeStamp() const {
    throw std::logic_error("getLastTimeStamp is invalid for empty trajectory points");
}

const std::vector<EmptyTrajectoryPoints::point> &EmptyTrajectoryPoints::getPoints() const {
    throw std::logic_error("getPoints is invalid for empty trajectory points");
}

const std::vector<double> &EmptyTrajectoryPoints::getTimeStamps() const {
    throw std::logic_error("getTimeStamps is invalid for empty trajectory points");
}

double EmptyTrajectoryPoints::getDuration() const {
    throw std::logic_error("getDuration is invalid for empty trajectory points");
}

TrajectoryPointsSimulation EmptyTrajectoryPoints::asSimulation() {
    throw std::logic_error("asSimulation is invalid for empty trajectory points");
}
