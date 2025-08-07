#include "../include/TrajectoryStorage.h"
#include "../include/TrajectorySimulation.h"

TrajectoryStorage TrajectoryStorage::createFlagSavingStorage(const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime) {
    return TrajectoryStorage(true, initPoint, duration, sysDyn, startingTime);
}

TrajectoryStorage TrajectoryStorage::createNonFlagSavingStorage(const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime) {
    return TrajectoryStorage(false, initPoint, duration, sysDyn, startingTime);
}

TrajectoryStorage::TrajectoryStorage(bool saveFlags, const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime) :
    points(initPoint, duration, sysDyn, startingTime),
    dimC(sysDyn->getDimC()),
    saveFlags(saveFlags)
{}

void TrajectoryStorage::addPoint(int cu, const double *newTrajPoint, double time, BitFlag flag) {
    controlIndexes.push_back(cu);
    if (saveFlags) {
        flags.push_back(flag);
    }
    points.addPoint(newTrajPoint, time);    
}

void TrajectoryStorage::addPoints(int cu, int nbPoints, double endTime, double *res, BitFlag flag) {
    const double startTime = points.getTimeStamps().back();
    const double dt = (endTime - startTime) / nbPoints;
    
    // La notation &v[0] permet de convertir un std::vector en array C
    const double *startPos = &(points.getPoints().back()[0]);
    const SysDyn *sysDyn = points.getSysDyn();
    const int dim = points.getDim();
    double *currentPos = new double[dim];
    double **controlCoords = sysDyn->getControlCoords();
    const double *control = controlCoords[cu];
    
    std::copy(startPos, startPos+dim, currentPos);    
    
    for (int i = 1; i <= nbPoints; ++i) {
        (sysDyn->*(sysDyn->discretDynamics))(currentPos, control, res, dt);
        controlIndexes.push_back(cu);
        if (saveFlags) {
            flags.push_back(flag);
        }
        
        double t = startTime + i*dt;
        
        points.addPoint(res, t);
        std::copy(res, res+dim, currentPos);
    }
    delete[] currentPos;
}

void TrajectoryStorage::writeToFile(const string &filename, const std::vector<string> &flagNames) {
    
    FILE *fi = fopen(filename.c_str(), "w");
    
    if (fi == nullptr) {
        spdlog::error("Impossible to open the file {}", filename);
    }
    else {
        const SysDyn *sysDyn = points.getSysDyn();
        double **controlCoords = sysDyn->getControlCoords();
        const std::vector<point> &trajPoints = points.getPoints();
        const std::vector<double> &timeStamps = points.getTimeStamps();
        
        // Pour qu'il y ait autant de contrôles et drapeaux que d'états
        // le dernier contrôle est dupliqué
        //
        // Puisque la trajectoire est stockée comme un ensemble s0 c1 s1 c2 s2 ... cn sn
        // Il y aura toujours un état de plus que de contrôles
        controlIndexes.push_back(controlIndexes.back());
        if (saveFlags) {
            flags.push_back(flags.back());
        }

        for (std::size_t i = 0; i < trajPoints.size(); ++i) {
            for (int l1 = 0; l1 < points.getDim(); l1++)
            {
                fprintf(fi, "%15.8f ", trajPoints[i][l1]);
            }
            fprintf(fi, "%15.8f ", timeStamps[i]);

            const double *control = controlCoords[controlIndexes[i]];
            for (int dc = 0; dc < dimC; dc++)
            {
                fprintf(fi, "%15.8f ", control[dc]);
            }
            if (saveFlags) {
                writeFlag(fi, flags[i], flagNames);
            }
            fprintf(fi, "\n");
        }
        fclose(fi);
        controlIndexes.pop_back();
        if (saveFlags) {
            flags.pop_back();
        }
    }
}

void TrajectoryStorage::writeFlag(FILE *fi, BitFlag flag, const std::vector<string> &flagNames) {
    int i = 0;
    while (flag > 1) {
        if (flag & 1) {
            fprintf(fi, "%s(%d),", flagNames[i].c_str(), i+1);
        }
        flag >>= 1;
        i++;
    }
    fprintf(fi, "%s(%d)", flagNames[i].c_str(), i+1);
}

const TrajectoryStorage::point &TrajectoryStorage::getLastPoint() const {
    return points.getLastPoint();
}

double TrajectoryStorage::getLastTimeStamp() const {
    return points.getLastTimeStamp();
}

int TrajectoryStorage::getLastControlIndex(int defaultIndex) const {
    return controlIndexes.empty()
        ? defaultIndex
        : *controlIndexes.rbegin();
}

double TrajectoryStorage::getDuration() const {
    return points.getDuration();
}

TrajectorySimulation TrajectoryStorage::asSimulation() {
    return TrajectorySimulation(this);
}

const std::vector<TrajectoryStorage::point> &TrajectoryStorage::getPoints() const {
    return points.getPoints();
}

const std::vector<double> &TrajectoryStorage::getTimeStamps() const {
    return points.getTimeStamps();
}

const std::vector<int> &TrajectoryStorage::getControlIndexes() const {
    return controlIndexes;
}



const EmptyTrajectory::point &EmptyTrajectory::getLastPoint() const {
    throw std::logic_error("getLastPoint is invalid for empty trajectory");
}

double EmptyTrajectory::getLastTimeStamp() const {
    throw std::logic_error("getLastTimeStamp is invalid for empty trajectory");
}

double EmptyTrajectory::getDuration() const {
    throw std::logic_error("getDuration is invalid for empty trajectory");
}

int EmptyTrajectory::getLastControlIndex(int defaultIndex) const {
    throw std::logic_error("getLastControlIndex is invalid for empty trajectory");
}

TrajectorySimulation EmptyTrajectory::asSimulation() {
    throw std::logic_error("asSimulation is invalid for empty trajectory");
}

const std::vector<EmptyTrajectory::point> &EmptyTrajectory::getPoints() const {
    throw std::logic_error("getPoints is invalid for empty trajectory");
}

const std::vector<double> &EmptyTrajectory::getTimeStamps() const {
    throw std::logic_error("getTimeStamps is invalid for empty trajectory");
}

const std::vector<int> &EmptyTrajectory::getControlIndexes() const {
    throw std::logic_error("getControlIndexes is invalid for empty trajectory");
}
