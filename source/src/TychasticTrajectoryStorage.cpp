#include "../include/TychasticTrajectoryStorage.h"
#include "../include/TychasticTrajectorySimulation.h"

TychasticTrajectoryStorage TychasticTrajectoryStorage::createFlagSavingStorage(const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime) {
    return TychasticTrajectoryStorage(true, initPoint, duration, sysDyn, startingTime);
}

TychasticTrajectoryStorage TychasticTrajectoryStorage::createNonFlagSavingStorage(const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime) {
    return TychasticTrajectoryStorage(false, initPoint, duration, sysDyn, startingTime);
}

TychasticTrajectoryStorage::TychasticTrajectoryStorage(bool saveFlags, const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime) :
    points(initPoint, duration, sysDyn, startingTime),
    dimC(sysDyn->getDimC()),
    dimTy(sysDyn->getDimTy()),
    saveFlags(saveFlags)
{}

void TychasticTrajectoryStorage::addPoint(int cu, int tycheIndex, const double *newTrajPoint, double time, BitFlag flag) {
    controlIndexes.push_back(cu);
    tycheIndexes.push_back(tycheIndex);
    if (saveFlags) {
        flags.push_back(flag);
    }
    points.addPoint(newTrajPoint, time);    
}

void TychasticTrajectoryStorage::addPoints(int cu, int tycheIndex, int nbPoints, double endTime, double *res, BitFlag flag) {
    const double startTime = points.getTimeStamps().back();
    const double dt = (endTime - startTime) / nbPoints;
    
    // La notation &v[0] permet de convertir un std::vector en array C
    const double *startPos = &(points.getPoints().back()[0]);
    const SysDyn *sysDyn = points.getSysDyn();
    const int dim = points.getDim();
    double *currentPos = new double[dim];
    double **controlCoords = sysDyn->getControlCoords();
    double **tycheCoords = sysDyn->getTychCoords();
    const double *control = controlCoords[cu];
    const double *tyche = tycheCoords[tycheIndex];
    
    std::copy(startPos, startPos+dim, currentPos);    
    
    for (int i = 1; i <= nbPoints; ++i) {
        (sysDyn->*(sysDyn->discretDynamics_tych))(currentPos, control, tyche, res, dt);
        controlIndexes.push_back(cu);
        tycheIndexes.push_back(tycheIndex);
        if (saveFlags) {
            flags.push_back(flag);
        }
        
        double t = startTime + i*dt;
        
        points.addPoint(res, t);
        std::copy(res, res+dim, currentPos);
    }
    delete[] currentPos;
}

void TychasticTrajectoryStorage::writeToFile(const string &filename, const std::vector<string> &flagNames) {
    
    FILE *fi = fopen(filename.c_str(), "w");
    
    if (fi == nullptr) {
        spdlog::error("Impossible to open the file {}", filename);
    }
    else {
        const SysDyn *sysDyn = points.getSysDyn();
        double **controlCoords = sysDyn->getControlCoords();
        double **tycheCoords = sysDyn->getTychCoords();
        const std::vector<point> &trajPoints = points.getPoints();
        const std::vector<double> &timeStamps = points.getTimeStamps();
        
        // Pour qu'il y ait autant de contrôles, tychés et drapeaux que d'états
        // le dernier contrôle est dupliqué
        //
        // Puisque la trajectoire est stockée comme un ensemble s0 c1 s1 c2 s2 ... cn sn
        // Il y aura toujours un état de plus que de contrôles
        controlIndexes.push_back(controlIndexes.back());
        tycheIndexes.push_back(tycheIndexes.back());
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

            const double *tyche = tycheCoords[tycheIndexes[i]];
            for (int dty = 0; dty < dimTy; dty++)
            {
                fprintf(fi, "%15.8f ", tyche[dty]);
            }
            if (saveFlags) {
                writeFlag(fi, flags[i], flagNames);
            }
            fprintf(fi, "\n");
        }
        fclose(fi);
        controlIndexes.pop_back();
        tycheIndexes.pop_back();
        if (saveFlags) {
            flags.pop_back();
        }
    }
}

void TychasticTrajectoryStorage::writeFlag(FILE *fi, BitFlag flag, const std::vector<string> &flagNames) {
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

const TychasticTrajectoryStorage::point &TychasticTrajectoryStorage::getLastPoint() const {
    return points.getLastPoint();
}

double TychasticTrajectoryStorage::getLastTimeStamp() const {
    return points.getLastTimeStamp();
}

int TychasticTrajectoryStorage::getLastControlIndex(int defaultIndex) const {
    return controlIndexes.empty()
        ? defaultIndex
        : *controlIndexes.rbegin();
}

double TychasticTrajectoryStorage::getDuration() const {
    return points.getDuration();
}

TychasticTrajectorySimulation TychasticTrajectoryStorage::asSimulation() {
    return TychasticTrajectorySimulation(this);
}

const std::vector<TychasticTrajectoryStorage::point> &TychasticTrajectoryStorage::getPoints() const {
    return points.getPoints();
}

const std::vector<double> &TychasticTrajectoryStorage::getTimeStamps() const {
    return points.getTimeStamps();
}

const std::vector<int> &TychasticTrajectoryStorage::getControlIndexes() const {
    return controlIndexes;
}

const std::vector<int> &TychasticTrajectoryStorage::getTycheIndexes() const {
    return tycheIndexes;
}
