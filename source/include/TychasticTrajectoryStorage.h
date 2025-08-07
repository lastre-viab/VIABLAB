#ifndef TYCHASTICTRAJECTORYSTORAGE_H
#define TYCHASTICTRAJECTORYSTORAGE_H

#include "TychasticTrajectory.h"
#include "TrajectoryPointsStorage.h"

class TychasticTrajectorySimulation;

class TychasticTrajectoryStorage final : public TychasticTrajectory {
public:    
    void writeToFile(const string &filename, const std::vector<string> &flagNames = std::vector<string>{});
    
    static TychasticTrajectoryStorage createFlagSavingStorage(const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime = 0);
    static TychasticTrajectoryStorage createNonFlagSavingStorage(const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime = 0);
    
    void addPoint(int cu, int tychIndex, const double *newTrajPoint, double time, BitFlag flag = 0) override;
    void addPoints(int cu, int tychIndex, int nbPoints, double endTime, double *res, BitFlag flag = 0) override;
    
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

    TychasticTrajectoryStorage(bool saveFlags, const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime);

    void writeFlag(FILE *fi, BitFlag flag, const std::vector<string> &flagNames);
    
    TrajectoryPointsStorage points;
    std::vector<BitFlag> flags;
    std::vector<int> controlIndexes;
    std::vector<int> tycheIndexes;
    int dimC;
    int dimTy;
    bool saveFlags;

    friend TychasticTrajectorySimulation;
};

#endif /* TYCHASTICTRAJECTORYSTORAGE_H */
