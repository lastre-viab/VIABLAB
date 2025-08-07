#ifndef TRAJECTORYSTORAGE_H
#define TRAJECTORYSTORAGE_H

#include "Trajectory.h"
#include "TrajectoryPointsStorage.h"

class TrajectorySimulation;

class TrajectoryStorage final : public Trajectory {
public:    
    void writeToFile(const string &filename, const std::vector<string> &flagNames = std::vector<string>{});
    
    static TrajectoryStorage createFlagSavingStorage(const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime = 0);
    static TrajectoryStorage createNonFlagSavingStorage(const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime = 0);
    
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

    TrajectoryStorage(bool saveFlags, const double *initPoint, double duration, const SysDyn *sysDyn, double startingTime);

    void writeFlag(FILE *fi, BitFlag flag, const std::vector<string> &flagNames);
    
    TrajectoryPointsStorage points;
    std::vector<BitFlag> flags;
    std::vector<int> controlIndexes;
    int dimC;
    bool saveFlags;

    friend TrajectorySimulation;
};

class EmptyTrajectory final : public Trajectory {
public:
    void addPoint(int cu, const double *newTrajPoint, double time, BitFlag flag = 0) override {};
    void addPoints(int cu, int nbPoints, double endTime, double *res, BitFlag flag = 0) override {};
    
    const point &getLastPoint() const override;
    double getLastTimeStamp() const override;
    /* Comme la trajectoire peut ne pas avoir de contrôles
       (on est encore au point initial), on ajoute une valeur "par défaut"
       si aucun contrôle n'est encore inscrit dans la liste */
    int getLastControlIndex(int defaultIndex) const override;

    const std::vector<point> &getPoints() const override;      
    const std::vector<double> &getTimeStamps() const override; 
    const std::vector<int> &getControlIndexes() const override;
    
    double getDuration() const override;

    TrajectorySimulation asSimulation() override;
};

#endif /* TRAJECTORYSTORAGE_H */
