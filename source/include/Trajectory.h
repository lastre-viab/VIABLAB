#ifndef TRAJECTORY_H
#define TRAJECTORY_H
#include <vector>

class TrajectorySimulation;

class Trajectory {
public:
    using point = std::vector<double>;
    using BitFlag = unsigned long long int;
public:
    virtual void addPoint(int cu, const double *newTrajPoint, double time, BitFlag flag = 0) = 0;
    virtual void addPoints(int cu, int nbPoints, double endTime, double *res, BitFlag flag = 0) = 0;
    
    virtual const point &getLastPoint() const = 0;
    virtual double getLastTimeStamp() const = 0;
    
    /* Comme la trajectoire peut ne pas avoir de contrôles
       (on est encore au point initial), on ajoute une valeur "par défaut"
       si aucun contrôle n'est encore inscrit dans la liste */
    virtual int getLastControlIndex(int defaultIndex) const = 0;
    virtual double getDuration() const = 0;

    virtual const std::vector<point> &getPoints() const = 0;
    virtual const std::vector<double> &getTimeStamps() const = 0;
    virtual const std::vector<int> &getControlIndexes() const = 0;
    
    virtual TrajectorySimulation asSimulation() = 0;

    virtual ~Trajectory() = default;    
};

#endif /* TRAJECTORY_H */
