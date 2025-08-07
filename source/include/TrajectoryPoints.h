#ifndef TRAJECTORYPOINTS_H
#define TRAJECTORYPOINTS_H

#include <vector>

class TrajectoryPointsSimulation;

class TrajectoryPoints {
public:
    typedef std::vector<double> point;
    static point toPoint(const double *tab, int len) {
        return point(tab, tab+len);
    }
public:
    virtual void addPoint(const double *newTrajPoint, double time) = 0;
    virtual const point &getLastPoint() const = 0;
    virtual double getLastTimeStamp() const = 0;
    virtual const std::vector<point> &getPoints() const = 0;
    virtual const std::vector<double> &getTimeStamps() const = 0;
    virtual double getDuration() const = 0;

    virtual TrajectoryPointsSimulation asSimulation() = 0;
    
    virtual ~TrajectoryPoints() = default;
};

#endif /* TRAJECTORYPOINTS_H */
