#include <numeric>
#include <functional> // Pour std::invoke

#include "../include/ControlPickCriteria.h"

#include "../include/SysDyn.h"
#include "../include/Trajectory.h"
#include "../include/TrajectoryPoints.h"
#include "../include/ControlPicker.h"
#include "../include/ControlPickStrategy.h"

ControlPickCriteria::ControlPickCriteria(ControlPicker *picker, Trajectory *traj, TrajectoryPoints *trajDiscrete, int strategyIndex,
	StrategyIndexBitFlag &flag) :
	picker(picker), traj(traj), trajDiscrete(trajDiscrete), flag(flag), strategyIndex(strategyIndex)
    {
    trajectoryHelper = picker->GetTrajectoryHelper();
    sysDyn = trajectoryHelper->GetDynSys();
    timeStep = sysDyn->calculRho_local(&(traj->getLastPoint()[0]));
    }

ControlPickCriteria::ControlPickCriteria(ControlPicker *picker, Trajectory *traj, TrajectoryPoints *trajDiscrete, int strategyIndex,
	StrategyIndexBitFlag &flag, double timeStep) :
	picker(picker), traj(traj), trajDiscrete(trajDiscrete), timeStep(timeStep), flag(flag), strategyIndex(strategyIndex)
    {
    trajectoryHelper = picker->GetTrajectoryHelper();
    sysDyn = trajectoryHelper->GetDynSys();
    }

ControlPicker* ControlPickCriteria::getPicker()
    {
    return picker;
    }

const ControlPicker* ControlPickCriteria::getPicker() const
    {
    return picker;
    }

void ControlPickCriteria::addContributionFromStrategy(int i)
    {
    flag |= (1 << i);
    }

void ControlPickCriteria::addContributionFromStrategies(StrategyIndexBitFlag contributorIndexes)
    {
    flag |= contributorIndexes;
    }

SysDyn* ControlPickCriteria::getSysDyn()
    {
    return sysDyn;
    }

Trajectory* ControlPickCriteria::getTrajectory()
    {
    return traj;
    }

TrajectoryPoints* ControlPickCriteria::getDiscreteTrajectory()
    {
    return trajDiscrete;
    }

const Grid* ControlPickCriteria::getGrid() const
    {
    return sysDyn->getGrid();
    }

const double* ControlPickCriteria::getCurrentPos() const
    {
    return &(traj->getLastPoint()[0]);
    }

const double* ControlPickCriteria::getCurrentDiscretePos() const
    {
    return &(trajDiscrete->getLastPoint()[0]);
    }

const std::vector<std::vector<double>>& ControlPickCriteria::getPoints() const
    {
    return traj->getPoints();
    }

unsigned long long int ControlPickCriteria::getCurrentPosGridIndex() const
    {
    return getGrid()->getNearestPointInSet(getCurrentPos());
    }

int ControlPickCriteria::getDim() const
    {
    return sysDyn->getDim();
    }

double ControlPickCriteria::getTime() const
    {
    return traj->getLastTimeStamp();
    }

double** ControlPickCriteria::getControlCoords() const
    {
    return sysDyn->getControlCoords();
    }

const int* ControlPickCriteria::sortPreferedControlIndexes() const
    {
    return picker->sortPreferedControlIndexes(getCurrentPos(), getTime(), strategyIndex);
    }

void ControlPickCriteria::resetPreferedControlIndexes()
    {
    int *preferedControlIndexes = picker->getPreferedControlIndexes();
    std::iota(preferedControlIndexes, preferedControlIndexes + getNbControlCoords(), 0);
    }

unsigned long long int ControlPickCriteria::getNbControlCoords() const
    {
    return sysDyn->getTotalNbPointsC();
    }

int ControlPickCriteria::getDimC() const
    {
    return sysDyn->getDimC();
    }

double ControlPickCriteria::getTimeStep() const
    {
    return timeStep;
    }

double ControlPickCriteria::getGridTimeStep() const
    {
    return sysDyn->calculRho_local(getCurrentDiscretePos());
    }

void ControlPickCriteria::setIndexSorter(indexSorter_t indexSorter)
    {
    picker->setIndexSorter(indexSorter);
    }

int ControlPickCriteria::getStrategyIndex() const
    {
    return strategyIndex;
    }

int ControlPickCriteria::getTrajectoryIndex() const
    {
    return picker->getTrajectoryIndex();
    }

bool ControlPickCriteria::isViableControl(pickedControl &pickedControl) const
    {
    const int dim = sysDyn->getDim();
    double imageVect[dim];
    return isViableControl(pickedControl, imageVect);
    }

bool ControlPickCriteria::isViableControl(pickedControl &pickedControl, double *imageVect) const
    {
    double rho = pickedControl.timeStep;
    int cu = pickedControl.controlIndex;
    const double *control = sysDyn->getControlCoords()[cu];
    return trajectoryHelper->isViableControl(getCurrentPos(), control, imageVect, rho);
    }

bool ControlPickCriteria::isSameGridPosition(pickedControl &pickedControl) const
    {
    const int dim = sysDyn->getDim();
    double imageVect[dim];
    if (!isViableControl(pickedControl, imageVect))
	{
	return false;
	}
    return isSameGridPosition(pickedControl, imageVect);
    }

bool ControlPickCriteria::isSameGridPosition(pickedControl &pickedControl, double *imageVect) const
    {
    return getGrid()->getNearestPointInSet(imageVect) == getCurrentPosGridIndex();
    }

void ControlPickCriteria::findViableTimeStep(pickedControl &p) const
    {
    while (!isViableControl(p))
	{
	// La décroissance du pas de temps est faite en divisant par 2 à chaque itération
	// Une autre décroissance est envisageable, mais celle-ci a l'air de fonctionner
	p.timeStep /= 2;
	}
    }
