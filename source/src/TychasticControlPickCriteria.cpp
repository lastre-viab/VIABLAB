#include <numeric>
#include <functional> // Pour std::invoke

#include "../include/TychasticControlPickCriteria.h"

#include "../include/SysDyn.h"
#include "../include/TychasticTrajectory.h"
#include "../include/TrajectoryPoints.h"
#include "../include/TychasticControlPicker.h"
#include "../include/TychasticControlPickStrategy.h"

TychasticControlPickCriteria::TychasticControlPickCriteria(TychasticControlPicker *picker, const trajectoryParams &tp, TychasticTrajectory *traj, TrajectoryPoints *trajDiscrete, int strategyIndex, StrategyIndexBitFlag &flag, int tycheIndex):
controlPicker(picker),
traj(traj),
trajDiscrete(trajDiscrete),
flag(flag),
strategyIndex(strategyIndex),
tycheIndex(tycheIndex)
{
	trajectoryHelper = controlPicker->GetTrajectoryHelper();
	sysDyn = trajectoryHelper->GetDynSys();
	timeStep = sysDyn->calculRho_local(&(traj->getLastPoint()[0]));
	if (tp.ARE_STRATEGIES_GUARANTEED) {
		TychasticControlPickCriteria::isViableControl = &TychasticControlPickCriteria::isViableGuaranteedControl;
		TychasticControlPickCriteria::findViableTimeStep = &TychasticControlPickCriteria::findViableGuaranteedTimeStep;
		TychasticControlPickCriteria::isSameGridPosition = &TychasticControlPickCriteria::isGuaranteedSameGridPosition;
	}
	else {
		TychasticControlPickCriteria::isViableControl = &TychasticControlPickCriteria::isViableControlForTyche;
		TychasticControlPickCriteria::findViableTimeStep = &TychasticControlPickCriteria::findViableTimeStepForTyche;
		TychasticControlPickCriteria::isSameGridPosition = &TychasticControlPickCriteria::isSameGridPositionForTyche;
	}
}

TychasticControlPickCriteria::TychasticControlPickCriteria(TychasticControlPicker *picker, const trajectoryParams &tp, TychasticTrajectory *traj, TrajectoryPoints *trajDiscrete, int strategyIndex, StrategyIndexBitFlag &flag, int tycheIndex, double timeStep) :
    				controlPicker(picker),
					traj(traj),
					trajDiscrete(trajDiscrete),
					timeStep(timeStep),
					flag(flag),
					strategyIndex(strategyIndex),
					tycheIndex(tycheIndex)
{
	trajectoryHelper = controlPicker->GetTrajectoryHelper();
	sysDyn = trajectoryHelper->GetDynSys();

	if (tp.ARE_STRATEGIES_GUARANTEED) {
		TychasticControlPickCriteria::isViableControl = &TychasticControlPickCriteria::isViableGuaranteedControl;
		TychasticControlPickCriteria::findViableTimeStep = &TychasticControlPickCriteria::findViableGuaranteedTimeStep;
		TychasticControlPickCriteria::isSameGridPosition = &TychasticControlPickCriteria::isGuaranteedSameGridPosition;
	}
	else {
		TychasticControlPickCriteria::isViableControl = &TychasticControlPickCriteria::isViableControlForTyche;
		TychasticControlPickCriteria::findViableTimeStep = &TychasticControlPickCriteria::findViableTimeStepForTyche;
		TychasticControlPickCriteria::isSameGridPosition = &TychasticControlPickCriteria::isSameGridPositionForTyche;
	}
}

TychasticControlPicker *TychasticControlPickCriteria::getPicker() {
	return controlPicker;
}

const TychasticControlPicker *TychasticControlPickCriteria::getPicker() const {
	return controlPicker;
}

void TychasticControlPickCriteria::addContributionFromStrategy(int i) {
	flag |= (1<<i);
}

void TychasticControlPickCriteria::addContributionFromStrategies(StrategyIndexBitFlag contributorIndexes) {
	flag |= contributorIndexes;
}

SysDyn *TychasticControlPickCriteria::getSysDyn() {
	return sysDyn;
}

TychasticTrajectory *TychasticControlPickCriteria::getTrajectory() {
	return traj;
}

TrajectoryPoints *TychasticControlPickCriteria::getDiscreteTrajectory() {
	return trajDiscrete;
}

const Grid *TychasticControlPickCriteria::getGrid() const {
	return sysDyn->getGrid();
}

const double *TychasticControlPickCriteria::getCurrentPos() const {
	return &(traj->getLastPoint()[0]);
}

const double *TychasticControlPickCriteria::getCurrentDiscretePos() const {
	return &(trajDiscrete->getLastPoint()[0]);
}

const std::vector<std::vector<double>> &TychasticControlPickCriteria::getPoints() const {
	return traj->getPoints();
}

unsigned long long int TychasticControlPickCriteria::getCurrentPosGridIndex() const {
	return getGrid()->getNearestPointInSet(getCurrentPos());
}

int TychasticControlPickCriteria::getDim() const {
	return sysDyn->getDim();
}

double TychasticControlPickCriteria::getTime() const {
	return traj->getLastTimeStamp();
}

double **TychasticControlPickCriteria::getControlCoords() const {
	return sysDyn->getControlCoords();
}

const int *TychasticControlPickCriteria::sortPreferedControlIndexes() const {
	return controlPicker->sortPreferedControlIndexes(
			getCurrentPos(),
			getTime(),
			strategyIndex);
}

void TychasticControlPickCriteria::resetPreferedControlIndexes() {
	int *preferedControlIndexes = controlPicker->getPreferedControlIndexes();
	std::iota(preferedControlIndexes, preferedControlIndexes+getNbControlCoords(), 0);
}

unsigned long long int TychasticControlPickCriteria::getNbControlCoords() const {
	return sysDyn->getTotalNbPointsC();
}

int TychasticControlPickCriteria::getDimC() const {
	return sysDyn->getDimC();
}

double TychasticControlPickCriteria::getTimeStep() const {
	return timeStep;
}

double TychasticControlPickCriteria::getGridTimeStep() const {
	return sysDyn->calculRho_local(getCurrentDiscretePos());
}

void TychasticControlPickCriteria::setIndexSorter(indexSorter_t indexSorter) {
	controlPicker->setIndexSorter(indexSorter);
}

int TychasticControlPickCriteria::getStrategyIndex() const {
	return strategyIndex;
}

int TychasticControlPickCriteria::getTrajectoryIndex() const {
	return controlPicker->getTrajectoryIndex();
}

bool TychasticControlPickCriteria::isViableGuaranteedControl(pickedControl &pickedControl, int) const {
	double rho = pickedControl.timeStep;
	int cu = pickedControl.controlIndex;
	const double *control = sysDyn->getControlCoords()[cu];
	return trajectoryHelper->isViableGuaranteedControl(getCurrentPos(), control, rho);
}

bool TychasticControlPickCriteria::isViableControlForTyche(pickedControl &pickedControl, int tycheIndex) const {
	const int dim = sysDyn->getDim();
	double *imageVect = new double[dim];

	double rho = pickedControl.timeStep;
	int cu = pickedControl.controlIndex;
	const double *control = sysDyn->getControlCoords()[cu];
	const double *tyche = sysDyn->getTychCoords()[tycheIndex];

	bool res = (trajectoryHelper->isViableControl_tych(getCurrentPos(), control, tyche, imageVect, rho));

	delete [] imageVect;

	return res;
}

void TychasticControlPickCriteria::getTychasticImage(pickedControl &pickedControl, int tycheIndex, double *imageVect) const {

	double rho = pickedControl.timeStep;
	int cu = pickedControl.controlIndex;
	const double *control = sysDyn->getControlCoords()[cu];
	const double *tyche = sysDyn->getControlCoords()[tycheIndex];

	return sysDyn->getTychasticImage(getCurrentPos(), control, tyche, imageVect, rho);
}

void TychasticControlPickCriteria::findViableTimeStepForTyche(pickedControl &p, int tyIndex) const {
	while (!isViableControlForTyche(p, tyIndex)) {
		// La décroissance du pas de temps est faite en divisant par 2 à chaque itération
		// Une autre décroissance est possible, mais celle-ci a l'air de fonctionner
		p.timeStep /= 2;
	}
}

void TychasticControlPickCriteria::findViableGuaranteedTimeStep(pickedControl &p, int) const {
	while (!isViableGuaranteedControl(p, 0)) {
		// La décroissance du pas de temps est faite en divisant par 2 à chaque itération
		// Une autre décroissance est possible, mais celle-ci a l'air de fonctionner
		p.timeStep /= 2;
	}
}

bool TychasticControlPickCriteria::isGuaranteedSameGridPosition(pickedControl &pickedControl, int) const {
	bool allSame = true;
	unsigned long long int nbPoints = sysDyn->getTotalNbPointsTy();
	unsigned long long int i = 0;

	const int dim = sysDyn->getDim();
	double *imageVect = new double[dim];

	const double *xCoordsDouble = getCurrentPos();
	double **tycheCoords = sysDyn->getTychCoords();

	while (allSame && i < nbPoints) {
		const double *tyche = tycheCoords[i];
		if (sysDyn->constraintsXV_tych(xCoordsDouble, tyche) < PLUS_INF) {
			getTychasticImage(pickedControl, i, imageVect);
			allSame = (getGrid()->getNearestPointInSet(imageVect) == getCurrentPosGridIndex());
		}
		++i;
	}

	delete [] imageVect;

	return allSame;
}

bool TychasticControlPickCriteria::isSameGridPositionForTyche(pickedControl &cu, int tyIndex) const {
	const int dim = sysDyn->getDim();
	double *imageVect = new double[dim];
	getTychasticImage(cu, tyIndex, imageVect);

	bool res = (getGrid()->getNearestPointInSet(imageVect) == getCurrentPosGridIndex());

	delete [] imageVect;

	return res;
}

double **TychasticControlPickCriteria::getTychCoords() {
	return sysDyn->getTychCoords();
}

int TychasticControlPickCriteria::getCurrentTycheIndex() const {
	return tycheIndex;
}
