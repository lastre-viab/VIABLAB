#ifndef VIABIBITSETTRAJECTORYHELPER_H
#define VIABIBITSETTRAJECTORYHELPER_H

#include "SysDyn.h"
#include "GridBitSet.h"
#include "Trajectory.h"
#include "TrajectoryPoints.h"
#include "Bubble.h"
#include "TrajectoryHelpers.h"
#include "ControlPickerBitSet.h"
#include "TychePicker.h"

class ViabiBitSetTrajectoryHelper {
public:
    ViabiBitSetTrajectoryHelper(Grid_BitSet *gr, SysDyn *ds, TrajectoryParametersManager *tpm);

    ~ViabiBitSetTrajectoryHelper();
    
    int computeViableTrajectory(double *initPosition, double finalTime, string fileName);
    int computeViableTrajectoryHeavy(double *initPosition, double *initControl, double finalTime, string fileName);
    
    int computeCautiousTrajectory(double *initPosition, double finalTime, string fileName);
    int computeCautiousTrajectoryHeavy(double *initPosition, double *initControl, double finalTime, string fileName);
    int computeStrategyTrajectory(double *initPosition, double finalTime, string fileName);
    int computeTychasticStrategyTrajectory(double *initPosition, double finalTime, string fileName, TychePicker &tychePicker);

    unsigned long long int findViableDiscreteSuccessor(unsigned long long int pos, double time, double dt, int &usedCu);
    unsigned long long int findViableDiscreteSuccessor_tych(unsigned long long int pos, double time, double dt, int tychIndex, int &usedCu);
    
    int findViabControl(double *currentPos, double &dt, int nbStepIter, double stepCoeff, double *resPos, int &nbViabVoisins, bool &succes);
    int findViabControl_bis(const double *currentPos, unsigned long long int optimDiscreteSuccessor, double &dt, int nbStepIter, double stepCoeff,
                            double *resPos, bool &succes);

    int findViabControl_bis_tych(const double *currentPos, unsigned long long int optimDiscreteSuccessor, double &dt, int tychIndex, int nbStepIter, double stepCoeff,
                                 double *resPos, bool &succes);
    /**
     * Trouve le contrôle "le plus proche" de u.
     * La distance choisie est ici (arbitrairement) la distance euclidienne.
     */
    int getClosestControlTo(const double *u);

    indexSorter_t sortIndexes;    
private:
    bool computeViableTrajectoryIteration(
        unsigned long long int& currentDiscreteTrajPos,
        double *xCoordsDouble,
        unsigned long long int *currentPosIntCoords, double *doubleCoordsOnDiscreteTraj, double *imageVect,
        double finalTime, double &time,
        Trajectory &traj, TrajectoryPoints &trajDiscrete,
        int realTimeStepsPerDiscreteStep);
    
    /*!
     * Renvoie si un point de grille repéré par ses coordonnées de grille
     * a dans son voisinage de Moore un point non contenu dans le noyau de viabilité
     */
    bool isOnBorder(unsigned long long int *currentPos);  

    bool findPositionAtBorder(double *xCoordsDouble,
                              Bubble &bubble,
                              double &futureTime, double &futureRho,
                              bool &onBorder, unsigned long long int &discreteBorderPos,
                              double *borderCoords, int realTimeStepsPerDiscreteStep);
    bool findPositionAtBorderHeavy(double *xCoordsDouble, const double *control,
                                   Bubble &bubble,
                                   double &futureTime, double &futureRho,
                                   bool &onBorder, unsigned long long int &discreteBorderPos, double *borderCoords);
    int applyClosestToControl(unsigned long long int currentDiscreteTrajPos,
                              int cu,
                              double rho, double normalizedTime,
                              double *xCoordsDouble,
                              double *imageVect);
    int *getPreferedControlIndexes();
    int *sortPreferedControlIndexes(const double *x, double t, int strategyIndex = 0);

private:
    Grid_BitSet *grid;
    int dim;
    int dimC;

    TrajectoryParametersManager *tpm;
    SysDyn *dynsys;

    int *preferedControlIndexes;    
    controlWeight_t controlWeight;
    int trajIndex;
};

#endif /* VIABIBITSETTRAJECTORYHELPER_H */
