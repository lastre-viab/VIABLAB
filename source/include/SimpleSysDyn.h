#ifndef SIMPLESYSDYN_H_
#define SIMPLESYSDYN_H_

#include <functional>
#include <string>

#include "SysdynBase.h"
#include "ControlGrid.h"

class SimpleSysDyn : public SysdynBase
    {
public:
    SimpleSysDyn();
    SimpleSysDyn(const systemParams &SP, int stateDim, const controlParams &cp, Grid *refGrid, ControlGrid *controlGrid = nullptr);
    virtual ~SimpleSysDyn();

    void initializeMethods(const systemParams &SP) override;

    double* getLimInfC();
    double* getLimSupC();
    double* getStepC();
    unsigned long long int getDimC() const;
    unsigned long long int* getNbPointsC();
    unsigned long long int getTotalNbPointsC() const;
    double** getControlCoords() const;
    unsigned long long int** getControlIntCoords();

    virtual double calculRho_local(const double *x) const;

    int getFDDynType() const;
    std::string getRetroFileName() const;
    DynType getDynType() const;
    void setDynamicsForward();
    void setDynamicsBackward();

    void FDiscret(const double *x, const double *u, double *res, double rho) const;
    void FDiscretEuler(const double *x, const double *u, double *res, double rho) const;
    void FDiscretRK2(const double *x, const double *u, double *res, double rho) const;
    void FDiscretRK4(const double *x, const double *u, double *res, double rho) const;

protected:
    double calculL_local_num(const double *x) const;
    double calculMF_local_num(const double *x) const;
    double calculL_local_ana(const double *x) const;
    double calculMF_local_ana(const double *x) const;
    double returnL_local_ana(const double *x) const;
    double returnMF_local_ana(const double *x) const;

    ControlGrid *controls;
    void (*dynamics)(const double*, const double*, double*);
    void (*dynamics_fd)(const unsigned long long int*, const unsigned long long int*, unsigned long long int*);

    void (*localDynBounds)(const double *x, double *res);
    void (*jacobian)(const double *x, const double *u, double **jacob);

    double (*constraintsXU)(const double *x, const double *u);
    double (*constraintsXU_fd)(const unsigned long long int *x, const unsigned long long int *u);
    double (*constraintsX)(const double*);
    double (*constraintsX_fd)(const unsigned long long int*);
    double (*controlEligibilityForTraj_fd)(const unsigned long long int *x, const unsigned long long int *u, const unsigned long long int *previousU);
    double (*dynConstraintsForTraj)(const double*, double*);
    double (*target)(const double*);
    double (*target_fd)(const unsigned long long int*);
    double (*lFunc)(const double *x, const double *u);
    double (*lFunc_fd)(const unsigned long long int *x, const unsigned long long int *u);
    double (*muFunc_fd)(const unsigned long long int *x, const unsigned long long int *u);
    double (*mFunc)(const double *x, const double *u);

    std::function<double(const double*)> calculLFunc;
    std::function<double(const double*)> calculMFunc;
    std::function<void(const double*, const double*, double*, double)> discretDynamics;

    int fd_dyn_type;
    double timeHorizon;
    std::string retroFileName;
    int discretisation;
    DynType dynType;
    double lfunc_L;
    double lfunc_MF;
    double dynSignFactor;
    bool ownsControls;
    friend class SysDyn;
    };

#endif

