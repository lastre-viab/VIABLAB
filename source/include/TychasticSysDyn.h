#ifndef TYCHASTICSYSDYN_H_
#define TYCHASTICSYSDYN_H_

#include "SimpleSysDyn.h"

class TychasticSysDyn : public SimpleSysDyn
    {
public:
    TychasticSysDyn();
    TychasticSysDyn(const systemParams &SP, int stateDim, const controlParams &cp, Grid *refGrid, ControlGrid *controlGrid = nullptr, ControlGrid *tychesGrid = nullptr);
    ~TychasticSysDyn() override;

    void initializeMethods(const systemParams &SP) override;

    double* getLimInfTy();
    double* getLimSupTy();
    double* getStepTy();
    unsigned long long int getDimTy() const;
    unsigned long long int* getNbPointsTy();
    unsigned long long int getTotalNbPointsTy() const;
    double** getTychCoords() const;
    unsigned long long int** getTychIntCoords();

    void getTychasticImage(const double *x, const double *u, const double *v, double *imageVect, double rho) const;

protected:
    ControlGrid *tyches;
    void (*dynamics_tych)(const double*, const double*, const double*, double*);
    void (*dynamics_tych_fd)(const unsigned long long int*, const unsigned long long int*, const unsigned long long int*, unsigned long long int*);
    double (*constraintsXV_tych)(const double *x, const double *v);
    double (*lFunc_tych)(const double *x, const double *u, const double *v);
    double (*lFunc_tych_fd)(const unsigned long long int *x, const unsigned long long int *u, const unsigned long long int *v);
    double (*mFunc_tych)(const double *x, const double *u, const double *v);
    void (*jacobian_tych)(const double *x, const double *u, const double *v, double **jacob);
    std::function<void(const double*, const double*, const double*, double*, double)> discretDynamicsTych;

    double calculL_local_num_tych(const double *x) const;
    double calculL_local_ana_tych(const double *x) const;
    double calculMF_local_num_tych(const double *x) const;

    void FDiscret_tych(const double *x, const double *u, const double *v, double *res, double rho) const;
    void FDiscretEuler_tych(const double *x, const double *u, const double *v, double *res, double rho) const;
    void FDiscretRK2_tych(const double *x, const double *u, const double *v, double *res, double rho) const;
    void FDiscretRK4_tych(const double *x, const double *u, const double *v, double *res, double rho) const;
    bool ownsTyches;
    friend class SysDyn;
    };

#endif

