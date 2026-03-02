#ifndef HYBRIDSYSDYN_H_
#define HYBRIDSYSDYN_H_

#include "SimpleSysDyn.h"

class HybridSysDyn : public SimpleSysDyn
    {
public:
    HybridSysDyn();
    HybridSysDyn(const systemParams &SP, int continuousStateDim, int discreteStateDim, const controlParams &cp, Grid *refGrid, ControlGrid *controlGrid = nullptr, ControlGrid *hybridControls = nullptr);
    ~HybridSysDyn() override;

    void initializeMethods(const systemParams &SP) override;

    double* getLimInfHybrid();
    double* getLimSupHybrid();
    double* getStepHybrid();
    unsigned long long int getDimHybrid() const;
    unsigned long long int* getNbPointsHybrid();
    unsigned long long int getTotalNbPointsHybrid() const;
    double** getHybridCoords() const;
    unsigned long long int** getHybridIntCoords();

    double calculRho_local_hybrid(const double *xc, const unsigned long long int *xd) const;

    void FDiscret_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int *resd, double *tempReset, unsigned long long int *tempResD, double rho) const;
    void FDiscretEuler_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int *resd, double *tempReset, unsigned long long int *tempResD, double rho) const;
    void FDiscretRK2_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int *resd, double *tempReset, unsigned long long int *tempResD, double rho) const;
    void FDiscretRK4_hybrid(const double *xc, const unsigned long long int *xd, const double *uc, const unsigned long long int *ud, double *resc,
	unsigned long long int *resd, double *tempReset, unsigned long long int *tempResD, double rho) const;

protected:
    int dimS_hc;
    int dimS_hd;
    ControlGrid *hybridTransistionControls;

    void (*dynamics_hybrid_d)(const double*, const unsigned long long int*, const unsigned long long int*, unsigned long long int*);
    void (*dynamics_hybrid_c)(const double*, const unsigned long long int*, const double*, double*);
    void (*resetmap_hybrid)(const double*, const unsigned long long int*, const unsigned long long int*, double*, const unsigned long long int*);

    double (*constraintsXU_hybrid)(const double*, const unsigned long long int*, const double*, const unsigned long long int*);
    double (*constraintsX_hybrid)(const double*, const unsigned long long int*);

    void (*localDynBounds_hybrid)(const double*, const unsigned long long int*, double*);
    void (*jacobian_hybrid)(const double*, const unsigned long long int*, const double*, double**);

    std::function<double(const double*, const unsigned long long int*)> calculLHybridFunc;
    std::function<double(const double*, const unsigned long long int*)> calculMHybridFunc;
    std::function<void(const double*, const unsigned long long int*, const double*, const unsigned long long int*, double*, unsigned long long int*, double*, unsigned long long int*, double)> discretDynamicsHybrid;

    double calculL_local_num_hybrid(const double *xc, const unsigned long long int *xd) const;
    double calculL_local_ana_hybrid(const double *xc, const unsigned long long int *xd) const;
    double calculMF_local_num_hybrid(const double *xc, const unsigned long long int *xd) const;
    double calculMF_local_ana_hybrid(const double *xc, const unsigned long long int *xd) const;
    double returnL_local_ana_hybrid(const double *xc, const unsigned long long int *xd) const;
    double returnMF_local_ana_hybrid(const double *xc, const unsigned long long int *xd) const;
    bool ownsHybridControls;
    friend class SysDyn;
    };

#endif

