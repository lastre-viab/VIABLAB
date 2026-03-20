#ifndef SYSDYNBASE_H_
#define SYSDYNBASE_H_

#include "defs.h"
#include "Params.h"
#include "Grid.h"

class SysdynBase
    {
public:
    SysdynBase();
    explicit SysdynBase(Grid *gridRef, int dimension);
    virtual ~SysdynBase();

    SysdynBase(const SysdynBase&) = delete;
    SysdynBase& operator=(const SysdynBase&) = delete;
    SysdynBase(SysdynBase&&) = delete;
    SysdynBase& operator=(SysdynBase&&) = delete;

    virtual void initializeMethods(const systemParams &SP) = 0;

    const Grid* getGrid() const;
    int getDim() const;
    bool isTimeStepGlobal() const;

protected:
    void configureGrid(Grid *gridRef);
    void allocateStateWorkspace(int dimension);
    void releaseStateWorkspace();

    int dimS;
    bool isTychastic;
    bool isHybrid;
    double L;
    double MF;
    double timeStepFactor;
    double *image;
    double *FXmoinsH;
    double *xTemp;
    double *FXplusH;
    Grid *grid;
    bool globalTimeStep;
    ComputeMethod computeMF;
    ComputeMethod computeLC;
    double **jacob;
    };

#endif

