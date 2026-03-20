#include "../include/SysdynBase.h"

SysdynBase::SysdynBase()
    : dimS(0), isTychastic(false), isHybrid(false), L(0.0), MF(0.0), timeStepFactor(1.0), image(nullptr), FXmoinsH(nullptr), xTemp(nullptr), FXplusH(nullptr),
      grid(nullptr), globalTimeStep(true), computeMF(ComputeMethod::ANALYTICAL), computeLC(ComputeMethod::ANALYTICAL), jacob(nullptr)
    {
    }

SysdynBase::SysdynBase(Grid *gridRef, int dimension)
    : SysdynBase()
    {
    configureGrid(gridRef);
    allocateStateWorkspace(dimension);
    }

SysdynBase::~SysdynBase()
    {
    releaseStateWorkspace();
    }

void SysdynBase::configureGrid(Grid *gridRef)
    {
    grid = gridRef;
    }

void SysdynBase::allocateStateWorkspace(int dimension)
    {
    releaseStateWorkspace();
    dimS = dimension;
    if (dimS <= 0)
	{
	return;
	}

    image = new double[dimS];
    FXmoinsH = new double[dimS];
    xTemp = new double[dimS];
    FXplusH = new double[dimS];

    jacob = new double*[dimS];
    for (int i = 0; i < dimS; ++i)
	{
	jacob[i] = new double[dimS];
	}
    }

void SysdynBase::releaseStateWorkspace()
    {
    if (jacob != nullptr)
	{
	for (int i = 0; i < dimS; ++i)
	    {
	    delete[] jacob[i];
	    }
	delete[] jacob;
	jacob = nullptr;
	}
    delete[] FXplusH;
    delete[] FXmoinsH;
    delete[] image;
    delete[] xTemp;
    FXplusH = nullptr;
    FXmoinsH = nullptr;
    image = nullptr;
    xTemp = nullptr;
    }

const Grid* SysdynBase::getGrid() const
    {
    return grid;
    }

int SysdynBase::getDim() const
    {
    return dimS;
    }

bool SysdynBase::isTimeStepGlobal() const
    {
    return globalTimeStep;
    }

