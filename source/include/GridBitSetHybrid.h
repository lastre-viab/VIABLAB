/*
 * GridBitSetHybrid.h
 *
 *  Created on: 17 août 2025
 *      Author: adesi
 */

#ifndef SRC_GRIDBITSETHYBRID_H_
#define SRC_GRIDBITSETHYBRID_H_

#include "GridBitSet.h"
#include "Params.h"

class GridBitSetHybrid final : public Grid_BitSet
{
public:
	GridBitSetHybrid() = default;
	GridBitSetHybrid(const gridParams &gp);
	virtual ~GridBitSetHybrid();
	unsigned long long int LocalizeHybridPoint(double * doubleCoords, unsigned long long int * intCoords);
	unsigned long long int * GetContinuousStateShifts();
	unsigned long long int getNbCellsContinuousState();
	bool isHybridPointInGrid(double *, unsigned long long int *);
	void saveHybridCoupe(const string &nomFichier) const;
	void saveValOnGridHybrid(const string &fileName) const;
	void saveValOnGridLightHybrid(const string &fileName) const;
	void loadSetHybrid(const string &fileName);

private:
	void ComputeShiftsForContinuousState();
	unsigned long long int * continuousStateShifts;
int dim_hc; // continuous state dimension dim = dim_hc + dim_hd
int dim_hd; // distrete state dimension
// convention the continuous state is encoded in axes from 0 to dim_hc -1
//Convention : discrete state is encoded in axes from dim_hc to dim-1
unsigned long long int nbCellPointsContinuousState;
};

#endif /* SRC_GRIDBITSETHYBRID_H_ */
