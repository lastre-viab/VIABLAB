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
	boost::dynamic_bitset<> analyseTrameMasque(unsigned long long int posX) const override;
  boost::dynamic_bitset<> analyseTrameMasqueWithVectorShifts(unsigned long long int posX) const override;
  bool isValidSubShift(const unsigned long long int *coords, const long long int *shift) const override;

	unsigned long long int LocalizeHybridPoint(double * doubleCoords, unsigned long long int * intCoords);
	unsigned long long int * GetContinuousStateShifts();
	unsigned long long int getNbCellsContinuousState();
	bool isHybridPointInGrid(double *, unsigned long long int *);
	void saveHybridCoupe(const string &nomFichier) const;
	void saveValOnGridHybrid(const string &fileName) const;
	void saveValOnGridLightHybrid(const string &fileName) const;
	void loadSetHybrid(const string &fileName);
	void periodizeHybridPoint(double *vect) const;

private:
	void ComputeShiftsForContinuousState();
	unsigned long long int * continuousStateShifts;
int dim_hc; // continuous state dimension dim = dim_hc + dim_hd
int dim_hd; // distrete state dimension
// convention the continuous state is encoded in axes from 0 to dim_hc -1
//Convention : discrete state is encoded in axes from dim_hc to dim-1
unsigned long long int nbCellPointsContinuousState;
  int nbPointsCubeSub_c;
  int pow3Sub_c;
  long long int *indicesDecalSub_c;
  long long int *indicesDecalCellSub_c;
};

#endif /* SRC_GRIDBITSETHYBRID_H_ */
