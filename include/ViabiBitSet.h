/*
 * ViabiBitSet.h
 *  *
 *    VIABLAB : a numerical library for Mathematical Viability Computations
 *    Copyright (C) <2020>  <Anna DESILLES, LASTRE>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *   
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Created on: 17 sept. 2013
 *      Author: ANYA
 */

#ifndef VIABIBITSET_H_
#define VIABIBITSET_H_

#include "Viabi.h"
#include "GridBitSet.h"
#include "ParametersManager.h"


class ViabiBitSet: public Viabi {
public:
	ViabiBitSet(ParametersManager *pm);
	virtual ~ViabiBitSet();
	virtual void printViabiInfo();

	virtual void ViabilityKernel( bool sortieOK,int nbArret);
	virtual void  CaptureBasin();
	virtual void GarantedViabilityKernel( bool sortieOK,int nbArret);

	virtual	void initialiseTarget();
	virtual	void initialiseConstraints();
	virtual void computeTrajectories();

	virtual void loadViableSets();
	virtual void saveViableSets();

	void setK0();
	void setK0_fd();
	void testK0();

	void  saveViabRetro(string fileName);
	void  saveViabGarantiRetro(string fileName);





	bool  findViabImagePoint(double *currentPos , bool print);
	bool findViabImagePoint_noControl(double *xCoordsDouble, bool print);
	void  initialiseTargetPointList();
	int  computeViableTrajectorySetVal(double *initPosition, double finalTime, string fileName);
	int computeViableTrajectoryHeavy(double *initPosition, double *initControl, double finalTime, string fileName);

	int computeViableTrajectory(double *initPosition, double finalTime, string fileName);
	int findViabControl(double *currentPos,
			double &dt,
			int nbStepIter,
			double stepCoeff,
			double *resPos,
			int & nbViabVoisins,
			bool &succes );
private:
	Grid_BitSet* grid;
	void InitViabiBitSet(algoViabiParams avbp);
	void ViabilityKernelSimple( bool sortieOK,int nbArret);
	void  computeDiscreteImageOfPoint(double *doublePointCoords, unsigned long long int * intPointCoords  );
	void  computeDiscreteImageOfPoint_noControl(double *doublePointCoords, unsigned long long int * intPointCoords  );

	void  noyauViabi_FD( bool sortieOK,int nbArret);
	void  noyauViabi( bool sortieOK,int nbArret);
	void  noyauViabi_sansControle( bool sortieOK,int nbArret);
	void noyauViabi_omp( bool sortieOK,int nbArret);
	void noyauViabi_sansControle_omp( bool sortieOK,int nbArret);
	void noyauViabiGaranti_FD( bool sortieOK,int nbArret);
	void CaptureBasin_ContinuousDynamics();
	void CaptureBasin_DiscreteDynamics();

		void computeViableTrajectories();
	/*!
	 *  \brief  Copie pour raisons de rapidt�  de la valeur de dimension d'�tat
	 */
	int dim;
	/*!
	 *  \brief  Copie pour raisons de rapidt�  de la valeur de dimension de contr�le
	 */
	int dimC;
	/*!
	 *  \brief Pointeur sur la base de donn�es servant � enregister la r�troaction optimale
	 */
	string filePrefix;



	/*!
	 *  \brief  Une structure servant � stocker les donn�es  de l'image discr�te
	 *  d'un point au cours de parcours  de calcul  de \f$ \Phi(C_{n}\setminus C_{n-1}) \f$ .
	 *  \see computeDiscreteImageOfPoint
	 *  \see computeCurrentImageGlobalRho
	 *
	 *  Attention! Variable globale dans les m�thodes de la classe!
	 */
	discretImageSet_simple pointDI;

	double   *doubleVect, *doubleVect1;
	unsigned long long int * imageCells;

	/*!
	 *  \brief La liste de structures repr�sentant  des mailles d'une image en construction
	 * Attention! Les m�thodes suivantes modifient  cette variable comme �tant globale
	 * \see addDataToCurrentImage
	 * \see computeCurrentImageGlobalRho
	 */
	imageCellsIntList currentImageList;
	/*!
	 * \brief La liste de structures repr�sentant  des points d'une image en construction
	 * Attention! Les m�thodes suivantes modifient  cette variable comme �tant globale
	 * \see createPointsListGlobalRho
	 * \see addDataToPointsList
	 */
	imagePointsIntList currentImagePointsList;



	void  addConvexCombinations(unsigned long long int posX, unsigned long long int numCell, unsigned long long int * tempImageCell, double rho ,list<unsigned long long int>::iterator *itStart);
	void  computeConvexifiedImage( int iter);
	void  addDataToCurrentImage(list<unsigned long long int >::iterator *startIt, unsigned long long int newCell,list<unsigned long long int>::iterator *resIt );
	int  addNewPoints();
	void  computeCurrentImage( int iter);
	void  createPointsList();
	void addDataToPointsList(list<unsigned long long int >::iterator *startIt, unsigned long long int newPoint,list<unsigned long long int>::iterator *resIt );

};



#endif /* VIABIBITSET_H_ */
