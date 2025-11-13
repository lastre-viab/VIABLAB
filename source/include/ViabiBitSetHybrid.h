/*
 * ViabiBitSetHybrid.h
 *
 *  Created on: 21 août 2025
 *      Author: adesi
 */

#ifndef SRC_VIABIBITSETHYBRID_H_
#define SRC_VIABIBITSETHYBRID_H_

#include "Viabi.h"
#include "GridBitSetHybrid.h"

class ViabiBitSetHybrid final : public Viabi
    {
public:
    ViabiBitSetHybrid();
    virtual ~ViabiBitSetHybrid();
    ViabiBitSetHybrid(ParametersManager *pm);


    virtual void printViabiInfo();

    virtual void ViabilityKernel(int nbArret);
    virtual void CaptureBasin();
    virtual void GarantedViabilityKernel(int nbArret);

    virtual void initialiseTarget();
    virtual void initialiseConstraints();
    virtual void computeTrajectories();

    virtual void loadViableSets();
    virtual void saveViableSets();

    virtual SysDyn* GetSysDynForViabProblem();

    void saveViableSets(const string &baseFilenameSuffix);
    void saveIntermediateViableSets(int refine);
    void setK0();
    void setK0_fd();
    void testK0();

    void saveViabRetro(string fileName);
    void saveViabGarantiRetro(string fileName);

    bool findViabImagePoint(double *currentPos, unsigned long long int * currentIntCoords, bool print);
    bool findGuarantedViabImagePoint(double *xCoordsDouble, bool print);
    bool findViabImagePoint_noControl(double *xCoordsDouble, bool print);
    void initialiseTargetPointList();

private:
    GridBitSetHybrid *grid;

    string getSetName(SetType type);
    void InitViabiBitSetHybrid(const algoViabiParams &avbp);
    void ViabilityKernelSimple(int nbArret);
    void GuarantedViabilityKernelSimple(int nbArret);
    void computeDiscreteImageOfPoint(double *doublePointCoords, unsigned long long int *intPointCoords);
    void computeDiscreteImageOfPoint_noControl(double *doublePointCoords, unsigned long long int *intPointCoords);

    void noyauViabi_FD(int nbArret);
    void noyauViabi(int nbArret);
    void noyauViabiGuaranti(int nbArret);
    void noyauViabi_sansControle(int nbArret);
    void noyauViabi_omp(int nbArret);
    void noyauViabi_sansControle_omp(int nbArret);
    void noyauViabiGaranti_FD(int nbArret);
    void CaptureBasin_ContinuousDynamics();
    void CaptureBasin_DiscreteDynamics();
    void computeViableTrajectories();
    void computeTychasticViableTrajectories();
    bool CheckViability(double *xCoordsDouble, unsigned long long int * xCoordsInt, double * uc, unsigned long long int * ud,
    	double *doubleImage,
    	unsigned long long int *intImage,
    	unsigned long long int *gridCoords,
	double * tempResetc, unsigned long long int * tempResetd, double rho);
    /*!
     *  \brief  Copie pour raisons de rapidt�  de la valeur de dimension d'�tat
     */
    int dim;
    int dim_c;
    int dim_d;
    /*!
     *  \brief  Copie pour raisons de rapidt�  de la valeur de dimension de contr�le
     */
    int dimC_c;
    int dimC_d;

    int nbPointsCube_c;
    /*!
     *  \brief Pointeur sur la base de donn�es servant � enregister la r�troaction optimale
     */
    string filePrefix;

    unsigned long long int * indicesDecalCell_c;

    /*!
     *  \brief  Une structure servant � stocker les donn�es  de l'image discr�te
     *  d'un point au cours de parcours  de calcul  de \f$ \Phi(C_{n}\setminus C_{n-1}) \f$ .
     *  \see computeDiscreteImageOfPoint
     *  \see computeCurrentImageGlobalRho
     *
     *  Attention! Variable globale dans les m�thodes de la classe!
     */
    discretImageSet_simple pointDI;

    map<unsigned long long int, double> currentCellsImage;
    /*!
     * \brief La liste de structures représentant  des points d'une image en construction
     * Attention! Les méthodes suivantes modifient  cette variable comme étant globale
     * \see createPointsListGlobalRho
     * \see addDataToPointsList
     */
    map<unsigned long long int, double> currentPointsImage;


    unsigned long long int *imageCells;

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

    void addConvexCombinations(unsigned long long int posX, double pointVal, double newCellVal, unsigned long long int numCell, double rho);
    void computeConvexifiedImage(int iter);
    void addDataToCurrentImage(list<unsigned long long int>::iterator *startIt, unsigned long long int newCell,
	    list<unsigned long long int>::iterator *resIt);
    int addNewPoints();
    void computeCurrentImage(int iter);
    void createPointsList();
    void addDataToPointsList(list<unsigned long long int>::iterator *startIt, unsigned long long int newPoint,
	    list<unsigned long long int>::iterator *resIt);

    };

#endif /* SRC_VIABIBITSETHYBRID_H_ */
