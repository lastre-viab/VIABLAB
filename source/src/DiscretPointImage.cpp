/*
 * DiscretPointImage.cpp
 *
 *  Created on: 2 nov. 2024
 *      Author: adesi
 */

#include "../include/DiscretPointImage.h"
DiscretPointImage::DiscretPointImage()
    {
    }

DiscretPointImage::DiscretPointImage(unsigned long long int pos, int d, SysDyn *sd, Grid *gd)
    {
    posX = pos;
    dim = d;
    dynsys = sd;
    grid = gd;
    pointIntCoords = new unsigned long long int[dim];
    pointDoubleCoords = new double[dim];
    grid->numToIntAndDoubleCoords(posX, pointIntCoords, pointDoubleCoords);
    imageCellsWithControls = new map<unsigned long long int, list<unsigned long long int>>();
    _isBuilt = false;
    }

void DiscretPointImage::Build()
    {
    //cout << "Start compute discrete image\n " << std::flush;
    imageCellsWithControls->clear();
    double **controlCoords = dynsys->getControlCoords();
    unsigned long long int nbCTotal = dynsys->getTotalNbPointsC();
    unsigned long long int nbCellsTotal = grid->getNbTotalPoints();
    unsigned long long int intPointCoords[dim];
    double doubleVect1[dim];
    double rho;
    try
    {
	rho = dynsys->calculRho_local(pointDoubleCoords);

	unsigned long long int cu;
	unsigned long long int cellnum;

	//cout << "Calcul discrete image ETAPE 1 \n " << std::flush;
	for (cu = 0; cu < nbCTotal; cu++)
	    {

	    /*!
	     * on calcule les coordonnées réelles de l'image du point par la dynamique discrete
	     * elles sont stockes dans le tableau doubleVect
	     * si \f$ u\in U(x)\f$  (contréle admissible) on calcule l'image réelle par la dynamique
	     *  discrétisée  en temps  du point x avec le controle u
	     */
	    if (dynsys->constraintsXU(pointDoubleCoords, controlCoords[cu]) < PLUS_INF)
		{
		try
		{
		    (dynsys->*(dynsys->discretDynamics))(pointDoubleCoords, controlCoords[cu], doubleVect1, rho);
		}
		catch (...)
		    {
		    std::cout << "Allocation failed dans discrete dyn: " << std::flush;
		    cout.flush();
		    throw;
		    }

		if (grid->isPointInGrid(doubleVect1))
		    {
		    /*!
		     * Si l'image est  dans les limites de la grille on étudie si elle vérifie les contraintes
		     */
		    if (dynsys->constraintsX(doubleVect1) < PLUS_INF)
			{
			/*!
			 * Si l'image  est dans l'ensemble de contraintes sur l'état \f$ K \f$
			 * on calcule le numéro de maille qui contient cett image
			 */
			cellnum = grid->localizePoint(doubleVect1);

			}
		    else
			{
			/*!
			 * Si l'image n'est pas dans  \f$ K \f$ on enregistre un numéro de maille factice qui signifie que
			 * cette image est rejetée
			 */
			cellnum = nbCellsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble
			}
		    }
		else
		    {
		    /*!
		     * Si l'image n'est pas dans  la grille on enregistre un numéro de maille factice qui signifie que
		     * cette image est rejetée
		     */
		    if (grid->unboundedDomain && grid->isPointInGridWithConstr(doubleVect1) && (dynsys->constraintsX(doubleVect1) < PLUS_INF))
			{
			cellnum = nbCellsTotal; // sinon on enregistre un nombre convenu reconnaissanble
			}
		    else
			{
			cellnum = nbCellsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble
			}
		    }
		}
	    else
		{
		/*!
		 * Si l'image n'est pas dans  la grille on enregistre un numéro de maille factice qui signifie que
		 * cette image est rejetée
		 */
		/*!
		 * \todo prévoir la gestion des autorisations de sortie  par axe pour les variables qui peuvent
		 * se trouver dans des intervalles non bornés
		 */
		cellnum = nbCellsTotal + 1; // sinon on enregistre un nombre convenu reconnaissanble

		}
	    if (this->testConstraintesForCell(cellnum))
		{
		if (auto result = imageCellsWithControls->find(cellnum); result != imageCellsWithControls->end())
		    {
		    result->second.push_back(cu);
		    }
		else
		    {
		    list<unsigned long long int> newList =
			    {
				    cu
			    };
		    (*imageCellsWithControls)[cellnum] = newList;
		    }

		}

	    }
	_isBuilt = true;
    }
    catch (...)
	{
	std::cout << "exception dans calcul discrete image " << std::flush;
	cout.flush();
	throw;
	}
   // cout << " fin calcul discrete image " << std::flush;
    }

bool DiscretPointImage::testConstraintesForCell(unsigned long long int numCell)
    {
    bool res = true;
    long long int *indicesDecalCell = grid->getIndicesDecalCell();
    double *tempDoublePointCoords = new double[dim];
unsigned long long int tempIntCoords[dim];
    int nbPointsCube = (int) pow(2.0, dim);	//pow(2.0, dim);
    unsigned long long int pos;
    int i = 0;

    while (res & (i < nbPointsCube))
	{
	pos = numCell + indicesDecalCell[i];
	grid->numToIntAndDoubleCoords(pos, tempIntCoords, tempDoublePointCoords);
	res = (dynsys->constraintsX(tempDoublePointCoords) < PLUS_INF);

	i++;

	}

    delete[] tempDoublePointCoords;

    return res;
    }
map<unsigned long long int, list<unsigned long long int> >* DiscretPointImage::GetImageCellsWithControls()
    {
    return imageCellsWithControls;
    }

void DiscretPointImage::setPointNumAndRebuild(unsigned long long int pointNum)
    {
    posX = pointNum;
    grid->numToIntAndDoubleCoords(posX, pointIntCoords, pointDoubleCoords);
    imageCellsWithControls->clear();
    _isBuilt = false;
    //cout<< " start de build "<<endl;
    Build();
    }

bool DiscretPointImage::IsBuilt()
    {
    return _isBuilt;
    }

DiscretPointImage::~DiscretPointImage()
    {
    delete[] pointDoubleCoords;
    delete[] pointIntCoords;
    imageCellsWithControls->clear();
    delete imageCellsWithControls;
    }

