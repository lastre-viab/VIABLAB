/*
 * GridBitSetHybrid.cpp
 *
 *  Created on: 17 août 2025
 *      Author: adesi
 */

#include "../include/GridBitSetHybrid.h"

GridBitSetHybrid::~GridBitSetHybrid()
    {
    // TODO Auto-generated destructor stub
    }

GridBitSetHybrid::GridBitSetHybrid(const gridParams &gp) :
						Grid_BitSet(gp)
    {
    dim_hc = gp.DIM_HC;
    dim_hd = gp.DIM_HD;
    maxStep = 0;//recompute here the maxStep only on continuous subgrid
    for (int i = 0; i < dim_hc; i++)
	{
	maxStep = max(maxStep, step[i]);
	}
    nbCellPointsContinuousState = pow(2, dim_hc);
  pow3Sub_c = 1;
  for (int p = 0; p < dim_hc - 1; p++)
  {
    pow3Sub_c *= 3;
  }
  indicesDecalSub_c = new long long int[pow3Sub_c];
  int indiceShiftSub = 0;
  unsigned long long int indiceTemp[dim - 1];
  //logVector("[GridBitSetHybrid] : indicesDecalSub ", indicesDecalSub, getPow3Sub());
  unsigned long long int shift = indicesDecalSub[(int) getPow3Sub() - 1];
  //spdlog::debug(" shift {}, pow3sub_c = {}, ",shift, pow3Sub_c);
  for (int k = 0; k < (int) getPow3Sub(); k++)
  {
    numToIntCoords_gen(indicesDecalSub[k] + shift, dim - 1, nbPointsSubGrid, indiceTemp);
    //logVector("[GridBitSetHybrid] : indiceTemp ", indiceTemp, dim - 1);
    bool zeroTest = true;
    for (int i = 0; i < dim_hd; i++) {
      zeroTest &= indiceTemp[dim_hc -1 + i] == 0;
    }
    if (zeroTest) {
      //spdlog::debug(" k {}, indicesDecalSub= {}, ",k, indicesDecalSub[k]);

      indicesDecalSub_c[indiceShiftSub] = indicesDecalSub[k];
      indiceShiftSub++;
    }

  }
    continuousStateShifts = new unsigned long long int[nbCellPointsContinuousState];
    int indiceShift = 0;
	for (int k = 0; k < (int) nbPointsCube; k++)
	{
		bool zeroTest = true;
		for (int i = 0; i < dim_hd; i++) {
			zeroTest &= lesDecalagesCell[k][dim_hc + i] == 0;
		}
		if (zeroTest) {
			continuousStateShifts[indiceShift] = indicesDecalCell[k];
			indiceShift++;
		}

	}
    }

bool GridBitSetHybrid::isValidSubShift(const unsigned long long int *coords, const long long int *shift) const
{

  bool result = true;
  for (int i = 0; i < dim-1 && result; i++)
  {
    result &= ((int) coords[i] + shift[i] < (int)nbPointsSubGrid[i] && 0 <= (int) coords[i] + shift[i]);
    if ( i >= dim_hc - 1)
    {
      result &=  shift[i]  == 0;
    }
  }

  return result;
}


boost::dynamic_bitset<> GridBitSetHybrid::analyseTrameMasqueWithVectorShifts(unsigned long long int posX) const
    {
    int i = 0;
    double x = limInf[0] + step[0];
    unsigned long long int *indice = new unsigned long long int[dim - 1];
    numToIntCoords_gen(posX, dim - 1, nbPointsSubGrid, indice);

    boost::dynamic_bitset<> laTrame(longTrame);
    unsigned long long int iFront;
    unsigned long long int iNext;
    boost::dynamic_bitset<> masqueDecale(longTrame);
    boost::dynamic_bitset<> masque(longTrame);

    bool testBord = false;
    i = 0;
    while ((i < dim_hc - 1) && !testBord)
	{
	testBord = testBord | ((indice[i] == nbPointsSubGrid[i] - 1) | (indice[i] == 0));
	i++;
	}

    if (testBord)
	{
	return ((*gridTab[posX]));
	}
    else
	{
	laTrame = (*gridTab[posX]);
	// on test si la trame n'est pas vide
	if (laTrame.none())
	    {
	    return laTrame;
	    }
	else
	    {
	    int k = 0;
	  bool toInit = true;

		while (k < getPow3Sub())
		    {
		    if (isValidSubShift(indice, neighborShiftsSub[k]))
			{
		      if (toInit)
		      {
		        masque = (((*gridTab[posX + indicesDecalSub[k]])));
		        toInit = false;
		      }
			else
			{
			  masque &= (((*gridTab[posX + indicesDecalSub[k]])));
			}
			//on d�cale la trame de 1 bit et on fait le ET pour chaque
			//bit �a donne : b[i] ET b[i+1]
			masqueDecale = ((*gridTab[posX + indicesDecalSub[k]])) >> (1);
			masqueDecale[longTrame - 1] = 1;
			masque &= (masqueDecale);
			//on d�cale la trame de 1 bite  et on fait le ET pour chaque
			//bit �a donne : b[i] ET b[i-1]
			masqueDecale = ((*gridTab[posX + indicesDecalSub[k]])) << (1);
			masqueDecale[0] = 1;
			masque &= (masqueDecale);
			}
		    k++;
		    }


	    masque ^= (laTrame);
	    masque &= (laTrame);

	    iFront = laTrame.find_first();
	    masque[iFront] = 1;

	    while (iFront < longTrame)
		{
		iNext = laTrame.find_next(iFront);
		while (iNext - iFront == 1 && iNext < longTrame)
		    {
		    iFront = iNext;
		    iNext = laTrame.find_next(iFront);
		    }
		if (iNext < longTrame)
		    {
		    masque[iFront] = 1;
		    masque[iNext] = 1;
		    iFront = iNext;
		    }
		else
		    {
		    masque[iFront] = 1;
		    iFront = iNext;
		    }
		}
	    return (masque);
	    }
	}
    }


boost::dynamic_bitset<> GridBitSetHybrid::analyseTrameMasque(unsigned long long int posX) const
    {
    int i = 0;
    double x = limInf[0] + step[0];
    unsigned long long int *indice = new unsigned long long int[dim - 1];
    numToIntCoords_gen(posX, dim - 1, nbPointsSubGrid, indice);
    boost::dynamic_bitset<> laTrame(longTrame);
    unsigned long long int iFront;
    unsigned long long int iNext;
    boost::dynamic_bitset<> masqueDecale(longTrame);
    boost::dynamic_bitset<> masque;
    //deuxi�me  test : si la trame n'est pas sur la fornti�re de l'espace

    //cout<<"\n";
    bool testBord = false;
    i = 0;
    while ((i < dim - 1) && !testBord)
	{
	testBord = testBord | ((indice[i] == nbPointsSubGrid[i] - 1) | (indice[i] == 0));
	i++;
	}
    delete[] indice;

    if (testBord)
	{
	return ((*gridTab[posX]));
	}
    else
	{
	laTrame = (*gridTab[posX]);
	// on test si la trame n'est pas vide
	if (laTrame.none())
	    { // cout<<"trame vide "<<laTrame;
	    return laTrame;
	    }
	else
	    {
	    int k = 0;
	    while ((indicesDecalSub_c[k] < 0 && posX < ((unsigned long long int) -indicesDecalSub_c[k])) && (k < pow3Sub_c))
		{
		//cout<< " k ="<<k<< " indices dcal="<<indicesDecalSub[k]<<endl;
		k++;
		}
	    if ((k < pow3Sub_c) & (posX + indicesDecalSub_c[k] < nbPointsTotalSubGrid))
		{
		masque = (*gridTab[posX + indicesDecalSub_c[k]]);

		//on d�cale la trame de 1 bite  et on fait le ET pour chaque
		//bit �a donne : b[i] ET b[i+1]
		masqueDecale = ((*gridTab[posX + indicesDecalSub_c[k]])) >> (1);

		masqueDecale[longTrame - 1] = 1;
		masque &= (masqueDecale);

		//on d�cale la trame de 1 bite  et on fait le ET pour chaque
		//bit �a donne : b[i] ET b[i-1]
		masqueDecale = ((*gridTab[posX + indicesDecalSub_c[k]])) << (1);

		masqueDecale[0] = 1;
		masque &= (masqueDecale);

		while (k < pow3Sub_c)
		    {
		    if (indicesDecalSub_c[k] < 0 && posX >= ((unsigned long long int) -indicesDecalSub_c[k]))
			{
			masque &= (((*gridTab[posX + indicesDecalSub_c[k]])));

			//on d�cale la trame de 1 bite  et on fait le ET pour chaque
			//bit �a donne : b[i] ET b[i+1]
			masqueDecale = ((*gridTab[posX + indicesDecalSub_c[k]])) >> (1);
			masqueDecale[longTrame - 1] = 1;
			masque &= (masqueDecale);
			//on d�cale la trame de 1 bite  et on fait le ET pour chaque
			//bit �a donne : b[i] ET b[i-1]
			masqueDecale = ((*gridTab[posX + indicesDecalSub_c[k]])) << (1);
			masqueDecale[0] = 1;
			masque &= (masqueDecale);

			}
		    k++;
		    }
		}

	    masque ^= (laTrame);
	    masque &= (laTrame);

	    iFront = laTrame.find_first();

	    masque[iFront] = 1;

	    while (iFront < longTrame)
		{
		iNext = laTrame.find_next(iFront);
		while (iNext - iFront == 1 && iNext < longTrame)
		    {
		    iFront = iNext;
		    iNext = laTrame.find_next(iFront);
		    }
		if (iNext < longTrame)
		    {
		    masque[iFront] = 1;
		    masque[iNext] = 1;
		    iFront = iNext;
		    }
		else
		    {
		    masque[iFront] = 1;
		    iFront = iNext;
		    }
		}
	    return (masque);
	    }
	}
    }

unsigned long long int* GridBitSetHybrid::GetContinuousStateShifts()
    {
    return continuousStateShifts;
    }

unsigned long long int GridBitSetHybrid::getNbCellsContinuousState()
    {
    return nbCellPointsContinuousState;
    }

bool GridBitSetHybrid::isHybridPointInGrid(double *continuousState, unsigned long long int *discreteState)
    {
    bool isInGrid = true;
    int i = 0;

    while (isInGrid & (i < dim_hc))
	{
	isInGrid &= ((continuousState[i] <= limSup[i]) & (continuousState[i] >= limInf[i]));
	i++;
	}
    if (isInGrid)
	{
	while (isInGrid & (i < dim))
	    {
	    isInGrid &= (discreteState[i - dim_hc] < nbPoints[i]);
	    i++;
	    }
	}

    return isInGrid;
    }

/*!
 * cette fonction corrige les coordonn�es r�elles d'un vecteur
 * qui sont p�riodiques en les ramenant le cas �ch�ant  dans l'intervalle de la p�riode
 */
void GridBitSetHybrid::periodizeHybridPoint(double *vect) const
{
	int i;
	double dist;
	bool testDepass;
	////cout<< " peridisation  vect=";
	//printVector(vect,dim);

	for (i = 0; i < dim_hc; i++)
	{
		/*
		 * On teste  si la variable est d�crar�e p�riodique
		 */
		////cout<< " i= "<< i<< " periodic["<<i<<"]= "<<periodic[i]<<endl;
		if (periodic[i])
		{
			testDepass = true;
			/*
			 * La periode
			 */
			dist = limSup[i] - limInf[i];
			////cout<< " dist= "<<dist<<endl;
			/*
			 * Correction par p�riode enti�re
			 */
			while (testDepass)
			{
				//	//cout<< " test depass="<<testDepass<<endl;
				if (vect[i] > limSup[i])
				{
					vect[i] -= dist;
				}
				else
				{
					if (vect[i] < limInf[i])
					{
						vect[i] += dist;
					}
					else
					{
						testDepass = false;
					}
				}
				//		//cout<< " vect= ";
				//		printVector(vect,dim);
			}
		}
	}
	//	//cout<< "Periodise a fini  vect= ";
	//				printVector(vect,dim);
	////system("pause");
}

unsigned long long int GridBitSetHybrid::LocalizeHybridPoint(double *doubleCoords, unsigned long long int *intCoords)
    {
    unsigned long long int numCell = (unsigned long long int) (((doubleCoords)[0] - limInf[0]) / step[0]);
    if (numCell == (nbPoints[0] - 1))
	numCell--;

    for (int i = 1; i < dim_hc; i++)
	{
	unsigned long long int indiceCell = (unsigned long long int) (((doubleCoords)[i] - limInf[i]) / step[i]);
	if (indiceCell == (nbPoints[i] - 1))
	    {
	    indiceCell--;
	    }
	numCell = numCell * (nbPoints[i]) + indiceCell;
	}

    for (int i = dim_hc; i < dim; i++)
	{
	numCell = numCell * (nbPoints[i]) + intCoords[i - dim_hc];
	}
    return numCell;
    }

void GridBitSetHybrid::saveHybridCoupe(const string &nomFichier) const
    {

    unsigned long long int pos;
    double coordReelles[dim];
    unsigned long long int coordsInt[dim];
    ofstream fichier(nomFichier.c_str());
    if (fichier)  // si l'ouverture a r�ussi
	{
	int compte = 0;
	bool test;
	pos = 0;
	while (pos < nbTotalPoints)
	    {
	    numToIntAndDoubleCoords(pos, coordsInt, coordReelles);
	    if (this->isInSet(coordsInt))
		{
		compte++;
		//si le point appartient � l'ensemble
		//on  va enregistrer ses corrdonn�es r�elles

		test = true;
		for (int k = 0; k < (int) dirs.size(); k++)
		    {
		    if (dirs[k] < dim_hc)
			{
			test = test & (abs(coordReelles[dirs[k]] - values[k]) <= step[dirs[k]]);
			}
		    else
			{
			test = test & (coordsInt[dirs[k]] == values_fd[k]);
			}

		    }
		if (test)
		    {

		    for (int i = 0; i < dim; i++)
			{
			fichier << " " << coordReelles[i] << " ";
			}
		    fichier << "\n";
		    }
		}
	    pos++;
	    }
	//cout<<" nbPoints ecrits="<<compte<<"\n";
	fichier.close();  // je referme le fichier
	}
    else
	// sinon
	cerr << "Erreur � l'ouverture !" << endl;
    //cout<<"fichier fini\n";
    }

void GridBitSetHybrid::saveValOnGridHybrid(const string &fileName) const
    {

    unsigned long long int indice[dim - 1];
    double xCoordsDouble[dim];
    unsigned long long int xCoordsInt[dim];
    ofstream fichierB(fileName.c_str());

    if (fichierB)
	{

	for (unsigned long long int posX = 0; posX < nbPointsTotalSubGrid; posX++)
	    {
	    //cout<< " posx = "<<posX<<endl;
	    numToIntCoords_gen(posX, dim - 1, nbPointsSubGrid, indice);

	    for (int j = 0; j < dirTramage; j++)
		{
		xCoordsDouble[j] = limInf[j] + indice[j] * step[j];
		xCoordsInt[j] = indice[j];
		}

	    for (int j = dirTramage + 1; j < (int) dim; j++)
		{
		xCoordsDouble[j] = limInf[j] + indice[j - 1] * step[j];
		xCoordsInt[j] = indice[j-1];
		}

	    for (unsigned long long int k = 0; k < longTrame; k++)
		{
		xCoordsDouble[dirTramage] = limInf[dirTramage] + k * step[dirTramage];
		xCoordsInt[dirTramage] = k;
		if ((*gridTab[posX])[k])
		    {
		    for (int l1 = 0; l1 < dim_hc; l1++)
			{
			fichierB << xCoordsDouble[l1] << " ";
			}
		    for (int l1 = 0; l1 < dim_hd; l1++)
			{
			fichierB << xCoordsInt[dim_hc + l1] << " ";
			}
		    fichierB << "1 \n";
		    }
		else
		    {
		    for (int l1 = 0; l1 < dim_hc; l1++)
			{
			fichierB << xCoordsDouble[l1] << " ";
			}
		    for (int l1 = 0; l1 < dim_hd; l1++)
			{
			fichierB << xCoordsInt[dim_hc + l1] << " ";
			}
		    fichierB << "0 \n";
		    }
		}
	    }  // fin de for de parcours de la trame
	fichierB.close();
	// je referme le fichier
	}
    else
	// sinon
	cerr << "Erreur � l'ouverture !" << endl;
    }

void GridBitSetHybrid::saveValOnGridLightHybrid(const string &fileName) const
    {
    unsigned long long int indice[dim - 1];
    double xCoordsDouble[dim];
    unsigned long long int xCoordsInt[dim];
    ofstream fichierB(fileName.c_str());

    if (fichierB)
	{

	for (unsigned long long int posX = 0; posX < nbPointsTotalSubGrid; posX++)
	    {
	    //cout<< " posx = "<<posX<<endl;
	    numToIntCoords_gen(posX, dim - 1, nbPointsSubGrid, indice);

	    for (int j = 0; j < dirTramage; j++)
		{
		xCoordsDouble[j] = limInf[j] + indice[j] * step[j];
		xCoordsInt[j] = indice[j];
		}

	    for (int j = dirTramage + 1; j < (int) dim; j++)
		{
		xCoordsDouble[j] = limInf[j] + indice[j - 1] * step[j];
		xCoordsInt[j] = indice[j-1];
		}

	    for (unsigned long long int k = 0; k < longTrame; k++)
		{
		xCoordsDouble[dirTramage] = limInf[dirTramage] + k * step[dirTramage];
		xCoordsInt[dirTramage] = k;
		if ((*gridTab[posX])[k])
		    {
		    for (int l1 = 0; l1 < dim_hc; l1++)
			{
			fichierB << xCoordsDouble[l1] << " ";
			}
		    for (int l1 = 0; l1 < dim_hd; l1++)
			{
			fichierB << xCoordsInt[dim_hc + l1] << " ";
			}
		    fichierB << "1 \n";
		    }
		}
	    }  // fin de for de parcours de la trame
	fichierB.close();
	// je referme le fichier
	}
    else
	// sinon
	cerr << "Erreur � l'ouverture !" << endl;
    }

void GridBitSetHybrid::loadSetHybrid(const string &fileName)
    {
    string line;
    ifstream userDataFile(fileName);
    double xCoords[dim_hc];
    bool val;
    unsigned long long int xCoordsInt[dim];
    istringstream sstr;
    int d = 0;

    if (!userDataFile)
	{
	spdlog::error("Error loading set, could not open set file : {}", fileName);
	}
    else
	{
	while (getline(userDataFile, line))
	    {
	    sstr.str(line);
	    for (d = 0; d < dim_hc; d++)
		{
		sstr >> xCoords[d];
		// Transformation inverse des coordonnées en double vers les coordonnées entières
		const double xCoordsIntI = (xCoords[d] - limInf[d]) / step[d] - 0.5 * arePointsGridCenters;
		xCoordsInt[d] = llround(xCoordsIntI);
		}
	    for(d = dim_hc; d< dim; d++)
		{
		sstr >> xCoordsInt[d];
		}
	    sstr >> val;

	    unsigned long long int posX = this->intCoordsToNum_dm1(xCoordsInt);
	    (*gridTab[posX])[xCoordsInt[dirTramage]] = val;
	    }
	}
    }
