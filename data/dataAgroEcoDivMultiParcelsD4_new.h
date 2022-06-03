/*
 * dataGAIA_phases1.0-2.h
 *
 *  Created on: 23 niv. 2015
 *      Author: Anna D�silles
 *
 *      Ce fichier d�crit le mod�le pour le calcul de la viabilit�
 *      d'une exploitation � une seule parcelle avec les contraintes budg�taires
 *      intra-cycle
 */

#ifndef DATAGAIA_MODEL1_H_
#define DATAGAIA_MODEL1_H_



#include "agroEcoDiv_params_defs.h"
/*!
 *  Ce mod�le d�crit l'exploitation d'une seule parcelle
 *   avec une variable suppl�mentatire : la saison
 *   La variable de saison permet de contraindre les choix de culture � planter
 *   en fonction de la p�riode de l'ann�e
 *   La nouvelel variable prend dont 12 valeurs correspondant aux 12 mois de l'ann�e
 */



string paramsFile = "agroEcoDiv_D4.json";

void optimTrajConvert(string fileNameBis, string fileName, int W2, double lScale);

unsigned long long int numParcelle =1;

unsigned long long int nbParcels=1;
double surf;
double cFixe;

void triBulle(unsigned long long int * permut, double * tab, unsigned long long int dim);



double * surfaces;
unsigned long long int TT;


inline void dynamics_fd(unsigned long long int  * x, unsigned long long int *u, unsigned long long int * image)
{

	// x[0]=ibqs1
	//x[1] = Ibqs2
	//x[2] = gamma : temps restant sur la parcelle inactive
	//x[3] = p : parcelle active : 0 si ( I1 actif, I2 passif)  1 si les deux actifs, 2 si (I2 actif et I1 passif)
	//x[3] = tau : temps

	//u[0] = cycle de parcelle 1
	//u[1]=cycle de parcelle 2
	// v : le controle actif

	unsigned long long int gamma, gammaNew ;
	unsigned long long int p , pNew;

	if(x[2] == 0)
	{
		p=1;
		gamma = 0;
	}
	else
	{
		if( x[2] % 2 == 0)
		{
			p =2;
			gamma = x[2]/2;
		}
		else
		{
			p=0;
			gamma = (x[2] +1)/2;
		}
	}
	unsigned long long int delta1, delta2;
	if(x[3]<TT)
	{
		switch(p)
		{
		case 0://I1 actif, I2 inactif
			delta1 = modelData->getCycleDuration(u[0]);

			image[0]=  modelData->getItransition(x[0], x[3], u[0]) ;
			// cout<< " image 0 = "<<image[0]<<endl;
			image[1]= x[1];
			gammaNew = abs((int) delta1 - (int) gamma);
			pNew = 1 + sign(delta1 - gamma);
			switch(pNew)
			{
			case 0: image[2] = 2*gammaNew-1; break;
			case 1: image[2] = 0; break;
			case 2: image[2] = 2*gammaNew; break;
			}
			image[3] = min(TT,x[3] + min(delta1, gamma));
			break;
			case 1://les deux actifs
				delta1 = modelData->getCycleDuration(u[0]);
				delta2 = modelData->getCycleDuration(u[1]);
				image[0]=  modelData->getItransition(x[0], x[3], u[0]) ;
				// cout<< " image 0 = "<<image[0]<<endl;
				image[1]= modelData->getItransition(x[1], x[3], u[1]) ;
				gammaNew = abs((int) delta1 - (int) delta2);
				pNew = 1 + sign(delta1 - delta2);
				switch(pNew)
				{
				case 0: image[2] = 2*gammaNew-1; break;
				case 1: image[2] = 0; break;
				case 2: image[2] = 2*gammaNew; break;
				}
				image[3] =  min(TT, x[3] + min(delta1, delta2) );
				break;
				case 2://I1 inactif, I2 actif
					delta2 = modelData->getCycleDuration( u[1]);
					image[0]=  x[0];
					// cout<< " image 0 = "<<image[0]<<endl;
					image[1]= modelData->getItransition(x[1], x[3], u[1]) ;
					gammaNew = abs((int) gamma - (int) delta2);
					pNew = 1 + sign(gamma - delta2);
					switch(pNew)
					{
					case 0: image[2] = 2*gammaNew-1; break;
					case 1: image[2] = 0; break;
					case 2: image[2] = 2*gammaNew; break;
					}
					image[3] = min(TT, x[3] + min(gamma, delta2) );
					break;

		}
		//	 cout<< " image 1 = "<<image[1]<<endl;
	}
	else
	{
		image[0]=x[0];
		image[1]=x[1];
		image[2]=x[2];
		image[3]=x[3];
	}
}
inline double l_fd(unsigned long long int  * x, unsigned long long int * u )
{
	double res;
	int s,  saisonalite1, saisonalite2, controlInactifCompatible;
	// x[0]=ibqs1
	//x[1] = Ibqs2
	//x[2] = gamma : temps restant sur la parcelle inactive
	//x[3] = p : parcelle active : 0 si ( I1 actif, I2 passif)  1 si les deux actifs, 2 si (I2 actif et I1 passif)
	//x[3] = tau : temps

	//u[0] = cycle de parcelle 1
	//u[1]=cycle de parcelle 2
	// v : le controle actif

	unsigned long long int gamma;
	unsigned long long int p;

	if(x[2] == 0)
	{
		p=1;
		gamma = 0;
		//cout<< " x[2] = "<< x[2] << " p = "<< p << " gamma = "<< gamma << endl;
	}
	else
	{
		if( x[2] % 2 == 0)
		{
			p =2;
			gamma = x[2]/2;
			//cout<< " x[2] = "<< x[2] << " p = "<< p << " gamma = "<< gamma << endl;
		}
		else
		{
			p=0;
			gamma = (x[2] +1)/2;
			//cout<< " x[2] = "<< x[2] << " p = "<< p << " gamma = "<< gamma << endl;
		}
	}

	unsigned long long int delta1, delta2;
	unsigned long long int mStart1, mEnd1, mStart2, mEnd2;
	unsigned long long int previousI;
	double bilan1, bilan2, bilan;
	if(x[3]<TT)
	{
		//cout<< " calcul de l_fd x[3] = "<< x[3]<<endl;
		//cout<< " controle "<<u[0]<<" "<< u[1]<<endl;
		switch(p)
		{
		case 0://I1 actif => u0 est soumis � la saisonnalit� et au d�ficit, u1 est soumis � la compatibilit� de dur�e avec le temps restant
			delta1 = modelData->getCycleDuration(u[0]);
			delta2 = modelData->getCycleDuration(u[1]);


			mStart1 = 1;
			mEnd1 = mStart1+min(delta1, gamma)-1;
			//cout<< " mStart1 = "<<mStart1<< " mEnd1 "<< mEnd1<<endl;
			bilan1 = modelData->getCycleBalanceForGivenPeriod(x[0], mStart1, mEnd1, u[0]);
			//cout<< " bilan 1 = "<<bilan1<<endl;

			mStart2 = delta2 - gamma +1;
			mEnd2 = mStart2+min(delta1, gamma)-1;
			//cout<< " mStart2 = "<<mStart1<< " mEnd2 "<< mEnd1<<endl;
			//cout<< " tau  = "<<x[3]<< "  I = "<< x[1] << " u1 = "<<u[1]<<endl;
			previousI = modelData->getItransitionPrevious(x[1], x[3]+gamma - delta2, u[1]);

			//cout << " CurrentI on parcel 2 : "<<x[1]<< "  previous I "<<previousI<<endl;
			if(previousI > 1000) exit(1);
			bilan2 = modelData->getCycleBalanceForGivenPeriod(previousI, mStart2, mEnd2, u[1]);
			//cout<< " bilan 2 = "<<bilan2<<endl;
			res = bilan1 * surfaces[0] + bilan2 * surfaces[1];

			break;
		case 1://les deux actifs
			delta1 = modelData->getCycleDuration(u[0]);
			delta2 = modelData->getCycleDuration(u[1]);
			// cout<< " delta1 = "<<delta1<< " delta2 "<< delta2<<endl;
			mStart1 = 1;
			mEnd1 = mStart1+min(delta1, delta2)-1;
			//cout<< " mStart1 = "<<mStart1<< " mEnd1 "<< mEnd1<<endl;
			bilan1 = modelData->getCycleBalanceForGivenPeriod(x[0], mStart1, mEnd1, u[0]);

			mStart2 = 1;
			mEnd2 = mStart2+min(delta1, delta2)-1;
			bilan2 = modelData->getCycleBalanceForGivenPeriod(x[1], mStart2, mEnd2, u[1]);

			res = bilan1 * surfaces[0] + bilan2 * surfaces[1];
			break;
		case 2://I2 actif
			delta1 = modelData->getCycleDuration(u[0]);
			delta2 = modelData->getCycleDuration(u[1]);

			mStart1 = delta1 - gamma +1;
			mEnd1 = mStart1+min(delta2, gamma)-1;
			previousI = modelData->getItransitionPrevious(x[0], x[3]+gamma - delta1, u[0]);
			//cout << " CurrentI on parcel 1 : "<<x[0]<< " previous I "<<previousI<<endl;
			bilan1 = modelData->getCycleMaxDeficitForGivenPeriod(previousI, mStart1, mEnd1, u[0]);

			mStart2 = 1;
			mEnd2 = mStart2+min(delta2, gamma)-1;
			bilan2 = modelData->getCycleMaxDeficitForGivenPeriod(x[1], mStart2, mEnd2, u[1]);

			res = bilan1 * surfaces[0] + bilan2 * surfaces[1];

			break;

		}
	}
	else
	{
		res=0.0;
	}
	return  res;
}

inline double  constraintsX_fd( unsigned long long int * x )
{
	double res;
	if( (x[3]>=TT))
	{
		if( (IRealValues[x[0]]>=IStar) && (IRealValues[x[1]]>=IStar))
		{
			// cout<< " init constrainte  x="<<x[0]<<" "<<x[1]<<endl;
			res= 0.0;//-deficitMax;//0.0;//minMaxDeficit[x[0]];//cFixe*x[1];
		}
		else
		{
			res=PLUS_INF;
		}
	}
	else
	{
		res=0.0-deficitMax;
	}
	return  res;
}

/*!
 * \brief Function  defining the mixed  constraints
 *
 * This function defines the set U(x) for admissible controls as function of the state
 * @param x state variable
 * @param u control variable
 * @return  value that caraterise the constraints set
 */
inline double  constraintsXU_fd( unsigned long long int * x, unsigned long long int * u )
{
	unsigned long long int gamma;
	unsigned long long int p;

	if(x[2] == 0)
	{
		p=1;
		gamma = 0;
		//cout<< " x[2] = "<< x[2] << " p = "<< p << " gamma = "<< gamma << endl;
	}
	else
	{
		if( x[2] % 2 == 0)
		{
			p =2;
			gamma = x[2]/2;
			//	cout<< " x[2] = "<< x[2] << " p = "<< p << " gamma = "<< gamma << endl;
		}
		else
		{
			p=0;
			gamma = (x[2] +1)/2;
			//cout<< " x[2] = "<< x[2] << " p = "<< p << " gamma = "<< gamma << endl;
		}
	}
	int s=(x[3])%12;
	int  saisonalite1, saisonalite2, controlInactifCompatible;
	// x[0]=ibqs1
	//x[1] = Ibqs2
	//gamma = gamma : temps restant sur la parcelle inactive
	//x[3] = p : parcelle active : 0 si ( I1 actif, I2 passif)  1 si les deux actifs, 2 si (I2 actif et I1 passif)
	//x[3] = tau : temps

	//u[0] = cycle de parcelle 1
	//u[1]=cycle de parcelle 2
	// v : le controle actif
	unsigned long long int delta1, delta2;

	double res;
	switch(p)
	{
	case 0://I1 actif => u0 est soumis � la saisonnalit� et au d�ficit, u1 est soumis � la compatibilit� de dur�e avec le temps restant
		delta1 = modelData->getCycleDuration(u[0]);
		delta2 = modelData->getCycleDuration(u[1]);
		saisonalite1=modelData->getCycleSeasonality(u[0], s);
		controlInactifCompatible = (delta2 >= gamma)?1:0;


		res=PLUS_INF*( (1.0-saisonalite1) + (1.0 - controlInactifCompatible ) )+1.0*(saisonalite1) * controlInactifCompatible;
		break;
	case 1://les deux actifs
		delta1 = modelData->getCycleDuration(u[0]);
		delta2 = modelData->getCycleDuration(u[1]);
		saisonalite1=modelData->getCycleSeasonality(u[0], s);
		saisonalite2=modelData->getCycleSeasonality(u[1], s);

		res=PLUS_INF*( (1.0-saisonalite1) + (1.0 - saisonalite2 ) )+1.0*(saisonalite1) * saisonalite2;

		break;
	case 2://I2 actif
		delta1 = modelData->getCycleDuration(u[0]);
		delta2 = modelData->getCycleDuration(u[1]);
		saisonalite2=modelData->getCycleSeasonality(u[1], s);
		controlInactifCompatible = (delta1 >= gamma)?1:0;
		res=PLUS_INF*( (1.0-saisonalite2) + (1.0 - controlInactifCompatible ) )+1.0*(saisonalite2) * controlInactifCompatible;

		break;

	}

	return res;
}
inline double  constraintsXUY_fd( unsigned long long int * x, unsigned long long int * u )
{

	unsigned long long int gamma;
	unsigned long long int p;

	if(x[2] == 0)
	{
		p=1;
		gamma = 0;
		//cout<< " x[2] = "<< x[2] << " p = "<< p << " gamma = "<< gamma << endl;
	}
	else
	{
		if( x[2] % 2 == 0)
		{
			p =2;
			gamma = x[2]/2;
			//cout<< " x[2] = "<< x[2] << " p = "<< p << " gamma = "<< gamma << endl;
		}
		else
		{
			p=0;
			gamma = (x[2] +1)/2;
			//cout<< " x[2] = "<< x[2] << " p = "<< p << " gamma = "<< gamma << endl;
		}
	}

	double res;
	int s,  saisonalite1, saisonalite2, controlInactifCompatible;
	// x[0]=ibqs1
	//x[1] = Ibqs2
	//x[2] = gamma : temps restant sur la parcelle inactive
	//x[3] = p : parcelle active : 0 si ( I1 actif, I2 passif)  1 si les deux actifs, 2 si (I2 actif et I1 passif)
	//x[3] = tau : temps

	//u[0] = cycle de parcelle 1
	//u[1]=cycle de parcelle 2
	// v : le controle actif
	unsigned long long int delta1, delta2;
	unsigned long long int mStart1, mEnd1, mStart2, mEnd2;
	unsigned long long int previousI;
	unsigned long long int  mStart[2], IStart[2];
	double maxDeficit1, maxDeficit2, maxDeficit;

	if(x[3]<TT)
	{
		s=(x[3])%12;

		switch(p)
		{
		case 0://I1 actif => u0 est soumis � la saisonnalit� et au d�ficit, u1 est soumis � la compatibilit� de dur�e avec le temps restant
			delta1 = modelData->getCycleDuration(u[0]);
			delta2 = modelData->getCycleDuration(u[1]);
			saisonalite1=modelData->getCycleSeasonality(u[0], s);
			controlInactifCompatible = (delta2 >= gamma)?1:0;

			mStart[0] = 1;
			IStart[0] = x[0];

			mStart[1] = delta2 - gamma +1;
			IStart[1] = modelData->getItransitionPrevious(x[1], x[3]+gamma - delta2, u[1]);

			maxDeficit =  modelData->getCycleMaxDeficitForGivenPeriodMultiParcel(IStart, mStart, min(delta1, gamma), u);

			res=PLUS_INF*( (1.0-saisonalite1) + (1.0 - controlInactifCompatible ) )+(maxDeficit-deficitMax)*(saisonalite1) * controlInactifCompatible;
			break;
		case 1://les deux actifs
			delta1 = modelData->getCycleDuration(u[0]);
			delta2 = modelData->getCycleDuration(u[1]);
			saisonalite1=modelData->getCycleSeasonality(u[0], s);
			saisonalite2=modelData->getCycleSeasonality(u[1], s);

			mStart[0] = 1;
			IStart[0] = x[0];
			mStart[1] = 1;
			IStart[1] = x[1];

			maxDeficit = modelData->getCycleMaxDeficitForGivenPeriodMultiParcel(IStart, mStart, min(delta1, delta2), u);;
			res=PLUS_INF*( (1.0-saisonalite1) + (1.0 - saisonalite2 ) )+(maxDeficit-deficitMax)*(saisonalite1) * saisonalite2;

			break;
		case 2://I2 actif
			delta1 = modelData->getCycleDuration(u[0]);
			delta2 = modelData->getCycleDuration(u[1]);
			saisonalite2=modelData->getCycleSeasonality(u[1], s);
			controlInactifCompatible = (delta1 >= gamma)?1:0;


			mStart[0] = 1;
			IStart[0] = modelData->getItransitionPrevious(x[0], x[3]+gamma - delta1, u[0]);

			mStart[1] = 1;
			IStart[1] = x[1];

			maxDeficit = modelData->getCycleMaxDeficitForGivenPeriodMultiParcel(IStart, mStart, min(gamma, delta2), u);;
			res=PLUS_INF*( (1.0-saisonalite2) + (1.0 - controlInactifCompatible ) )+(maxDeficit-deficitMax)*(saisonalite2) * controlInactifCompatible;

			break;

		}
	}

	else
	{
		res=0.0-deficitMax;
	}
	//double res=PLUS_INF*(1.0-saisonalite[u[0]][u[1]][s])+0.0*(saisonalite[u[0]][u[1]][s]);

	//cout<< "  mu func  res = "<<res<<endl;
	return res;
}



/*!
 * \brief Function  defining the state   constraints, corresonds  to k(x)
 *
 * This function defines the set K for admissible  states
 * @param x state variable
 * @return  value that caraterise the constraints set
 */

void loadModelData(ParametersManager *PM)
{
	ostringstream os;

	int i,j,k,l,xx,yy,zz;

	int nbCosts=7;
	double costs[nbCosts];

	double val;
	FILE * iFile,   *firstDataFile;
	string line;
	ifstream * iDataFile= new ifstream();

	istringstream sstr;

	/*   *****************************************************************************
	 * Lecture des params du projet
	 * *******************************************************************************
	 */
	gridParams * gp = PM->getGridParameters();
	algoViabiParams * avp = PM->getAlgoParameters();
	controlParams * cp = PM->getControlParameters();
	systemParams * sysp = PM->getSystemParameters();



	string firstDataFileName = "SimulationGaia.json";


	modelData=new FarmData_Gardener(firstDataFileName);
	deficitMax = modelData->getAdmissibleDeficit();
	sysp->maxTime=(double) (modelData->getTimeHorizon());
	TT = modelData->getTimeHorizon();
	/*!
	 * \var nbPointsControl[dimc]
	 * number of  discretization points for the control variable
	 * \see controlPa;rams
	 */
	cp->NBPOINTSC[0] =  modelData->getNbControls();
	cp->NBPOINTSC[1] =  modelData->getNbControls();
	/*!
	 * \var nbPointsState[dim]
	 * number of  discretization points for the state variable
	 * \see gridParams
	 */

	gp->NBPOINTS[0]    =  modelData->getNbIPoints();
	gp->NBPOINTS[1]    =  modelData->getNbIPoints();
	gp->NBPOINTS[2]    =  1+2*modelData->getMaxCycleDuration()+1;
	gp->NBPOINTS[3]    =  modelData->getTimeHorizon()+1;


	gp->LIMINF[0] = 0.0;
	gp->LIMINF[1] = 0.0;
	gp->LIMINF[2] = 0.0;
	gp->LIMINF[3] = 0.0;


	gp->LIMSUP[0] = 1.0;
	gp->LIMSUP[1] = 1.0;
	gp->LIMSUP[2] = (double) (1+2*modelData->getMaxCycleDuration());
	gp->LIMSUP[3] = (double) modelData->getTimeHorizon();

	IStar=modelData->getIStar();
	IRealValues = modelData->getIgrid();
	prefix = modelData->getPrefix();

	nbParcels = modelData->getNbParcels();
	surfaces = modelData->getParcelsSurfaces();

	cout<< " fini loadModelData T = "<<sysp->maxTime<<endl;
	int dim = gp->DIM;
	nbTrajs = modelData->getNbTrajs();
	avp->INIT_POINTS_FD = new unsigned long long int [dim * nbTrajs];
	avp->INIT_VALUES_FD = modelData->getInitValues();
	unsigned long long int *tempinitPoints = modelData->getInitPointsCoords();

	for(int j = 0; j<nbTrajs; j++)
	{
		avp->INIT_POINTS_FD[j*dim + 0] = tempinitPoints[2*j + 0];
		avp->INIT_POINTS_FD[j*dim + 1] = tempinitPoints[2*j + 0];
		avp->INIT_POINTS_FD[j*dim + 2] = 0;
		avp->INIT_POINTS_FD[j*dim + 3] = tempinitPoints[2*j + 1];

	}
}
//system("pause");

void postProcess();
void postProcess()
{
	cout<< " post process nettoyage memoire \n";



}


inline double constraintsX( double * x )
{

	// x[0]=ibqs
	//x[1]=t
	bool res1;

	if( (x[3]>=TT))
		res1= (x[0]>=IStar && x[1]>=IStar);
	else
		res1=true;
	return (double)res1;
}


inline  double target_fd (unsigned long long int * x)
{


	double res;

	if( (x[3]>=TT))
	{
		if(IRealValues[x[0]]>=IStar && IRealValues[x[1]]>=IStar)
		{
			//cout<< " init constrainte  x="<<x[0]<<" "<<x[1]<<endl;
			res= 0.0;//-deficitMax;//0.0;//minMaxDeficit[x[0]];//cFixe*x[1];
		}
		else
		{
			res=PLUS_INF;
		}
	}
	else
	{
		res=PLUS_INF;
	}
	return  res;
}


/* *   *************************************************************************
 *      Definition of value function
 ***************************************************************************/





double interpolationLin(int nbDataPoints, double *tabT, double *tabX,double t){
	int i=1;
	double res;

	if (t<tabT[0])
	{
		return(tabX[0]);
	}
	while (i<nbDataPoints && tabT[i]<t)
	{
		i++;
	}
	if (i==nbDataPoints)
	{
		return(tabX[nbDataPoints-1]);
	}
	res=tabX[i-1] + ((t - tabT[i-1])/(tabT[i] - tabT[i-1])) * (tabX[i] - tabX[i-1]);

	return(res);
}

double interpolationLinPhase(int nbDataPoints, double *tabT, double *tabX,double t){
	/*
	 * Cette fonction est r�serv�e � l'interpolation du cap
	 * � partir du cap r�eel, \phi \in\R.
	 * La fonction interpolle la valeur du cap et le ram�ne � [0,2pi]
	 */
	int i=1;
	double res;

	if (t<tabT[0])
	{
		res=tabX[0];
	}
	while (i<nbDataPoints && tabT[i]<t)
	{
		i++;
	}
	if (i==nbDataPoints)
	{
		res=tabX[nbDataPoints-1];
	}
	else
	{
		res=tabX[i-1] + ((t - tabT[i-1])/(tabT[i] - tabT[i-1])) * (tabX[i] - tabX[i-1]);
	}
	while(res>2.0*pi)
	{
		res=res-2.0*pi;
	}
	while(res<0)
	{
		res=res+2.0*pi;
	}
	return(res);
}



void triBulle(unsigned long long int * permut, double * tab, unsigned long long int dim)
{
	for(unsigned long long int i=0;i<dim;i++)
	{
		permut[i]=i;
	}
	double tmp;
	unsigned long long int tmpi;
	unsigned long long int cpt=1;
	while(cpt>0)
	{
		cpt=0;
		//cout<< " tab = ";
		//printVector(tab,dim);
		//cout<< " permut= "; printVector(permut,dim);
		for(unsigned long long int i=0;i<dim-1;i++)
		{
			if(tab[i]>tab[i+1])
			{
				tmp=tab[i];
				tab[i]=tab[i+1];
				tab[i+1]=tmp;

				tmpi=permut[i];
				permut[i]=permut[i+1];
				permut[i+1]=tmpi;
				cpt++;
			}
		}
		//cout<<" cpt de permut = "<<cpt<< " tab = ";
		//printVector(tab,dim);
		//cout<< " permut= "; printVector(permut,dim);
	}

}

#endif /* DATAGAIA_MODEL1_H_ */
