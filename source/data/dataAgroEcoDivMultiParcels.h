/*
 * dataGAIA_phases1.0-2.h
 *
 *  Created on: 23 niv. 2015
 *      Author: Anna Dï¿½silles
 *
 *      Ce fichier dï¿½crit le modï¿½le pour le calcul de la viabilitï¿½
 *      d'une exploitation ï¿½ une seule parcelle avec les contraintes budgï¿½taires
 *      intra-cycle
 */

#ifndef DATAGAIA_MODEL1_H_
#define DATAGAIA_MODEL1_H_



#include "agroEcoDiv_params_defs.h"
/*!
 *  Ce modï¿½le dï¿½crit l'exploitation d'une seule parcelle
 *   avec une variable supplï¿½mentatire : la saison
 *   La variable de saison permet de contraindre les choix de culture ï¿½ planter
 *   en fonction de la pï¿½riode de l'annï¿½e
 *   La nouvelel variable prend dont 12 valeurs correspondant aux 12 mois de l'annï¿½e
 */

/*! \var dim
 *  \brief State  dimension
 */
const int dim=5;
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=2;
/*! \var dicret_type
 *  \brief the  type of discretization scheme used  for the dynamics
 *  EE or 1 = Euler Explicit  scheme
 *  EI or 2 = Euler Implicit
 *  RK2I or 3 = RK2 Implicit (RK2 for -F)
 */
const int discret_type=RK2I;

/*!
 * \var T maximum time horizon for the study
 */
double T=3*12;

/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MIN[dimc]={0,0};
/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX[dimc]={1,1};

/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */
double STATE_MINreel[dim]={0.0,0.0,0.0,0.0,0.0};
double Treel;
/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAXreel[dim]={1.0, 10.0, 10.0, 10.0, 10.0};
double step[dim];


/*
 * le vecteur scaling indique pour chaque variable s'il faut ou pas mettre Ã  l'Ã©chelle
 *   1 = mise Ã  l'Ã©chelle
 *   0 = pas de mise Ã  l'Ã©chelle
 */
const int scaling[dim+1]          = {0,0,0,0,0,0};    //- SCALED VARIABLES (1:scaled, 0:not scaled)



void optimTrajConvert(string fileNameBis, string fileName, int W2, double lScale);

unsigned long long int numParcelle =1;

unsigned long long int nbParcels=1;


/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */

double  STATE_MIN[dim]         = {0,0,0,0,0};
/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */

double  STATE_MAX[dim]         = {1,1,1,1,1};


double  STATE_MIN_V[dim+1]         = {0,0,0,0,0,0};
/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */

double  STATE_MAX_V[dim+1]         = {1.0,1.0,1.0,1.0,1.0,1.0};



/*
 * target = 1
 * departure =0;
 * Ce parametre determine le sens des trajectoires
 */
int target_or_departure_problem=1;
int setType=VIAB;
int compute_tmin= 0;


/*!
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */
unsigned long long int nbPointsControl[dimc] =  {1,1};
/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int nbPointsState[dim]    =  {10,10,(unsigned long long int)T+1,(unsigned long long int)T+1,(unsigned long long int)T+1};
unsigned long long int nbPointsStateV[dim+1]    =  {10,(unsigned long long int)T+1,10,(unsigned long long int)T+1,(unsigned long long int)T+1,(unsigned long long int)T+1};

double surf;
double cFixe;






/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0,0,0,0};
int periodicV[dim+1]={0,0,0,0,0,0};


int compteOnlyOptimalRetro=1;

int computeSet=1;
/*
 * sÃ©lection de type de reconstruction de trajectoire
 *  VD= 1, viable par dÃ©faut, on sÃ©lectionne le premier contrÃ´le viable trouvÃ©
 *  VL= 2, viable lourd: le contrÃ´le reste constant tant qu'il viable;
 *  cette mÃ©thode nÃ©cessite une initialisation de contrÃ´le
 */
int typeTraj=VD;

double initPoint[dim]={ 1.5, 3.8, 3.8 , 3.8 , 3.8  };
double initPointV[dim+1]={ 1.5, 3.8, 9.0, 3.8 , 3.8 , 3.8  };

unsigned long long int initPoint_fd_V[dim+1]={12,0,0,0,0,0};

void triBulle(unsigned long long int * permut, double * tab, unsigned long long int dim);

unsigned long long int projection[dim]={0,0,0,0,0};

unsigned long long int trajProjection[dim]={0,0,0,0,0};

unsigned long long int dirTramage=0;


/*!
 * \var globalDeltaT
 *  boolï¿½en indique si le pas de temps  doit ï¿½re choisi globalement
 *  ou localement pour les algorithmes de viabilitï¿½
 */
bool globalDeltaT=false;

int ompThreads=1;


void dynamics_fd(unsigned long long int  * x, unsigned long long int *u, unsigned long long int * image);


/*!
 * \var dynType  defines the type of dynamics in the model
 *      1 or CC : continuous
 *      2 or DD : discrete
 *      3 or HD : hybrid
 *      4 or FDD full discret
 */
int dynType=DD;
/*!
 * \var fd_type  defines the type of representation for full discret dynamics in the model
 *      1 or FUNC : usign dynamics_fd function
 *      2 or RETRO : using exaustive representation in the RETRO action file
 *
 */
int fd_type=FUNC;

void dynamics_fd(unsigned long long int  * x, unsigned long long int *u, unsigned long long int * image)
{

	// x[0]=ibqs1
	//x[1] = Ibqs2
	//x[2] = gamma : temps restant sur la parcelle inactive
	//x[3] = p : parcelle active : 0 si ( I1 actif, I2 passif)  1 si les deux actifs, 2 si (I2 actif et I1 passif)
	//x[4] = tau : temps

	//u[0] = cycle de parcelle 1
	//u[1]=cycle de parcelle 2
	// v : le controle actif
	unsigned long long int delta1, delta2;
	if(x[4]<T)
	{
		switch(x[3])
		{
		case 0:
			delta1 = modelData->getTimeTransition(x[4], u[0]);

			image[0]=  modelData->getItransition(x[0], x[4], u[0]) ;
			// cout<< " image 0 = "<<image[0]<<endl;
			image[1]= x[1];
			image[2] = abs((int) delta1 - (int) x[2]);
			image[3] = 1 + sign(delta1 - x[2]);
			image[4] = x[4] + min(delta1, x[2]);
			break;
		case 1:
			delta1 = modelData->getTimeTransition(x[4], u[0]);
			delta2 = modelData->getTimeTransition(x[4], u[1]);
			image[0]=  modelData->getItransition(x[0], x[4], u[0]) ;
			// cout<< " image 0 = "<<image[0]<<endl;
			image[1]= modelData->getItransition(x[1], x[4], u[1]) ;
			image[2] = abs((int) delta1 - (int) delta2);
			image[3] = 1 + sign(delta1 - delta2);
			image[4] = x[4] + min(delta1, delta2);
			break;
		case 2:
			delta2 = modelData->getTimeTransition(x[4], u[1]);
			image[0]=  modelData->getItransition(x[0], x[4], u[0]) ;
			// cout<< " image 0 = "<<image[0]<<endl;
			image[1]= modelData->getItransition(x[1], x[4], u[1]) ;
			image[2] = abs((int) x[2] - (int) delta2);
			image[3] = 1 + sign(x[2] - delta2);
			image[4] = x[4] + min(x[2], delta2);
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
		image[4]=x[4];
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
	//x[4] = tau : temps

	//u[0] = cycle de parcelle 1
	//u[1]=cycle de parcelle 2
	// v : le controle actif
	unsigned long long int delta1, delta2;
	unsigned long long int mStart1, mEnd1, mStart2, mEnd2;
	double bilan1, bilan2, bilan;
	if(x[4]<T)
	{
		//cout<< " calcul de l_fd x[3] = "<< x[3]<<endl;
		//cout<< " controle "<<u[0]<<" "<< u[1]<<endl;
		switch(x[3])
		{
		case 0://I1 actif => u0 est soumis à la saisonnalité et au déficit, u1 est soumis à la compatibilité de durée avec le temps restant
			delta1 = modelData->getTimeTransition(x[4], u[0]);
			delta2 = modelData->getTimeTransition(x[4], u[1]);


			mStart1 = 1;
			mEnd1 = mStart1+min(delta1, x[2])-1;
			//cout<< " mStart1 = "<<mStart1<< " mEnd1 "<< mEnd1<<endl;
			bilan1 = modelData->getCycleBalanceForGivenPeriod(x[0], mStart1, mEnd1, u[0]);
//cout<< " bilan 1 = "<<bilan1<<endl;
//cout<< " mStart2 = "<<mStart1<< " mEnd2 "<< mEnd1<<endl;
			mStart2 = delta2 - x[2] +1;
			mEnd2 = delta2;
			bilan2 = modelData->getCycleBalanceForGivenPeriod(x[1], mStart2, mEnd2, u[1]);
			//cout<< " bilan 2 = "<<bilan2<<endl;
			res = bilan1 + bilan2;

			break;
		case 1://les deux actifs
			delta1 = modelData->getTimeTransition(x[4], u[0]);
			delta2 = modelData->getTimeTransition(x[4], u[1]);
			//cout<< " delta1 = "<<delta1<< " delta2 "<< delta2<<endl;
			mStart1 = 1;
			mEnd1 = mStart1+min(delta1, delta2)-1;
			//cout<< " mStart1 = "<<mStart1<< " mEnd1 "<< mEnd1<<endl;
			bilan1 = modelData->getCycleBalanceForGivenPeriod(x[0], mStart1, mEnd1, u[0]);

			mStart2 = 1;
			mEnd2 = mStart2+min(delta2, delta2)-1;
			bilan2 = modelData->getCycleBalanceForGivenPeriod(x[1], mStart2, mEnd2, u[1]);

			res = bilan1 + bilan2;
			break;
		case 2://I2 actif
			delta1 = modelData->getTimeTransition(x[4], u[0]);
			delta2 = modelData->getTimeTransition(x[4], u[1]);
			mStart1 = delta1 - x[2] +1;
			mEnd1 = delta1;
			bilan1 = modelData->getCycleMaxDeficitForGivenPeriod(x[0], mStart1, mEnd1, u[0]);

			mStart2 = 1;
			mEnd2 = mStart2+min(delta2, x[2])-1;
			bilan2 = modelData->getCycleMaxDeficitForGivenPeriod(x[1], mStart2, mEnd2, u[1]);

			res = bilan1 + bilan2;

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

	if( (x[4]>=T))
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
		if(x[3]==1)
		{
			if(x[2]>0)
			{
				// cout<< " init constrainte  x="<<x[2]<<" "<<x[3]<<endl;
				res=PLUS_INF;
			}
			else
			{
				res=0.0-deficitMax;
			}
		}
		else
		{
			if(x[2]==0)
			{
				// cout<< " init constrainte  x="<<x[2]<<" "<<x[3]<<endl;
				res=PLUS_INF;
			}
			else
			{
				res=0.0-deficitMax;
			}
		}
	}
	return  res;
}

unsigned long long int maxNbRetro=1;


/*      *****************************************
 *  Definition of constraints and target
 *************************************************** */



/*!
 * \brief Function  defining the mixed  constraints
 *
 * This function defines the set U(x) for admissible controls as function of the state
 * @param x state variable
 * @param u control variable
 * @return  value that caraterise the constraints set
 */
inline double constraintsXU( double * x, double * u )
{
	return 1.0;
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

	int s=(x[2])%12;
	int  saisonalite1, saisonalite2, controlInactifCompatible;
	// x[0]=ibqs1
	//x[1] = Ibqs2
	//x[2] = gamma : temps restant sur la parcelle inactive
	//x[3] = p : parcelle active : 0 si ( I1 actif, I2 passif)  1 si les deux actifs, 2 si (I2 actif et I1 passif)
	//x[4] = tau : temps

	//u[0] = cycle de parcelle 1
	//u[1]=cycle de parcelle 2
	// v : le controle actif
	unsigned long long int delta1, delta2;

	double res;
	switch(x[3])
	{
	case 0://I1 actif => u0 est soumis à la saisonnalité et au déficit, u1 est soumis à la compatibilité de durée avec le temps restant
		delta1 = modelData->getTimeTransition(x[4], u[0]);
		delta2 = modelData->getTimeTransition(x[4], u[1]);
		saisonalite1=modelData->getCycleSeasonality(u[0], s);
		controlInactifCompatible = (delta2 >= x[2])?1:0;


		res=PLUS_INF*( (1.0-saisonalite1) + (1.0 - controlInactifCompatible ) )+1.0*(saisonalite1) * controlInactifCompatible;
		break;
	case 1://les deux actifs
		delta1 = modelData->getTimeTransition(x[4], u[0]);
		delta2 = modelData->getTimeTransition(x[4], u[1]);
		saisonalite1=modelData->getCycleSeasonality(u[0], s);
		saisonalite2=modelData->getCycleSeasonality(u[1], s);

		res=PLUS_INF*( (1.0-saisonalite1) + (1.0 - saisonalite2 ) )+1.0*(saisonalite1) * saisonalite2;

		break;
	case 2://I2 actif
		delta1 = modelData->getTimeTransition(x[4], u[0]);
		delta2 = modelData->getTimeTransition(x[4], u[1]);
		saisonalite2=modelData->getCycleSeasonality(u[1], s);
		controlInactifCompatible = (delta1 >= x[2])?1:0;


		res=PLUS_INF*( (1.0-saisonalite2) + (1.0 - controlInactifCompatible ) )+1.0*(saisonalite2) * controlInactifCompatible;

		break;

	}

	return res;
}
inline double  constraintsXU_fd_DP1( unsigned long long int * x, unsigned long long int * u )
{

	int s=(x[2])%12;
	int saisonalite1, saisonalite2, controlInactifCompatible;
	// x[0]=ibqs1
	//x[1] = Ibqs2
	//x[2] = gamma : temps restant sur la parcelle inactive
	//x[3] = p : parcelle active : 0 si ( I1 actif, I2 passif)  1 si les deux actifs, 2 si (I2 actif et I1 passif)
	//x[4] = tau : temps

	//u[0] = cycle de parcelle 1
	//u[1]=cycle de parcelle 2
	// v : le controle actif
	unsigned long long int delta1, delta2;

	double res;
	switch(x[3])
	{
	case 0://I1 actif => u0 est soumis à la saisonnalité et au déficit, u1 est soumis à la compatibilité de durée avec le temps restant
		delta1 = modelData->getTimeTransition(x[4], u[0]);
		delta2 = modelData->getTimeTransition(x[4], u[1]);
		saisonalite1=modelData->getCycleSeasonality(u[0], s);
		controlInactifCompatible = (delta2 >= x[2])?1:0;


		res=PLUS_INF*( (1.0-saisonalite1) + (1.0 - controlInactifCompatible ) )+1.0*(saisonalite1) * controlInactifCompatible;
		break;
	case 1://les deux actifs
		delta1 = modelData->getTimeTransition(x[4], u[0]);
		delta2 = modelData->getTimeTransition(x[4], u[1]);
		saisonalite1=modelData->getCycleSeasonality(u[0], s);
		saisonalite2=modelData->getCycleSeasonality(u[1], s);

		res=PLUS_INF*( (1.0-saisonalite1) + (1.0 - saisonalite2 ) )+1.0*(saisonalite1) * saisonalite2;

		break;
	case 2://I2 actif
		delta1 = modelData->getTimeTransition(x[4], u[0]);
		delta2 = modelData->getTimeTransition(x[4], u[1]);
		saisonalite2=modelData->getCycleSeasonality(u[1], s);
		controlInactifCompatible = (delta1 >= x[2])?1:0;


		res=PLUS_INF*( (1.0-saisonalite2) + (1.0 - controlInactifCompatible ) )+1.0*(saisonalite2) * controlInactifCompatible;

		break;

	}

	return res;
}
inline double  constraintsXUY_fd( unsigned long long int * x, unsigned long long int * u )
{


	double res;
	int s,  saisonalite1, saisonalite2, controlInactifCompatible;
	// x[0]=ibqs1
	//x[1] = Ibqs2
	//x[2] = gamma : temps restant sur la parcelle inactive
	//x[3] = p : parcelle active : 0 si ( I1 actif, I2 passif)  1 si les deux actifs, 2 si (I2 actif et I1 passif)
	//x[4] = tau : temps

	//u[0] = cycle de parcelle 1
	//u[1]=cycle de parcelle 2
	// v : le controle actif
	unsigned long long int delta1, delta2;
	unsigned long long int mStart1, mEnd1, mStart2, mEnd2;
	double maxDeficit1, maxDeficit2, maxDeficit;
	if(x[4]<T)
	{
		s=(x[4])%12;

		switch(x[3])
		{
		case 0://I1 actif => u0 est soumis à la saisonnalité et au déficit, u1 est soumis à la compatibilité de durée avec le temps restant
			delta1 = modelData->getTimeTransition(x[4], u[0]);
			delta2 = modelData->getTimeTransition(x[4], u[1]);
			saisonalite1=modelData->getCycleSeasonality(u[0], s);
			controlInactifCompatible = (delta2 >= x[2])?1:0;

			mStart1 = 1;
			mEnd1 = mStart1+min(delta1, x[2])-1;
			maxDeficit1 = modelData->getCycleMaxDeficitForGivenPeriod(x[0], mStart1, mEnd1, u[0]);

			mStart2 = delta2 - x[2] +1;
			mEnd2 = delta2;
			maxDeficit2 = modelData->getCycleMaxDeficitForGivenPeriod(x[1], mStart2, mEnd2, u[1]);

			maxDeficit = maxDeficit1 + maxDeficit2;

			res=PLUS_INF*( (1.0-saisonalite1) + (1.0 - controlInactifCompatible ) )+(maxDeficit-deficitMax)*(saisonalite1) * controlInactifCompatible;
			break;
		case 1://les deux actifs
			delta1 = modelData->getTimeTransition(x[4], u[0]);
			delta2 = modelData->getTimeTransition(x[4], u[1]);
			saisonalite1=modelData->getCycleSeasonality(u[0], s);
			saisonalite2=modelData->getCycleSeasonality(u[1], s);
			mStart1 = 1;
			mEnd1 = mStart1+min(delta1, x[2])-1;
			maxDeficit1 = modelData->getCycleMaxDeficitForGivenPeriod(x[0], mStart1, mEnd1, u[0]);

			mStart2 = 1;
			mEnd2 = mStart2+min(delta2, x[2])-1;
			maxDeficit2 = modelData->getCycleMaxDeficitForGivenPeriod(x[1], mStart2, mEnd2, u[1]);

			maxDeficit = maxDeficit1 + maxDeficit2;
			res=PLUS_INF*( (1.0-saisonalite1) + (1.0 - saisonalite2 ) )+(maxDeficit-deficitMax)*(saisonalite1) * saisonalite2;

			break;
		case 2://I2 actif
			delta1 = modelData->getTimeTransition(x[4], u[0]);
			delta2 = modelData->getTimeTransition(x[4], u[1]);
			saisonalite2=modelData->getCycleSeasonality(u[1], s);
			controlInactifCompatible = (delta1 >= x[2])?1:0;


			mStart1 = delta1 - x[2] +1;
			mEnd1 = delta1;
			maxDeficit1 = modelData->getCycleMaxDeficitForGivenPeriod(x[0], mStart1, mEnd1, u[0]);

			mStart2 = 1;
			mEnd2 = mStart2+min(delta2, x[2])-1;
			maxDeficit2 = modelData->getCycleMaxDeficitForGivenPeriod(x[1], mStart2, mEnd2, u[1]);

			maxDeficit = maxDeficit1 + maxDeficit2;
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

void loadModelData();
void loadModelData()
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



	string firstDataFileName = "SimulationGaia.json";


	modelData=new FarmData_Gardener(firstDataFileName);
	deficitMax = modelData->getAdmissibleDeficit();
	T=(double)modelData->getTimeHorizon();
	/*!
	 * \var nbPointsControl[dimc]
	 * number of  discretization points for the control variable
	 * \see controlParams
	 */
	nbPointsControl[0] =  modelData->getNbControls();
	nbPointsControl[1] =  modelData->getNbControls();
	/*!
	 * \var nbPointsState[dim]
	 * number of  discretization points for the state variable
	 * \see gridParams
	 */
	nbPointsState[0]    =  modelData->getNbIPoints();
	nbPointsState[1]    =  modelData->getNbIPoints();
	nbPointsState[2]    =  modelData->getMaxCycleDuration();
	nbPointsState[3]    =  3;
	nbPointsState[4]    =  modelData->getTimeHorizon()+1;
	STATE_MIN[0] = 0.0;
	STATE_MIN[1] = 0.0;
	STATE_MIN[2] = 0.0;
	STATE_MIN[3] = 0.0;
	STATE_MIN[4] = 0.0;

	STATE_MAX[0] = 1.0;
	STATE_MAX[1] = 1.0;
	STATE_MAX[2] = (double) modelData->getMaxCycleDuration();
	STATE_MAX[3] = 3.0;
	STATE_MAX[4] = (double) modelData->getTimeHorizon();

	STATE_MINreel[0] = 0.0;
	STATE_MINreel[1] = 0.0;
	STATE_MINreel[2] = 0.0;
	STATE_MINreel[3] = 0.0;
	STATE_MINreel[4] = 0.0;

	STATE_MAXreel[0] = 1.0;
	STATE_MAXreel[1] = 1.0;
	STATE_MAXreel[2] = (double) modelData->getMaxCycleDuration();
	STATE_MAXreel[3] = 3.0;
	STATE_MAXreel[4] = (double) modelData->getTimeHorizon();

	IStar=modelData->getIStar();
	IRealValues = modelData->getIgrid();
	prefix = modelData->getPrefix();

	cout<< " fini loadModelData T = "<<T<<endl;

	nbTrajs = modelData->getNbTrajs();
	initPoints_fd = new unsigned long long int [dim * nbTrajs];
	initValues_fd = modelData->getInitValues();
	 unsigned long long int *tempinitPoints = modelData->getInitPointsCoords();

	 for(int j = 0; j<nbTrajs; j++)
	 {
		 initPoints_fd[j*dim + 0] = tempinitPoints[2*j + 0];
		 initPoints_fd[j*dim + 1] = tempinitPoints[2*j + 0];
		 initPoints_fd[j*dim + 2] = 0;
		 initPoints_fd[j*dim + 3] = 1;
		 initPoints_fd[j*dim + 4] = tempinitPoints[2*j + 1];

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

	if( (x[4]>=T))
		res1= (x[0]>=IStar && x[1]>=IStar);
	else
		res1=true;
	return (double)res1;
}


inline  double target_fd (unsigned long long int * x)
{


	double res;

	if( (x[4]>=T))
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
	 * Cette fonction est rï¿½servï¿½e ï¿½ l'interpolation du cap
	 * ï¿½ partir du cap rï¿½eel, \phi \in\R.
	 * La fonction interpolle la valeur du cap et le ramï¿½ne ï¿½ [0,2pi]
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



#include "agroEcoDiv_unused.h"

#endif /* DATAGAIA_MODEL1_H_ */
