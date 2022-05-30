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

/*! \var dim
 *  \brief State  dimension
 */
const int dim=2;
/*! \var dimc
 *  \brief Control variable   dimension
 */
const int dimc=1;
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
double CONTROL_MIN[dimc]={0};
/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX[dimc]={1};

/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */
double STATE_MINreel[dim]={0.0,0.0};
double Treel;
/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAXreel[dim]={1.0, 10.0};
double step[dim];


/*
 * le vecteur scaling indique pour chaque variable s'il faut ou pas mettre à l'échelle
 *   1 = mise à l'échelle
 *   0 = pas de mise à l'échelle
 */
const int scaling[dim+1]          = {0,0,0};    //- SCALED VARIABLES (1:scaled, 0:not scaled)



void optimTrajConvert(string fileNameBis, string fileName, int W2, double lScale);

unsigned long long int numParcelle =1;

unsigned long long int nbParcels=1;


/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */

double  STATE_MIN[dim]         = {0,0};
/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */

double  STATE_MAX[dim]         = {1,1};


double  STATE_MIN_V[dim+1]         = {0,0,0};
/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */

double  STATE_MAX_V[dim+1]         = {1.0,1.0,1.0};



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
unsigned long long int nbPointsControl[dimc] =  {1};
/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
unsigned long long int nbPointsState[dim]    =  {10,(unsigned long long int)T+1};
unsigned long long int nbPointsStateV[dim+1]    =  {10,(unsigned long long int)T+1,10};

double surf;
double cFixe;






/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0};
int periodicV[dim+1]={0,0,0};


int compteOnlyOptimalRetro=1;

int computeSet=1;
/*
 * sélection de type de reconstruction de trajectoire
 *  VD= 1, viable par défaut, on sélectionne le premier contrôle viable trouvé
 *  VL= 2, viable lourd: le contrôle reste constant tant qu'il viable;
 *  cette méthode nécessite une initialisation de contrôle
 */
int typeTraj=VD;

double initPoint[dim]={ 1.5, 3.8 };
double initPointV[dim+1]={ 1.5, 3.8, 9.0 };

unsigned long long int initPoint_fd_V[dim+1]={12,0,0};

void triBulle(unsigned long long int * permut, double * tab, unsigned long long int dim);

unsigned long long int projection[dim]={0,0};

unsigned long long int trajProjection[dim]={0,0};

unsigned long long int dirTramage=0;


/*!
 * \var globalDeltaT
 *  bool�en indique si le pas de temps  doit �re choisi globalement
 *  ou localement pour les algorithmes de viabilit�
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

	// x[0]=ibqs
	//x[1] = t
	//u[0] = spec
	//u[1]=pratique
	//u[2]= fertilisant
	double I0=x[0];
	if(x[1]<T)
	{
		image[0]=  modelData->getItransition(x[0], x[1], u[0]) ;
		// cout<< " image 0 = "<<image[0]<<endl;
		image[1]= x[1]+ modelData->getTimeTransition(x[1], u[0]);
		//	 cout<< " image 1 = "<<image[1]<<endl;
	}
	else
	{
		image[0]=x[0];
		image[1]=x[1];
	}
}
inline double l_fd(unsigned long long int  * x, unsigned long long int * u )
{

	double res;
	if(x[1]<T)
	{
		res= modelData->getCycleBalance(x[0], x[1], u[0]);
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

	if( (x[1]>=T))
	{
		if(IRealValues[x[0]]>=IStar)
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
		res=0.0-deficitMax;
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

	int s=(x[1])%12;
	int saisonalite=modelData->getCycleSeasonality(u[0], s);
	double res=PLUS_INF*(1.0-saisonalite)+1.0*(saisonalite);
	return res;
}
inline double  constraintsXU_fd_DP1( unsigned long long int * x, unsigned long long int * u )
{

	int s=(x[2])%12;
	int saisonalite=modelData->getCycleSeasonality(u[0], s);
	double res=PLUS_INF*(1.0-saisonalite)+1.0*(saisonalite);


	return res;
}
inline double  constraintsXUY_fd( unsigned long long int * x, unsigned long long int * u )
{


	double res;
	if(x[1]<T)
	{
		int s=(x[1])%12;
		int saisonalite=modelData->getCycleSeasonality(u[0], s);
		double maxDeficit = modelData->getCycleMaxDeficit(x[0], x[1], u[0]);
		res=PLUS_INF*(1.0-saisonalite)+(maxDeficit-deficitMax)*(saisonalite);
		//if(x[1]<120)
		//res= maxDeficit[u[0]][u[1]][x[0]]-deficitMax;
		//else
		//res= maxDeficit[u[0]][u[1]][x[0]];
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
	/*!
	 * \var nbPointsState[dim]
	 * number of  discretization points for the state variable
	 * \see gridParams
	 */
	nbPointsState[0]    =  modelData->getNbIPoints();
	nbPointsState[1]    =  modelData->getTimeHorizon()+1;
	STATE_MIN[0] = 0.0;
	STATE_MIN[1] = 0.0;
	STATE_MAX[0] = 1.0;
	STATE_MAX[1] = (double) modelData->getTimeHorizon();

	STATE_MINreel[0] = 0.0;
	STATE_MINreel[1] = 0.0;
	STATE_MAXreel[0] = 1.0;
	STATE_MAXreel[1] = (double) modelData->getTimeHorizon();

	IStar=modelData->getIStar();
	IRealValues = modelData->getIgrid();
	prefix = modelData->getPrefix();

	cout<< " fini loadModelData T = "<<T<<endl;

	nbTrajs = modelData->getNbTrajs();
	initPoints_fd = modelData->getInitPointsCoords();
	initValues_fd = modelData->getInitValues();
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

	if( (x[1]>=T))
		res1= (x[0]>=IStar);
	else
		res1=true;
	return (double)res1;
}


inline  double target_fd (unsigned long long int * x)
{


	double res;

	if( (x[1]>=T))
	{
		if(IRealValues[x[0]]>=IStar)
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



#include "agroEcoDiv_unused.h"

#endif /* DATAGAIA_MODEL1_H_ */
