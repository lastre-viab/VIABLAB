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


void optimTrajConvert(string fileNameBis, string fileName, int W2, double lScale);

unsigned long long int numParcelle =1;

unsigned long long int nbParcels=1;
double * surfaces;


double cFixe;

unsigned long long int TT;
void dynamics_fd(unsigned long long int  * x, unsigned long long int *u, unsigned long long int * image)
{

	// x[0]=ibqs
	//x[1] = t
	//u[0] = spec
	//u[1]=pratique
	//u[2]= fertilisant
	double I0=x[0];
	if(x[1]<TT)
	{
		image[0]=  modelData->getItransition(x[0], 0, u[0]) ;
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
	if(x[1]<TT)
	{
		res= modelData->getCycleBalance(x[0], x[1], u[0])*surfaces[0]+modelData->getAnnexCycleBalance(x[1], modelData->getTimeTransition(x[1], u[0]));
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

	if( (x[1]>=TT))
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
	//res = 1.0;
	return res;
}
inline double  constraintsXUY_fd( unsigned long long int * x, unsigned long long int * u )
{


	double res;
	if(x[1]<TT)
	{
		int s=(x[1])%12;
		int saisonalite=modelData->getCycleSeasonality(u[0], s);
		double maxDeficit = modelData->getCycleMaxDeficitWithAnnexe(x[0], x[1], u[0]);
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
	/*!
	 * \var nbPointsState[dim]
	 * number of  discretization points for the state variable
	 * \see gridParams
	 */

	gp->NBPOINTS[0]    =  modelData->getNbIPoints();
	gp->NBPOINTS[1]    =  modelData->getTimeHorizon()+1;


	gp->LIMINF[0] = 0.0;
	gp->LIMINF[1] = 0.0;


	gp->LIMSUP[0] = 1.0;
	gp->LIMSUP[1] = (double) modelData->getTimeHorizon();

gp->FILE_PREFIX = modelData->getPrefix();
avp->FILE_PREFIX = modelData->getPrefix();


	IStar=modelData->getIStar();
	IRealValues = modelData->getIgrid();


	prefix = modelData->getPrefix();

	nbParcels = modelData->getNbParcels();
	surfaces = modelData->getParcelsSurfaces();

	cout<< " fini loadModelData T = "<<sysp->maxTime<<endl;
	int dim = gp->DIM;
	nbTrajs = modelData->getNbTrajs();
	avp->NB_TRAJS = nbTrajs;
	avp->INIT_POINTS_FD = new unsigned long long int [dim * nbTrajs];
	avp->INIT_VALUES_FD = modelData->getInitValues();
	unsigned long long int *tempinitPoints = modelData->getInitPointsCoords();

	for(int j = 0; j<nbTrajs; j++)
	{
		avp->INIT_POINTS_FD[j*dim + 0] = tempinitPoints[2*j + 0];
		avp->INIT_POINTS_FD[j*dim + 1] =  tempinitPoints[2*j + 1];
	}
}
//system("pause");


void postProcess(ParametersManager *PM)
{
	cout<< " post process nettoyage memoire \n";
}


inline  double target_fd (unsigned long long int * x)
{


	double res;

	if( (x[1]>=TT))
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

#endif /* DATAGAIA_MODEL1_H_ */
