#include "../include/ParametersManager.h"
#include "../include/utilities.h"

extern "C" {

std::string paramsFile = "ARKerU_params.json";

//Paramètres des dynamiques

#define NB_POSSIBLE_CROPS 9 // crops available
#define NB_ASSOLEMENTS 126 //9*8*7*6*5/5/4/3/2  // Crops vs possible
#define NB_CROPS 5 // crops planted each year
const float POSSIBLE_CROPS_PARAM[NB_POSSIBLE_CROPS][18] = {2.5,9.2,0.8, 0.2 ,0.12, 0,0,0,175,269,240,548,502,501,445,334,155,0,2, 0.3,0.65, 0.14 ,0.07, 0,0,0,0,0,0,27,10,135,161,13,1,212,1.3, 0.8,0.5, 0.21, 0.15,0,0,0,0,0,0,23,50,22,16,97,0,0,2, 5.3,0.8, 0.28, 0.15 ,0,0,0,178,419,848,918,70,0,0,0,0,0,3, 0.4,0.7, 0.17, 0.099, 0,0,0,0,16,40,4,5,5,115,0,0,0,1.5, 5.4,0.6, 0.09 ,0.01, 0.002,0,0,457,19,55,2075,2044,156,0,0,0,0,2.16, 1.7,0.6, 0.18, 0.12, 0,2,16,185,24,664,0,0,0,1,0,0,0,2, 1.1,0.7, 0.17, 0.1, 0,0,4,11,145,27,17,255,35,4,0,0,0,4 ,5.2,0.8, 0.2, 0.13, 0.002,0,26,143,296,557,705,586,411,235,326,0,0};

double X1MAX = 100000.0; //
double X1MIN = 0.0;
double X2MAX = 1.0;
double X2MIN = 0.0;
double X3MAX = 1.0;
double X3MIN = 0.0;
double U2MAX = 2.0;
double U2MIN = 0.05;
double U3MAX = 15.0;
double U3MIN = 2.0;

//Pour les grilles
int DX1 = 31; //151
int DX2 = 21; //101
int DX3 = 41; //201
int DU2 = 21; //51
int DU3 = 11; //31



//Pour x1
double GAIN = 0.7; // share of the raw income that do not include variable costs (excluding fixed costs)
double FIXED_COSTSR = 5000*12; // includes cookers salaries and equipment depreciation  # initial value = 500*365
double FOOD_PRICE = 1.5; // mean price in € of 1 kg of wholesale food  # initial value = 4
double PRICE_MTX = 2.0; // mean price at which the restaurant buys food from local production

double FIXED_COSTSA = 8000.0 + 2000.0*12;



//Pour x2
double NB_CLIENTS_MAX = 75*180.0; // maximal capacity of the number of meals served each year
double FPP = 0.8 ;// food per plate = the quantity of food in one plate (in kg)




// Pour x3


double THR_WORKLOAD =  4*35.0; //temps de travail max par mois (il y a une personne ici)
double RAP_FALLOW = 0.025;
double WASTE = 0.05;  // proportion of food waste between field and fork (storage issues...)
double R_MID = 0.8; // Production potentielle relative, pour laquelle l'effet sur l'attractivité du restaurant est neutre.



// Paramètres pour les contraintes
double X3QUAL = 0.2;
double X1V = 0.0;

double* XMAX ;
double* XMIN ;
long long unsigned int* DX;
double* UMAX;
double* UMIN;
unsigned long long int *DU;


void loadModelData(const ParametersManager *PM)
{
	const gridParams * gp = PM->getGridParameters();
	XMAX = gp->LIMSUP ;
	X1MAX = XMAX[0];
	X2MAX = XMAX[1];
	X3MAX = XMAX[2];
	XMIN = gp->LIMINF ;
	X1MIN = XMIN[0];
	X2MIN = XMIN[1];
	X3MIN = XMIN[2];
	DX = gp-> NBPOINTS;
	DX1 = DX[0];
	DX2 = DX[1];
	DX3 = DX[2];


	const controlParams *cp = PM->getControlParameters();
	UMAX = cp->LIMSUPC ;
	U2MAX = UMAX[1];
	U3MAX = UMAX[2];
	UMIN = cp->LIMINFC ;
	U2MIN = UMIN[1];
	U3MIN = UMIN[2];
	DU = cp-> NBPOINTSC;
	DU2  = DU[1];
	DU3  = DU[2];


}

float Valeur(int ix,float xmax, float xmin,int Dx){
	float x;
	x=xmin+(xmax-xmin)*(float)(ix)/(Dx-1);
	return x;
}

int CreateAssolements(float Assolements[NB_ASSOLEMENTS][6*NB_CROPS+12]){
	// 1 fonction qui liste les informations nécessaires pour les assolements
	int i,i1,i2,i3,i4,i5,j;
	i = 0;
	//printf("%d\n ",i);

	for (i1=0;i1<NB_POSSIBLE_CROPS-4;i1++){
		for (i2=i1+1;i2<NB_POSSIBLE_CROPS-3;i2++){
			for (i3=i2+1;i3<NB_POSSIBLE_CROPS-2;i3++){
				for (i4=i3+1;i4<NB_POSSIBLE_CROPS-1;i4++){
					for (i5=i4+1;i5<NB_POSSIBLE_CROPS;i5++){
						for (j=0;j<6;j++){
							//printf("%d %d\n ",i,j);
							Assolements[i][5*j] = POSSIBLE_CROPS_PARAM[i1][j];
							Assolements[i][5*j+1] = POSSIBLE_CROPS_PARAM[i2][j];
							Assolements[i][5*j+2] = POSSIBLE_CROPS_PARAM[i3][j];
							Assolements[i][5*j+3] = POSSIBLE_CROPS_PARAM[i4][j];
							Assolements[i][5*j+4] = POSSIBLE_CROPS_PARAM[i5][j];
						}
						for (j=6;j<18;j++){
							Assolements[i][j+24] = (POSSIBLE_CROPS_PARAM[i1][j] + POSSIBLE_CROPS_PARAM[i2][j] + POSSIBLE_CROPS_PARAM[i3][j] + POSSIBLE_CROPS_PARAM[i4][j] + POSSIBLE_CROPS_PARAM[i5][j])/NB_CROPS;
						}
						i = i+1;
					}
				}
			}

		}

	}
	return 0;
}

//Fonctions de calcul de la dynamique de x1

//Fonction qui calcule l'évolution du capital du restaurant

float GainRepas(float x2,float u3){
    return NB_CLIENTS_MAX*x2*u3*GAIN;
}

float CoutAchatFoodExterieur(float x2,float Rn_x3_u1){
	float c;
	float food_quant;
	///Computes expenses linked to food bought from other producers.
	food_quant = NB_CLIENTS_MAX*x2*FPP;
	c = fmaxf(food_quant-Rn_x3_u1,0)*FOOD_PRICE;
	return c;
}


float Cout(int u1,float Assolements[NB_ASSOLEMENTS][6*NB_CROPS+12],float u2){
    int i;
    float expenses=0.0;
    for (i=NB_CROPS;i<2*NB_CROPS;i++){
        expenses = expenses +Assolements[u1][i];
    }
    expenses = expenses*10000 ; // Passage d'€/m² en €/ha
    return expenses* u2 / NB_CROPS + FIXED_COSTSA;
}

float Capital_evolutionA(float Rn_x3_u1,int u1,float Assolements[NB_ASSOLEMENTS][6*NB_CROPS+12],float u2){
    /// Unused !
    float expenses=0.0;
    expenses = Cout(u1,Assolements,u2);
    return Rn_x3_u1*PRICE_MTX-expenses-FIXED_COSTSA;
}

//Fonctions de calcul de la dynamique de x2

// Fonction qui calcule d'évolution du coefficent d'attractivité


float Attractiveness_evolution(float x2,float u3,float Rn_x3_u1){
    float x2next,newx2,u3mid;
    // Variation liée au prix : plus le prix est bas, plus les clients sont contents (avec un effet neutre pour le prix moyen)
    u3mid = U3MIN + (U3MAX-U3MIN)/2;
    newx2 = x2 + (u3mid-u3) / u3mid;
    //Variation liée à l'utilisation des produits de Montrieux : plus une part importante du repas est produite localement, plus les clients sont satisfaits

    newx2 = newx2 + (Rn_x3_u1/(NB_CLIENTS_MAX*x2*FPP)-R_MID); // Il faudrait trouver un moyen d'avoir une idée de la production "moyenne"
    x2next = fmaxf(X2MIN, fminf(X2MAX,newx2));
    return x2next;
}

//Fonctions de calcul de la dynamique de x3


//fonction qui donne la production pour une espece caractérisée par TM et ri et une qualité du sol X3 en kg/m2
float R_n(float x3,float RM,float ri){
        float Rn;
        if (x3==0.0) Rn = 0.0;
        else if (ri!=0.5){
            if (x3<0.5){

                Rn = 2*RM*(x3*(1-2*ri)+ri-sqrtf(ri*ri+2*x3*(1-2*ri)))/(2*ri-1);
            }
            else {
                Rn = RM*(1+2*ri*((1-2*ri)*(x3-0.5)-(1-ri)+sqrtf((1-ri)*(1-ri)+2*(x3-0.5)*(2*ri-1)))/(2*ri-1));

            }

       //Rn = fmax(Rn,0);

       if(Rn<0){
       	printf("Rn négatif (1)\n");
       	printf("Rm x3 ri Rn %f %f %f %f\n",RM,x3,ri,Rn);

       	getchar();
       }
       else if(std::isnan(Rn)){
       	printf("Rn nan (1)\n");
       	printf("Rm ri Rn %f %f %f %f\n",RM,ri,Rn,sqrtf(((1-ri)*(1-ri)+2*(x3-0.5)*(2*ri-1)))/(2*ri-1));

       	getchar();
       }

        }
        else {
            if (x3<0.5){
                Rn = 2*RM*x3;
            }
            else{
                Rn = RM*(x3+0.5);
           }
       if(Rn<0){
       	printf("Rn négatif (2)");
       	printf("Rm ri %f %f\n",RM,ri);
       	getchar();
       }

       }
       return Rn;
}



float R_nS(float x3, int u1,float Assolements[NB_ASSOLEMENTS][6*NB_CROPS+12],float u2){
    ///u2 en hectares, R_nS en kg
    int i,j;
    float rwnum,rwden; // relative_work numerateur / dénominateur
    float RnS=0;
    for (i=0;i<NB_CROPS;i++){
            RnS = RnS+R_n(x3,Assolements[u1][i],Assolements[u1][i+2*NB_CROPS]);
    }
    RnS = RnS*10000;//passage m2 à hectare
    RnS = RnS*u2/NB_CROPS;
    rwnum = 0;
    rwden = 0;
    for (j=0;j<12;j++){
            rwnum = rwnum+fminf(THR_WORKLOAD,Assolements[u1][6*NB_CROPS+j]*u2);
            rwden = rwden+Assolements[u1][6*NB_CROPS+j]*u2;
    }
    if (rwden>0) {
	    RnS = RnS*rwnum/rwden;
    }

    return RnS*(1-WASTE);
}

//fonction qui calcule une borne supérieure de la production maximale de l activité agricole

float production_maximale(){
	float Assolements[NB_ASSOLEMENTS][6*NB_CROPS+12];
	float rcurrent,rres;
	int i,j,b;
	b = CreateAssolements(Assolements);
	rres = 0;
	for(i=0;i<NB_ASSOLEMENTS;i++) {
		for(j=0;j<DU2;j++) {
			rcurrent = R_nS(X3MAX,i,Assolements,Valeur(j,U2MAX,U2MIN,DU2));
			//printf("i %d j %d u2 %f %f\n",i,j,Valeur(j,U2MAX,U2MIN,DU2),rcurrent);
			if (rres < rcurrent) rres = rcurrent;
		}
	}
	return rres;
}

float cout_maximal(){
	float Assolements[NB_ASSOLEMENTS][6*NB_CROPS+12];
	float ecurrent,eres;
	int i,j,b;
	b = CreateAssolements(Assolements);
	eres = 0;
	for(i=0;i<NB_ASSOLEMENTS;i++) {
		for(j=0;j<DU2;j++) {
			ecurrent = Cout(i,Assolements,Valeur(j,U2MAX,U2MIN,DU2));
			if (eres < ecurrent) eres = ecurrent;
		}
	}
	return eres;
}


float gain_maximal(float Rmax){
	float ecurrent,eres,x2min,x2max,x2,u3,u3min,u3max;
	int i,j;
	eres = 0;

	x2min = X2MIN;
	x2max= X2MAX;
    u3min = U3MIN;
	u3max = U3MAX;

	for(i=0;i<DX2;i++) {
            x2 = Valeur(i,x2max,x2min,DX2);
		for(j=0;j<DU3;j++) {
            u3 = Valeur(j,u3max,u3min,DX3);
			ecurrent = GainRepas(x2,u3)-CoutAchatFoodExterieur(x2,Rmax)-FIXED_COSTSR;
			if (eres < ecurrent) eres = ecurrent;
		}
	}
	return eres;
}


float BISQ_evol(float x3, int u1,float Assolements[NB_ASSOLEMENTS][6*NB_CROPS+12],float u2,float x3min,float x3max){
	float dIb = 0;
	float dIp = 0;
	float yn;
	int i;
	for (i=0;i<NB_CROPS;i++){
        yn = R_n(x3,Assolements[u1][i],Assolements[u1][2*NB_CROPS+i]);
		dIb = dIb + yn/(2*Assolements[u1][i])*Assolements[u1][3*NB_CROPS+i]*u2 / NB_CROPS;
		dIp = dIp + (x3 * Assolements[u1][4*NB_CROPS+i] + Assolements[u1][5*NB_CROPS+i])*u2/NB_CROPS;
	}
	return fminf(fmaxf(x3 - dIb + dIp+RAP_FALLOW*(U2MAX-u2),x3min),x3max);
}



void dynamics(const double *x, const double *u, double *image)
{ 
	// Importance que cette fonction s'appelle dynamics !
	// x[0] = x1, x[1] = x2, x[2] = x3

	//printf("x3 %f", x[2]);
	float Assolements[NB_ASSOLEMENTS][6*NB_CROPS+12];
	CreateAssolements(Assolements);
	float Rn_x3_u1 = R_nS(x[2], u[0],Assolements,u[1]);
	
	image[0] = x[0] + GainRepas(x[1],u[2]) - CoutAchatFoodExterieur(x[1],Rn_x3_u1)-FIXED_COSTSR-Cout(u[0],Assolements,u[1]) ;
	image[1] = Attractiveness_evolution(x[1],u[2],Rn_x3_u1);
	image[2] = BISQ_evol(x[2],u[0],Assolements,u[1],X3MIN,X3MAX);
	
}



int IsinConstraintSet3D(double x1, double x2, double x3){
    // Renvoie 1 si viable, 0 sinon
	int b=1;
	if (x1<0) b = 0;
	else if (x3<X3QUAL) b = 0;
	return b;
}


 double constraintsX( const double *x )
{
	return (IsinConstraintSet3D(x[0], x[1], x[2])) ? 1.0 : PLUS_INF;
}

}
