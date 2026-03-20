#include "../include/utilities.h"
#include "../include/ParametersManager.h"
#include "../include/defs.h"
#include <cstdlib>

extern "C"
{

std::string paramsFile = "VALIUM_hybrid.json";

//------------------------------------------------------------------------------------------------------ 
//  Model description source for the problem VALIUM
//------------------------------------------------------------------------------------------------------ 

//Paramètres des dynamiques

#define NB_CULTURES 6    // Nombre de spéculations
#define NB_PARTIES 1     // Nombre de parcelles
#define NV_PUITS 9       // Niveau de puits
#define MT_IRRI 3        // Méthod d'irrigation 

// Table 2 : prix du surcreusement
const double SURCREUSE_PARAM[NV_PUITS] =
	{
		0, 2, 4, 10, 20, 32, 45, 60, 77
	};

// Table 3 : besoins en eau pour chaque spéculation M,L,G
const double BES_PARAM[NB_CULTURES + 1][MT_IRRI] =
	{
		0, 0, 0,
		11.0, 13.0, 8.0,
		15.0, 18.0, 11.0,
		10.0, 12.0, 7.0,
		 5.0, 7.0, 4.0,
		 6.0, 8.0, 5.0,
		5.0, 6.0, 4.0
	};

// Table 4 : coût de production mécanique M,L,G et prix de vente
const double CP_PV_PARAM[NB_CULTURES + 1][MT_IRRI + 1] =
	{
		0, 0, 0, 0,
		0, 10, 10, 25,   // 25,15  8.3
		0,
		8,
		8,
		15,     // 15,7
		0,
		6,
		6,
		12,     // 12,6
		0,
		5,
		5,
		9,      // 9,4
		0,
		6,
		6,
		10,     // 10,4
		0,
		4,
		4,
		6
		// 6,2
	};

double gamma_para;
double FIXED_COSTSR;
double FIXED_PLUIE;
double BES_para;
double *CC_LimInf;
double *CC_LimSup;
double *CC_Step;
int CC_dim;
unsigned long long int *CC_nbPoints;
unsigned long long int * continuousControlIntCoords;
void loadModelData(const ParametersManager *PM)
    {
    spdlog::info("[LoadModelData] : STARt");
    const modelParams *modelParams = PM->getModelParameters();
    gamma_para = modelParams->getDouble("GAMMA_PARAM");
    FIXED_COSTSR = modelParams->getDouble("COST_PARAM");
    FIXED_PLUIE = modelParams->getDouble("PLUIE_PARAM");
    BES_para = modelParams->getDouble("BES_PARAM");
    spdlog::info("[LoadModelData] : model params ok");
    const controlParams *cp = PM->getControlParameters();
    CC_dim = cp->DIMC;
    CC_nbPoints = cp->NBPOINTSC;
    CC_LimInf = cp->LIMINFC;
    CC_LimSup = cp->LIMSUPC;
    CC_Step = new double[CC_dim];
    for (int dc = 0; dc < CC_dim; dc++)
	{
	CC_Step[dc] = (CC_LimSup[dc] - CC_LimInf[dc]) / (CC_nbPoints[dc] - 1);
	}
    continuousControlIntCoords = new unsigned long long int [CC_dim];
    spdlog::info("[LoadModelData] : Finished");
    }

inline void getIntControlCoords(const double *coords, unsigned long long int *ud)
    {
    for (int i = 0; i < CC_dim; i++)
	{
	ud[i] = (unsigned long long int) std::floor(((coords)[i] - CC_LimInf[i]) / CC_Step[i]);
	}
    }

// Fonctions de calcul de la dynamique de x1=C

// CS(x2,u1,T2)=CS(P/5,S) Coût Surcreusement
double CoutSurcreuse(const unsigned long long int currentLevel, const unsigned long long int targetLevel)
    {
    double CS = 0.0;
    for (int i = currentLevel + 1; i < targetLevel + 1; i++)
	{
	CS += SURCREUSE_PARAM[i];
	}
    return CS;
    }

// CI(u2)=CI(I+) Coût Investissement
double CoutInvest(const unsigned long long int method_add)
    {
    double CI = 0.0;
    if (method_add == 1)
	{
	CI = 5.0;    // + Lance
	}
    else if (method_add == 2)
	{
	CI = 23.0;   // + Goutte à goutte
	}
    else if (method_add == 0)
	{
	CI = 0.0;
	}
    return CI;
    }

// Stocker le nombre de parties où chaque spéculation est plantée
int StockerNBCulture(const double *Sp, int NBCulture[NB_CULTURES])
    {
    for (int i = 0; i < NB_PARTIES; i++)
	{
	int Sp_i = (int) Sp[i];
	if (Sp_i != 0)
	    {
	    NBCulture[Sp_i - 1] += 1;
	    }
	}
    return 0;
    }

// Calculer le ombre de parties non vide
int NonEmptyParcelles(const double *Sp)
    {
    int N = 0;
    for (int i = 0; i < NB_PARTIES; i++)
	{
	int Sp_i = (int) Sp[i];
	if (Sp_i != 0)
	    {
	    N++;
	    }
	}
    return N;
    }

// PV(u678,A,N,T4)=PV(Sp) Prix Vente
double PrixVente(const unsigned long long int * ud)
    {
    return CP_PV_PARAM[ud[0]][3];
    }


// CP(u345,u678,T4)=CP(Sp,M) Coût de Production
double CoutProduction(const unsigned long long int *ud)
    {
    return CP_PV_PARAM[ud[0]][ud[1]];
    }

// Fonctions de calcul de la dynamique de x2=P

// P=5*S
double ProfondPuits_Evolution(const double level_surc)
    {
    return 5.0 * level_surc;
    }

//Fonctions de calcul de la dynamique de x3=E

// E+=I+
double Equipement_Evolution(const double method_add)
    {
    return method_add;
    }

//Fonctions de calcul de la dynamique de x4=h

// BES(u345,u678,T3)=BES(Sp,M) Besoin en Eau de la Spéculation
double BesoinEnEau(const unsigned long long int * ud)
    {
    return BES_PARAM[ud[0]][ud[1]];
    }

void dynamics_hybrid_c(const double *x, const unsigned long long int *xd, const double *u, double *image)
    {
    //u[0] = spec
    //u[1] = irridg
    //x[0] : capital
    //x[1] : niveau de nappe
    getIntControlCoords(u, continuousControlIntCoords);

    double PV_f = PrixVente(continuousControlIntCoords);

    double CP_f = CoutProduction(continuousControlIntCoords);

    double BES_f = BesoinEnEau(continuousControlIntCoords);

    image[0] = gamma_para * x[0] - CP_f + PV_f - FIXED_COSTSR;

    image[1] = x[1] + BES_f - FIXED_PLUIE;
    }

void dynamics_hybrid_d(const double *xc, const unsigned long long int *xd, const unsigned long long int *ud, unsigned long long int *image)
    {
    //ud[0] = surcreusement
    //ud[1] = investissement
    //x[0] : profondeur puit
    //x[1] : equipement

    image[0] = ud[0]; // nouveau niveau
    image[1] = xd[1] + ud[1]; // n
    }


	void resetMap_hybrid(const double *x, const unsigned long long int *currentDiscreteState, const unsigned long long int *discreteControl,
		double *image, unsigned long long int *nextDiscreteState)
    {
    // discrete control : [0] = surcreusement
    // discrete control : [1] = investissement
    // x[0] = capital
    // x[1] = niveau de la nappe
    // discreteState [0] = niveau puit
    // DiscreteState [1] = equipement
    // Cette méthode détermine le saut de l'état continue qui accompagne l'évolution de l'état discret


    // double CP_c = CoutProduction(M_c, Sp_c, CP_PV_PARAM);

    // Importance que cette fonction s'appelle dynamics !
    // x[0] = x1 = C(t), x[1] = x2 = P(t), x[2] = x3 = E(t), x[3] = x4 = h(t)

    image[0] = x[0] - CoutSurcreuse(currentDiscreteState[0], discreteControl[0]) - CoutInvest(discreteControl[1]);
    image[1] = x[1];
	nextDiscreteState[0] = currentDiscreteState[0];
	nextDiscreteState[1] = currentDiscreteState[1];
    }

// Constraints

// Check M[](t) constraint: M[](t) in {0,E(t),I+(t)}
bool CheckMConstraints(const unsigned long long int M, const unsigned long long int  E, const unsigned long long int Iplus)
    {

    if (E == 3)
	{
	return true;
	}

    if (M != 0 && M != E && M != Iplus)
	{
	return false;
	}
    return true;
    }

// Check I+(t) constraint
bool CheckIConstraint(const unsigned long long int  E, const unsigned long long int Iplus)
    {
    if (Iplus < 0)
	{
	return false;
	}
    // si E(t)=0，I+(t)=0,1,2
    if (E == 0)
	{
	return (Iplus == 0 || Iplus == 1 || Iplus == 2);
	}
    // si E(t)=1, I+(t)=0,2
    else if (E == 1)
	{
	return (Iplus == 0 || Iplus == 2);
	}
    // si E(t)=2, I+(t)=0,1
    else if (E == 2)
	{
	return (Iplus == 0 || Iplus == 1);
	}
    // si E(t)=3, I+(t)=0
    else if (E == 3)
	{
	return (Iplus == 0);
	}
    else
	{
	return false;
	}
    }


double constraintsXU_hybrid(const double *x, const unsigned long long int *xd, const double *u, const unsigned long long int *ud)
    {
    // discrete control : ud[0] = surcreusement
    // discrete control : ud[1] = investissement
    // x[0] = capital
    // x[1] = niveau de la nappe
    // discreteState xd[0] = niveau puit
    // DiscreteState xd[1] = equipement
    // cont control M : mode irrigation = u[1]
    // cont control S : spec = u[0]

    getIntControlCoords(u, continuousControlIntCoords);

    if (!CheckMConstraints(continuousControlIntCoords[1], xd[1], ud[1]))
	{
	//spdlog::info("[contraints XU ] : M pas OK");
	return PLUS_INF;
	}

    // Check S(t) constraint: S(t) in [P(t)/5, 8]
    if (xd[0] > ud[0])
	{
	return PLUS_INF;
	}

    // Check I+(t) constraint
    if (!CheckIConstraint(xd[1], ud[1]))
	{
	return PLUS_INF;
	}

    // Check the constraint: h(t)+BES_f(t)*a+BES_c(t)*a<=S(t)*5
    double BES_f = BesoinEnEau(continuousControlIntCoords);
    // double BES_c = BesoinEnEau(M_c, Sp_c, BES_PARAM);

    if (x[1] + BES_f * BES_para > ud[0] * 5)
	{
	// spdlog::info("[contraints XU ] :niveau de la nappe {}", x[1]);
	// spdlog::info("[contraints XU ] :besoin eau {}", BES_f * BES_para);
	// spdlog::info("[contraints XU ] :niveau puit {}", ud[0] * 5);
	return PLUS_INF;
	}

    // Check the constraint: h>=30 and M_k=0, then Sp_k=0

    if (x[1] >= 30 && continuousControlIntCoords[1] == 0 && continuousControlIntCoords[0] != 0)
	{
	// spdlog::info("[contraints XU ] : Pas droit d'arroser à la main");
	return PLUS_INF;
	}


    // Check the constraint: C(t)

    double CS = CoutSurcreuse(xd[0], ud[0]);
    double CI = CoutInvest(ud[1]);
    double PV_f = PrixVente(continuousControlIntCoords);

    double CP_f = CoutProduction(continuousControlIntCoords);

    // Constraint 1: γ*C(t)-CS(t)-CI(t)-CP_f(t)-CP_c(t)-FC>=0
    if (gamma_para * (x[0] - CS - CI) - CP_f - FIXED_COSTSR < 0)
	{
//	 spdlog::info("[contraints XU ] : pas assez de capital");
//	 spdlog::info("[contraints XU ] : capital {}", x[0]);
//	 spdlog::info("[contraints XU ] : cout creuser {}", CS);
//	 spdlog::info("[contraints XU ] : cout investir {}", CI);
//	 spdlog::info("[contraints XU ] : cout produire {}", CP_f);
//	 spdlog::info("[contraints XU ] : cout fixe {}", FIXED_COSTSR);
	return PLUS_INF;
	}
    if (gamma_para * x[0] < 0)
	{
	return PLUS_INF;
	}

    // If all constraints are satisfied, return 1.0 to indicate feasibility
   // spdlog::info("[contraints XU ] : OK");
    return 1.0;
    }
}
