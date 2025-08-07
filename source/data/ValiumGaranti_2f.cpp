#include "../include/utilities.h"
#include "../include/ParametersManager.h"

#include <cstdlib>

extern "C" {

std::string paramsFile = "ValiumGaranti_2f.json";

//------------------------------------------------------------------------------------------------------ 
//  Model description source for the problem VALIUM
//------------------------------------------------------------------------------------------------------ 

//Paramètres des dynamiques

#define NB_SAISONS 1    //Nombre de saisons
#define NB_CULTURES 6    // Nombre de spéculations
#define NB_PLAYERS 4
#define NB_PARTIES_P1 3     // Nombre de parcelles de P1
#define NB_PARTIES_P2 1     // Nombre de parcelles de P2-->Moi
#define NB_PARTIES_P3 2     // Nombre de parcelles de P3
#define NB_PARTIES_P4 3     // Nombre de parcelles de P4
#define NB_PARTIES_PLAYER (NB_PARTIES_P2)   // Système pour player N° 2
#define NB_PARTIES_TOTAL (NB_PARTIES_P1+NB_PARTIES_P2+NB_PARTIES_P3+NB_PARTIES_P4)
#define NB_PARTIES_OTHERS (NB_PARTIES_TOTAL-NB_PARTIES_PLAYER)
#define NV_PUITS 9       // Niveau de puits
#define MT_IRRI 3        // Méthod d'irrigation 

// Table 2 : prix du surcreusement
const float SURCREUSE_PARAM[NV_PUITS] = {
0,2,4,10,20,32,45,60,77
}; 

// Table 3 : besoins en eau pour chaque spéculation M,L,G
const float BES_PARAM[NB_CULTURES+1][MT_IRRI] = {
0,0,0,
11,13,8,
15,18,11,  
10,12,7,
5,7,4,
6,8,5,
5,6,4
}; 

// Table 4 : coût de production mécanique M,L,G et prix de vente
const float CP_PV_PARAM[NB_CULTURES+1][MT_IRRI+1] = {
0,0,0,0,
0,10,10,25,
0,8,8,15,
0,6,6,12,
0,5,5,9,
0,6,6,10,
0,4,4,6
}; 

double gamma_para; 
double FIXED_COSTSR; 
double FIXED_PLUIE;
double BES_para;
void loadModelData(const ParametersManager *PM)
{
    const modelParams *modelParams = PM->getModelParameters();
    gamma_para = modelParams->getDouble("GAMMA_PARAM");
    FIXED_COSTSR = modelParams->getDouble("COST_PARAM");
    FIXED_PLUIE = modelParams->getDouble("PLUIE_PARAM");
    BES_para = modelParams->getDouble("BES_PARAM");
}


// Fonctions de calcul de la dynamique de x1=C

// CS(x2,u1,T2)=CS(P/5,S) Coût Surcreusement
float CoutSurcreuse(const double profon_puis,const double level_surc,const float SURCREUSE_PARAM[NV_PUITS]){
	float CS = 0.0;
	 for(int i=(int)(profon_puis/5)+1;i<(int)(level_surc)+1;i++){
		 CS += SURCREUSE_PARAM[i];
	 }
	return CS;
}


// CI(u2)=CI(I+) Coût Investissement
float CoutInvest(const double Iplus){
	float CI = 0.0;
	if(Iplus==1){
		CI = 5.0;    // + Lance
	}
	else if(Iplus==2){
		CI = 23.0;   // + Goutte à goutte
    }
	return CI;
}


// Stocker le nombre de parties où chaque spéculation est plantée
int StockerNBCulture(const double * Sp, int NBCulture[NB_CULTURES]){
    for(int i=0;i<NB_PARTIES_TOTAL;i++){
        int Sp_i = (int)Sp[i];
        if (Sp_i!=0){
            NBCulture[Sp_i-1] += 1;
        }
    }
	return 0;
}

// Calculer le nombre de parties non vide
int NonEmptyParcelles(const double * Sp){
	int N=0;
	for(int i=0;i<NB_PARTIES_TOTAL;i++){
	    int Sp_i = (int)Sp[i];
        if (Sp_i!=0){
           N++;
        }
    }
    return N;
}


// PV(u678,A,N,T4)=PV(Sp) Prix Vente
float PrixVente(const double * Sp,int NB_PARTS,int NBCulture[NB_CULTURES],int NonEmptyParts,const float CP_PV_PARAM[NB_CULTURES+1][MT_IRRI+1]){
	float PV = 0.0;
	
	if(NonEmptyParts==0){
        return 0.0f;
    }

	for(int i=0;i<NB_PARTS;i++){
      int Sp_i = (int)Sp[i];
      if(Sp_i!=0){
        float ratio = (float) NBCulture[Sp_i-1] / NonEmptyParts;
        // initial 1/3
        if (ratio<=1){
          PV += CP_PV_PARAM[Sp_i][3];
        }
        else{
          PV += CP_PV_PARAM[Sp_i][3] / 3;
        }
      }
    }	
	return PV;
}


// CP(u345,u678,T4)=CP(Sp,M) Coût de Production
float CoutProduction(const double M,const double * Sp,int NB_PARTS,const float CP_PV_PARAM[NB_CULTURES+1][MT_IRRI+1]){
	float CP = 0.0;
  int M_int = (int)M;
	for(int i=0;i<NB_PARTS;i++){
	    int Sp_i = (int)Sp[i];
		CP += CP_PV_PARAM[Sp_i][M_int];
	}
	return CP;
}


// Fonctions de calcul de la dynamique de x2=P

// P=5*S
float ProfondPuits_Evolution(const double level_surc){
    return 5.0*level_surc;
}


//Fonctions de calcul de la dynamique de x3=E

// E+=I+
float Equipement_Evolution(const double Iplus){
    return Iplus;
}


//Fonctions de calcul de la dynamique de x4=h

// BES(u345,u678,T3)=BES(Sp,M) Besoin en Eau de la Spéculation
float BesoinEnEau(const double M,const double * Sp,int NB_PARTS,const float BES_PARAM[NB_CULTURES+1][MT_IRRI]){
	float BES = 0.0;
  int M_int = (int)M;
	for(int i=0;i<NB_PARTS;i++){
	    int Sp_i = (int)Sp[i];
		BES += BES_PARAM[Sp_i][M_int];
	}
	return BES;
}

// void (*dynamics_tych)(const double*, const double*, const double*, double*);
// Système dynamique pour joueur N°2
void dynamics_tych(const double * x, const double * u, const double * v, double * image)
{
    const double* Sp_player_f = &u[2];
    const double* Sp1_f = &v[NB_PLAYERS-1];
    const double* Sp2_f = &v[NB_PLAYERS-1 +NB_PARTIES_P1];
    const double* Sp3_f = &v[NB_PLAYERS-1 +NB_PARTIES_P1+NB_PARTIES_P2];
    const double* Sp4_f = &v[NB_PLAYERS-1 +NB_PARTIES_P1+NB_PARTIES_P2+NB_PARTIES_P3];
    const double* Sp_total_f = &v[NB_PLAYERS-1];
    int NBCulture_f[NB_CULTURES] = {0};
    StockerNBCulture(Sp_total_f,NBCulture_f);
    int NonEmptyParts_f = NonEmptyParcelles(Sp_total_f);

    float PV_player_f = PrixVente(Sp_player_f, NB_PARTIES_PLAYER, NBCulture_f, NonEmptyParts_f, CP_PV_PARAM);  
    float CP_player_f = CoutProduction(u[2], Sp_player_f, NB_PARTIES_PLAYER, CP_PV_PARAM);  
    float BES_player_f = BesoinEnEau(u[2], Sp_player_f, NB_PARTIES_PLAYER, BES_PARAM); 
    float BES1_f = BesoinEnEau(v[0], Sp1_f, NB_PARTIES_P1, BES_PARAM); 
    float BES3_f = BesoinEnEau(v[1], Sp3_f, NB_PARTIES_P3, BES_PARAM); 
    float BES4_f = BesoinEnEau(v[2], Sp4_f, NB_PARTIES_P4, BES_PARAM); 
    
	// Importance que cette fonction s'appelle dynamics !
	// x[0] = x1 = C1(t), x[1] = x2 = P1(t), x[2] = x3 = E1(t), x[3] = x7 = h(t)

	// Etat P1
	image[0] = gamma_para*x[0] - CoutSurcreuse(x[1],u[0],SURCREUSE_PARAM) - CoutInvest(u[1]) - CP_player_f + PV_player_f - FIXED_COSTSR;
	image[1] = ProfondPuits_Evolution(u[0]);
	image[2] = x[2] + Equipement_Evolution(u[1]);
	// Etat global
	image[3] = x[3] + BES_player_f + BES1_f + BES3_f + BES4_f - FIXED_PLUIE;
}


// Constraints

// Check M[](t) constraint: M[](t) in {0,E(t),I+(t)}
bool CheckMConstraints(const double M, const double E, const double Iplus) {
  if (E == 3){
    return true;
  }
  if(!(M == 0 || M == E || M == Iplus)){
      return false;
  }
  return true;
}

// Check S(t) constraint: S(t) in [P(t)/5, 8]
bool CheckSConstraint(const double S, const double P){
  return (S>=P/5) && (S<=8);
}

// Check I+(t) constraint
bool CheckIConstraint(const double E, const double Iplus) {
  if (Iplus < 0) {
    return false; 
  }
  // si E(t)=0，I+(t)=0,1,2
  if (E == 0) {
    return (Iplus == 0 || Iplus == 1 || Iplus == 2);
  }
  // si E(t)=1, I+(t)=0,2
  else if (E == 1) {
    return (Iplus == 0 || Iplus == 2);
  }
  // si E(t)=2, I+(t)=0,1
  else if (E == 2) {
    return (Iplus == 0 || Iplus == 1);
  }
  // si E(t)=3, I+(t)=0
  else if (E == 3) {
    return (Iplus == 0);
  }
  else {
    return false;
  }
}

//Constraints pour joueur N°2
double constraintsXU_tych(const double * x, const double * u, const double * v)
{
    const double* Sp_player_f = &u[2];
    const double* Sp1_f = &v[NB_PLAYERS-1];
    const double* Sp2_f = &v[NB_PLAYERS-1 +NB_PARTIES_P1];
    const double* Sp3_f = &v[NB_PLAYERS-1 +NB_PARTIES_P1+NB_PARTIES_P2];
    const double* Sp4_f = &v[NB_PLAYERS-1 +NB_PARTIES_P1+NB_PARTIES_P2+NB_PARTIES_P3];
    const double* Sp_total_f = &v[NB_PLAYERS-1];
    int NBCulture_f[NB_CULTURES] = {0};
    StockerNBCulture(Sp_total_f,NBCulture_f);
    int NonEmptyParts_f = NonEmptyParcelles(Sp_total_f);
    
    // Constraint for tych: &v[6]=&u[2]
    for(int i=0;i<NB_PARTIES_PLAYER;i++){
      if(Sp2_f[i] != Sp_player_f[i]){
        return PLUS_INF;
      }
    }

    // Check M[](t) constraint
    if(!CheckMConstraints(u[2], x[2], u[1])){
      return PLUS_INF;
    }

    // Check S(t) constraint
    if(!CheckSConstraint(u[0], x[1])){
      return PLUS_INF;
    }
  
	  // Check I+(t) constraint
	  if (!CheckIConstraint(x[2], u[1])) {
	    return PLUS_INF;
	  }
    
    float BES_player_f = BesoinEnEau(u[2], Sp_player_f, NB_PARTIES_PLAYER, BES_PARAM); 
	  float BES1_f = BesoinEnEau(v[0], Sp1_f, NB_PARTIES_P1, BES_PARAM); 
	  float BES3_f = BesoinEnEau(v[1], Sp3_f, NB_PARTIES_P3, BES_PARAM); 
	  float BES4_f = BesoinEnEau(v[2], Sp4_f, NB_PARTIES_P4, BES_PARAM); 
	  float BES_others_f = BES1_f + BES3_f + BES4_f;
	  
	  // Check the constraint: h(t)+BES_f(t)*a<=S(t)*5
    if(x[3]+BES_player_f*BES_para+BES_others_f*BES_para>u[2]*5){
      return PLUS_INF;
    }
    
    // Check the constraint: h>=30 and M=0, then Sp=0
    for(int i=0;i<NB_PARTIES_PLAYER;i++){
      if(x[3]+BES_others_f*BES_para>=30 && u[2]==0 && Sp_player_f[i]!=0){
        return PLUS_INF;
      }
    }
    //P1
    for (int i=0;i<NB_PARTIES_P1;i++){
      if(x[3]>=30 && v[0]==0 && Sp1_f[i]!=0){
        return PLUS_INF;
      }
    }
    //P3
    for (int i=0;i<NB_PARTIES_P3;i++){
      if(x[3]+BES1_f*BES_para>=30 && v[1]==0 && Sp3_f[i]!=0){
        return PLUS_INF;
      }
    }
    //P4
    for (int i=0;i<NB_PARTIES_P4;i++){
      if(x[3]+BES1_f*BES_para+BES3_f*BES_para>=30 && v[2]==0 && Sp4_f[i]!=0){
        return PLUS_INF;
      }
    }
    	
    // Check the constraint C(t): γ*C(t)-CS(t)-CI(t)-CP_f(t)-CP_c(t)-FC>=0
    float CS_player = CoutSurcreuse(x[1],u[0],SURCREUSE_PARAM);  
    float CI_player = CoutInvest(u[1]);                       
    float CP_player_f = CoutProduction(u[2], Sp_player_f, NB_PARTIES_PLAYER, CP_PV_PARAM);    
    if(gamma_para*x[0]-CS_player-CI_player-CP_player_f-FIXED_COSTSR<0){
  	    return PLUS_INF;
    }

    // If all constraints are satisfied, return 1.0 to indicate feasibility
    return 1.0;
}

}
