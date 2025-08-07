#include "../include/utilities.h"
#include "../include/ParametersManager.h"

#include <cstdlib>

extern "C" {

std::string paramsFile = "VALIUM.json";

//------------------------------------------------------------------------------------------------------ 
//  Model description source for the problem VALIUM
//------------------------------------------------------------------------------------------------------ 

//Paramètres des dynamiques

#define NB_CULTURES 6    // Nombre de spéculations
#define NB_PARTIES 1     // Nombre de parcelles
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
0,10,10,25,   // 25,15  8.3
0,8,8,15,     // 15,7
0,6,6,12,     // 12,6
0,5,5,9,      // 9,4
0,6,6,10,     // 10,4
0,4,4,6       // 6,2
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
float CoutInvest(const double method_add){
  float CI = 0.0;
  if(method_add==1){
    CI = 5.0;    // + Lance
  }
  else if(method_add==2){
    CI = 23.0;   // + Goutte à goutte
  }
  else if(method_add==0){
    CI = 0.0;
  }
  return CI;
}


// Stocker le nombre de parties où chaque spéculation est plantée
int StockerNBCulture(const double * Sp,int NBCulture[NB_CULTURES]){
    for(int i=0;i<NB_PARTIES;i++){
        int Sp_i = (int)Sp[i];
        if (Sp_i!=0){
            NBCulture[Sp_i-1] += 1;
        }
    }
	return 0;
}

// Calculer le ombre de parties non vide
int NonEmptyParcelles(const double * Sp){
	int N=0;
	for(int i=0;i<NB_PARTIES;i++){
	    int Sp_i = (int)Sp[i];
        if (Sp_i!=0){
           N++;
        }
    }
    return N;
}


// PV(u678,A,N,T4)=PV(Sp) Prix Vente
float PrixVente(const double * Sp,int NBCulture[NB_CULTURES],int NonEmptyParts,const float CP_PV_PARAM[NB_CULTURES+1][MT_IRRI+1]){
	float PV = 0.0;
	
	if(NonEmptyParts==0){
        return 0.0f;
    }

    for(int i=0;i<NB_PARTIES;i++){
        int Sp_i = (int)Sp[i];
        if(Sp_i!=0){
            PV += CP_PV_PARAM[Sp_i][3];
        }
    }

	// for(int i=0;i<NB_PARTIES;i++){
    //     int Sp_i = (int)Sp[i];
    //     if(Sp_i!=0){
    //         float ratio = (float) NBCulture[Sp_i-1] / NonEmptyParts;
    //         // initial 1/3
    //         if (ratio<=1){
    //             PV += CP_PV_PARAM[Sp_i][3];
    //         }
    // 	    else{
    //    	        PV += CP_PV_PARAM[Sp_i][3] / 3;
    //         }
    //     }
    // }	
	return PV;
}


// CP(u345,u678,T4)=CP(Sp,M) Coût de Production
float CoutProduction(const double * M,const double * Sp,const float CP_PV_PARAM[NB_CULTURES+1][MT_IRRI+1]){
	float CP = 0.0;
	for(int i=0;i<NB_PARTIES;i++){
	    int Sp_i = (int)Sp[i];
	    int M_i = (int)M[i];
		CP += CP_PV_PARAM[Sp_i][M_i];
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
float Equipement_Evolution(const double method_add){
  return method_add;
}


//Fonctions de calcul de la dynamique de x4=h

// BES(u345,u678,T3)=BES(Sp,M) Besoin en Eau de la Spéculation
float BesoinEnEau(const double * M,const double * Sp,const float BES_PARAM[NB_CULTURES+1][MT_IRRI]){
	float BES = 0.0;
	for(int i=0;i<NB_PARTIES;i++){
	    int Sp_i = (int)Sp[i];
	    int M_i = (int)M[i];
		BES += BES_PARAM[Sp_i][M_i];
	}
	return BES;
}


// Système dynamique
void dynamics_fd(const double * x, const double * u, double * image)
{
    const double* M_f = &u[2];
    const double* Sp_f = &u[2+NB_PARTIES];
    // const double* M_c = &u[2+2*NB_PARTIES];
    // const double* Sp_c = &u[2+3*NB_PARTIES];
	
    int NBCulture_f[NB_CULTURES] = {0};
    // int NBCulture_c[NB_CULTURES] = {0};
    	
    StockerNBCulture(Sp_f,NBCulture_f);
    // StockerNBCulture(Sp_c,NBCulture_c);
    	
    int NonEmptyParts_f = NonEmptyParcelles(Sp_f);
    // int NonEmptyParts_c = NonEmptyParcelles(Sp_c);
    	
    float PV_f = PrixVente(Sp_f, NBCulture_f, NonEmptyParts_f, CP_PV_PARAM);  
    // float PV_c = PrixVente(Sp_c, NBCulture_c, NonEmptyParts_c, CP_PV_PARAM); 
    	
    float BES_f = BesoinEnEau(M_f, Sp_f, BES_PARAM);  
    // float BES_c = BesoinEnEau(M_c, Sp_c, BES_PARAM);  

    float CP_f = CoutProduction(M_f, Sp_f, CP_PV_PARAM);  
    // float CP_c = CoutProduction(M_c, Sp_c, CP_PV_PARAM);  

	// Importance que cette fonction s'appelle dynamics !
	// x[0] = x1 = C(t), x[1] = x2 = P(t), x[2] = x3 = E(t), x[3] = x4 = h(t)
  

	image[0] = gamma_para*x[0] - CoutSurcreuse(x[1],u[0],SURCREUSE_PARAM) - CoutInvest(u[1]) - CP_f + PV_f - FIXED_COSTSR;
	    
	image[1] = u[0]*5.0;
	
	image[2] = x[2] + u[1];
	
	image[3] = x[3] + BES_f - FIXED_PLUIE;
}


// Constraints

// Check M[](t) constraint: M[](t) in {0,E(t),I+(t)}
bool CheckMConstraints(const double * M, int size, const double E, const double Iplus) {
  if (E == 3){
    return true;
  }

  for (int i = 0; i < size; i++) {
    if (M[i] != 0 && M[i] != E && M[i] != Iplus) {
      return false;
    }
  }
  return true;
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

double constraintsXU_fd(const double * x, const double * u)
{
    const double* M_f = &u[2];
    const double* Sp_f = &u[2+NB_PARTIES];
    // const double* M_c = &u[2+2*NB_PARTIES];
    // const double* Sp_c = &u[2+3*NB_PARTIES];
    int NBCulture_f[NB_CULTURES] = {0};
    // int NBCulture_c[NB_CULTURES] = {0};
    int NonEmptyParts_f = NonEmptyParcelles(Sp_f);
    // int NonEmptyParts_c = NonEmptyParcelles(Sp_c);
  
    // Check M[](t) constraint
    if (!CheckMConstraints(M_f, NB_PARTIES, x[2], u[1])) {
      return PLUS_INF;
    }
    
    // Check S(t) constraint: S(t) in [P(t)/5, 8]
    if(x[1]>u[0]*5){
      return PLUS_INF;
    }
  
	  // Check I+(t) constraint
	  if (!CheckIConstraint(x[2], u[1])) {
	    return PLUS_INF;
	  }

    // Check the constraint: h(t)+BES_f(t)*a+BES_c(t)*a<=S(t)*5
    float BES_f = BesoinEnEau(M_f, Sp_f, BES_PARAM);  
	  // float BES_c = BesoinEnEau(M_c, Sp_c, BES_PARAM); 
  
    if(x[3]+BES_f*BES_para>u[0]*5){
        return PLUS_INF;
    }

    // Check the constraint: h>=30 and M_k=0, then Sp_k=0
    for (int i=0;i<NB_PARTIES;i++){
      if(x[3]>=30 && M_f[i]==0 && Sp_f[i]!=0){
        return PLUS_INF;
      }
    }
    	
    // Check the constraint: C(t)
    float CS = CoutSurcreuse(x[1],u[0],SURCREUSE_PARAM);  
    float CI = CoutInvest(u[1]);                       
    float CP_f = CoutProduction(M_f, Sp_f, CP_PV_PARAM);    
    // float CP_c = CoutProduction(M_c, Sp_c, CP_PV_PARAM);
    float PV_f = PrixVente(Sp_f, NBCulture_f, NonEmptyParts_f, CP_PV_PARAM);
    // float PV_c = PrixVente(Sp_c, NBCulture_c, NonEmptyParts_c, CP_PV_PARAM);
  
    // Constraint 1: γ*C(t)-CS(t)-CI(t)-CP_f(t)-CP_c(t)-FC>=0
    if(gamma_para*x[0]-CS-CI-CP_f-FIXED_COSTSR<0){
  	    return PLUS_INF;
    }
    if(gamma_para*x[0]<0){
      return PLUS_INF;
    }
	
    // If all constraints are satisfied, return 1.0 to indicate feasibility
    return 1.0;
}

}
