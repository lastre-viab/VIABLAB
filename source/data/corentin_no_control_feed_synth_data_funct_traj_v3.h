/*! \file  Corentin4D_data.h
 *
 *
 *  \author: A. D�silles, LASTRE
 *  \brief  Fichier contant la d�finition d'un mod�le de viabilit� pour les tests
 *
 *  Ce fichier  contient les d�finitions  de tous les �l�ments qui d�crivent un
 *  un probl�me de viabilit� concret et  qui sont consid�r�s comme
 *  param�tres par le code; Ce fichier repr�sente ainsi une interface d'entr�e primitive pour le code.
 *
 *  Les �l�ments principaux d�finis ici  sont :
 *    -la  dynamique f(t,x,u)
 *    -les contraintes:
 *      -# sur l'�tat k(t,x)
 *      -# sur le controles U(x)
 *    - la cible c(t,x)
 *    -les fonctions d�finissant l'objectif  d'un probl�me d'optimisation
 *    	-# fonction l(t,x,u)
 *    	-# fonction m(t,x,u)
 *    - diff�rentes options d�finissant plus pr�cis�ment la nature du probl�me � r�soudre
 *    - diff�rentes options d�finissance la m�thode num�rique � utiliser
 *
 *    Voir plus loins dans les commentaires  du fichier la d�finition de chaque fonction et de chaque param�tre
 *
 *
 *
 */
//#include <algorithm>

double printVectorofVectors(std::vector<std::vector<double>> vv){
	cout << vv.size() << endl;
	for (int i = 0; i < vv.size(); i++) {
		for (
				auto it = vv[i].begin();
				it != vv[i].end(); it++)
			cout << *it << " ";
		cout << endl;
	}
	return 1;
}
double printVector(std::vector<double> v,std::string str){
	cout << "" << endl;
	cout << str << endl;
	for (unsigned int i=0; i< v.size(); i++)
		cout << v[i] << endl;
	return 1;
}


double printVectorString(std::vector<string> v){
	for (unsigned int i=0; i< v.size(); i++)
		cout << v[i] << endl;
	return 1;
}

double printVectorBool(std::vector<bool> v){
	for (unsigned int i=0; i< v.size(); i++)
		cout << v[i] << endl;
	return 1;
}

double print(double d){

	cout << d << endl;
	return 1;
}

double printInt(int i){

	cout << i << endl;
	return 1;
}

double printDoubleStar(double* vec){
	int size=sizeof(vec)+3;
	for (unsigned int i=0; i< size; i++)
		cout << vec[i] << endl;
	return 1;
}



double sumOf(std::vector<double> v, std::string str ){
	double sum=0.0;
	int size =0;
	size= v.size() ;
	int m=0;
	//cout << "La taille du vecteur est de : " << size << endl;
	for (m =0;m<size;++m){
		sum=sum+v[m];
	}

	cout << "La somme de " << str << " est égale à:" << sum <<endl;

	return sum;
}

double sumOfStar(double* vec, std::string str ){
	double sum=0.0;
	int size =sizeof(vec)+3;
	//double pr=printInt(size);
	//size= v.size() ;
	int m=0;
	//cout << "La taille du vecteur est de : " << size << endl;
	for (m =0;m<size;++m){
		sum=sum+vec[m];
	}

	//cout << "La somme de " << str << " est égale à:" << sum <<endl;

	return sum;
}

//Parametres;


//REGIONAUX
//double A_tot=238000.0;
//double tfinal=30.0;
double mu=0.07;



// double a_ol=0.041;
// double a_prot=0.007;
// double a_fod=0.412;
// double a_tub=0.001;
// double a_past=0.060; 

//Animal production(in kgN/livestock unit);

double C_mono_meat=9.9;
double C_mono_egg=2.3;
double C_rum_meat=3.8;
double C_rum_milk=19.5;
//FIXES
double phi_H=0.7;
double phi_R=0.19;
double phi_A_B=0.46;
double CN_H=10.0;
double CN_sol=10.0;
double CN_R=40.0;
double CN_A_B=16.0;
double xi_A=0.7;
double dep=10.0;
double percTimeGraz=0.5;
//Emissionfactors;
double epsilon_v_org=0.2;
double epsilon_v_syn=0.1;
double epsilon_Lapp=0.1;
double epsilon_N2Oapp=0.01;
double epsilon_N2OSM=0.03;


//Area;
//double a_cer=1-a_fod-a_past-a_tub-a_ol-a_prot;
//double A_cro=(a_cer+a_fod+a_ol+a_prot+a_tub)*A_tot;
//double A_F=(a_cer+a_fod+a_ol+a_tub)*A_tot;

double tau_M_init=0.0;

//double HI[]={0.46,0.8,0.35,0.58,0.8,0.6};
//double SR[]={4.4,5.0,5.0,5.3,0.8,5.0};

//double eta_max[]={50.0,100.0,25.0,30.0,50.0,340.0};
//double delta[]={0.321,0.292,0.154,0.191,0.296,1.754};//to compute
//double sigma_H[]={1.8,2.4,3.6,4.4,1.2,0.35};
//double sigma_R_R[]={1.0,1.6,1.2,0.8,1.2,0.14};
//double sigma_R_A[]={0.7,2.5,1.2,0.8,1.2,0.26};
//double psi[]={0.0,63.0,27.0,77.0,14.0,0.0};
//double a[]={0.47,0.41,0.04,0.01,0.06,0.01};

//double nu_H_H[]={0.4,0.0,0.15,0.2,0.0,0.75};
//double nu_T_H[]={0.0,0.0,0.1,0.0,0.0,0.0};
//double nu_H_A[]={0.4,0.95,0.4,0.7,1.0,0.0};
//double nu_T_A[]={0.05,0.0,0.2,0.0,0.0,0.1};

//double nu_T_H_A[6][7]={{0.4,0.05,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.95,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.2,0.4,0.0},{0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.7,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.1}};

//std::vector<double>  beta_mono={33.0,0.01,0.0,0.0,39.0,0.5,2.0};
//std::vector<double>  beta_rum={8.2,0.01,43.0,15.0,39.0,0.5,0.0};

//Created for data import
////Crop
struct Crop {

	std::string name;
	double area;
	double yield;
	double SR;
	double HI;
	double sigma_H;
	double sigma_R_A;
	double sigma_R_R;
	double BNF;
	double phi_R_A;
	double phi_R_R;
	double CN_R_A;
	double CN_R_R;
	bool fertilizedOrNot;
	bool annualCropOrNot;
	bool croplandOrNot;
	bool cultivatedOrNot;
	bool controlOrNot;
};

std::vector<Crop> importCropData (std::string stringName){


	string file = stringName;
	ifstream fin(file); //opening the file.

	const char delim = ',';
	std::vector<Crop> crops;// what will be returned

	//ignore first line
	string line;
	getline(fin, line);
	if (fin) {
		cout << "" << "\n";
		cout << "File opened!" << endl;
	}
	int j=0;

	while(getline(fin,line)){
		istringstream ss(line);
		Crop crop;

		string str;
		std::getline(ss, crop.name,delim);

		std::getline(ss, str,delim);
		crop.area=stod(str);
		std::getline(ss, str,delim);
		crop.yield=stod(str);
		std::getline(ss, str,delim);
		crop.SR=stod(str);
		std::getline(ss, str,delim);
		crop.HI=stod(str);
		std::getline(ss, str,delim);
		crop.sigma_H=stod(str);
		std::getline(ss, str,delim);
		crop.sigma_R_A=stod(str);
		std::getline(ss, str,delim);
		crop.sigma_R_R=stod(str);
		std::getline(ss, str,delim);
		crop.BNF=stod(str);
		std::getline(ss, str,delim);
		crop.phi_R_A=stod(str);
		std::getline(ss, str,delim);
		crop.phi_R_R=stod(str);
		std::getline(ss, str,delim);
		crop.CN_R_A=stod(str);
		std::getline(ss, str,delim);
		crop.CN_R_R=stod(str);
		std::getline(ss, str,delim);
		istringstream(str) >> std::boolalpha >> crop.fertilizedOrNot;
		//crop.fertilizedOrNot=boost::lexical_cast<bool>(str);
		std::getline(ss, str,delim);
		istringstream(str) >> std::boolalpha >> crop.annualCropOrNot;
		//crop.annualCropOrNot=boost::lexical_cast<bool>(str);
		std::getline(ss, str);
		istringstream(str) >> std::boolalpha >> crop.croplandOrNot;
		//crop.croplandOrNot=boost::lexical_cast<bool>(str);
		std::getline(ss, str);
		istringstream(str) >> std::boolalpha >> crop.cultivatedOrNot;
		//crop.croplandOrNot=boost::lexical_cast<bool>(str);
		std::getline(ss, str); 
		istringstream(str) >> std::boolalpha >> crop.controlOrNot;
		//crop.croplandOrNot=boost::lexical_cast<bool>(str);


		crops.push_back (crop);
		j += 1; //increment number of lines
		//cout << j << endl;

	}
	cout << "Table crop data" << "\n";
	for (unsigned int i=0; i< crops.size(); i++)
		cout << crops[i].name << ' '
		<< crops[i].area << ' '
		<< crops[i].yield << ' '
		<< crops[i].SR << ' '
		<< crops[i].HI << ' '
		<< crops[i].sigma_H << ' '
		<< crops[i].sigma_R_A << ' '
		<< crops[i].sigma_R_R << ' '
		<< crops[i].BNF << ' '
		<< crops[i].phi_R_A << ' '
		<< crops[i].phi_R_R << ' '
		<< crops[i].CN_R_A << ' '
		<< crops[i].CN_R_R << ' '
		<< crops[i].fertilizedOrNot << ' '
		<< crops[i].annualCropOrNot << ' '
		<< crops[i].croplandOrNot << ' '
		<< crops[i].cultivatedOrNot << ' '
		<< crops[i].controlOrNot << ' '
		<< "\n";


	fin.close(); //closing the file
	std::cout << "Number of entries: " << j << endl;

	return crops;

}

std::string REG_ID="22359";
std::string fileNameCrop="cropData_" + REG_ID + ".txt";
std::vector<Crop> cropData=importCropData(fileNameCrop);

std::vector<double> collect(double Crop::* f, std::vector<Crop> const& v) {
	std::vector<double> output;
	for (auto const& elem : v) {
		output.push_back(elem.*f);
	}
	return output;
}

std::vector<string> collectString(string Crop::* f, std::vector<Crop> const& v) {
	std::vector<string> output;
	for (auto const& elem : v) {
		output.push_back(elem.*f);
	}
	return output;
}

std::vector<bool> collectBool(bool Crop::* f, std::vector<Crop> const& v) {
	std::vector<bool> output;
	for (auto const& elem : v) {
		output.push_back(elem.*f);
	}
	return output;
}

std::vector<string> cropName=collectString(&Crop::name,cropData);
std::vector<double> A_init=collect(&Crop::area,cropData);
std::vector<double> eta_max=collect(&Crop::yield,cropData);
std::vector<double> SR=collect(&Crop::SR,cropData);
std::vector<double> HI=collect(&Crop::HI,cropData);
std::vector<double> sigma_H=collect(&Crop::sigma_H,cropData);
std::vector<double> sigma_R_A=collect(&Crop::sigma_R_A,cropData);
std::vector<double> sigma_R_R=collect(&Crop::sigma_R_R,cropData);
std::vector<double> psi=collect(&Crop::BNF,cropData);
std::vector<double> phi_R_A=collect(&Crop::phi_R_A,cropData);
std::vector<double> phi_R_R=collect(&Crop::phi_R_R,cropData);
std::vector<double> CN_R_A=collect(&Crop::CN_R_A,cropData);
std::vector<double> CN_R_R=collect(&Crop::CN_R_R,cropData);
std::vector<bool> fertilizedOrNot=collectBool(&Crop::fertilizedOrNot,cropData);
std::vector<bool> annualCropOrNot=collectBool(&Crop::annualCropOrNot,cropData);
std::vector<bool> croplandOrNot=collectBool(&Crop::croplandOrNot,cropData);
std::vector<bool> cultivatedOrNot=collectBool(&Crop::cultivatedOrNot,cropData);
std::vector<bool> controlOrNot=collectBool(&Crop::controlOrNot,cropData);

vector<double> changePercentGrassland(std::vector<double> Area){
	vector<double> A=Area;
	double A_tot=sumOf(Area," somme de surface total initiale ");
	std::list<std::string> grassland = { "HPNatGrassland"};
	double share=0.1;
	for (unsigned int i=0; i< Area.size(); i++){

		if (std::find(grassland.begin(), grassland.end(), cropName[i]) != grassland.end()) {
			A[i]=share/(1-share)*(A_tot-Area[i]);
			cout << "New grassland area : " << A[i] << endl;
		}

	}



	return A;
}

//std::vector<double> A=changePercentGrassland(A_init);
std::vector<double> A=A_init;
double computeWithBool (std::vector<double> v,std::vector<bool> b,std::string ofWhat){
	double d=0;

	for (unsigned int i=0; i< v.size(); i++)
		if (b[i] == true)
			d=d+v[i];

	cout << ofWhat << ": " << d << endl;
	return d;
}

double A_cro=computeWithBool(A,croplandOrNot,"Cropland area");
double A_F=computeWithBool(A,fertilizedOrNot,"Fertilized area");


///Diet
struct Diet {

	std::string feedCat;
	double mono;
	double rum;
};

std::vector<Diet> importDietData (std::string stringName){


	string file = stringName;
	ifstream fin(file); //opening the file.

	const char delim = ',';
	std::vector<Diet> diets;// what will be returned

	//ignore first line
	string line;
	getline(fin, line);
	if (fin) {
		cout << "" << "\n";
		cout << "File opened!" << endl;
	}
	int j=0;

	while(getline(fin,line)){
		istringstream ss(line);
		Diet diet;

		string str;
		std::getline(ss, diet.feedCat,delim);

		std::getline(ss, str,delim);
		diet.mono=stod(str);
		std::getline(ss, str);
		diet.rum=stod(str);

		diets.push_back (diet);
		j += 1; //increment number of lines


	}

	for (unsigned int i=0; i< diets.size(); i++)
		cout << diets[i].feedCat << ' '
		<< diets[i].mono << ' '
		<< diets[i].rum << ' '
		<< "\n";


	fin.close(); //closing the file
	std::cout << "Number of entries: " << j << endl;

	return diets;

}
//std::string fileNameDiet="coefFeedLvstck_" + REG_ID + ".txt";
std::string fileNameDiet="coefFeedLvstck.txt";
std::vector<Diet> dietData=importDietData(fileNameDiet);

std::vector<double> collectDiet(double Diet::*f, std::vector<Diet> const& v) {
	std::vector<double> output;
	for (auto const& elem : v) {
		output.push_back(elem.*f);
	}
	return output;
}

std::vector<double> beta_mono=collectDiet(&Diet::mono,dietData);
std::vector<double> beta_rum=collectDiet(&Diet::rum,dietData);


///Allocation coefficients livestock

/*
struct Nu {

    std::string cropCode;
    double concentrates;
    double forages;
};

std::vector<Nu> importNuData (char stringName[]){


 	string file = stringName;
  	ifstream fin(file); //opening the file.

  	const char delim = ',';
  	std::vector<Nu> nus;// what will be returned

  	//ignore first line
  	string line;
  	getline(fin, line);
	if (fin) {
	cout << "File opened!" << endl;
	}
    	int j=0;

	while(getline(fin,line)){
        istringstream ss(line);
        Nu nu;

	string str;
	std::getline(ss, nu.cropCode,delim);

	std::getline(ss, str,delim); 
	nu.concentrates=stod(str);
	std::getline(ss, str); 
	nu.forages=stod(str);

    nus.push_back (nu);
    j += 1; //increment number of lines


  	}

  		for (unsigned int i=0; i< nus.size(); i++)
        		cout << nus[i].cropCode << ' '
             			<< nus[i].concentrates << ' ' 
             			<< nus[i].forages << ' ' 
             			<< "\n";


  		fin.close(); //closing the file
  		std::cout << "Number of entries: " << j << endl;

  	return nus;

}

std::vector<Nu> nuData=importNuData("allocationCoefLivestock.txt");

std::vector<double> collectNu(double Nu::*f,  std::vector<Nu> const& v) {
    std::vector<double> output;
    for (auto const& elem : v) {
        output.push_back(elem.*f);
    }
    return output;
}

std::vector<double> nu_concentrates=collectNu(&Nu::concentrates,nuData);
std::vector<double> nu_forages=collectNu(&Nu::forages,nuData);

//double length = A.size(), width = beta_mono.size();
//vector<vector<double>> nu_lvstck(length, vector<double> (width, 0));


vector<vector<double>> asMatrix(std::vector<double> v1, std::vector<double> v2){
std::vector<std::vector<double> > output;
std::vector<std::vector<double> > transpose(v1.size(), vector<double> (2, 0));
output.push_back(v1);
output.push_back(v2);

for (int i = 0; i < v1.size(); ++i) {
      for (int j = 0; j < 2; ++j) {
			transpose[i][j] = output[j][i];
      }
   }

return transpose;
}

vector<vector<double>> nu_lvstck=asMatrix(nu_concentrates,nu_forages);
 */


struct feedCatList {

	std::string feedCat;
	bool importableOrNot;
	bool feedFoodCompOrNot;
};

struct Nu {

	std::string cropName;
	std::string feedCat;
	double coefAlloc;
};


std::vector<Nu> importNuData (char stringName[]){


	string file = stringName;
	ifstream fin(file); //opening the file.

	const char delim = ',';
	std::vector<Nu> nus;// what will be returned

	//ignore first line
	string line;
	getline(fin, line);
	if (fin) {
		cout << "" << "\n";
		cout << "File opened!" << endl;
	}
	int j=0;

	while(getline(fin,line)){
		istringstream ss(line);
		Nu nu;

		string str;
		std::getline(ss, nu.cropName,delim);

		std::getline(ss, nu.feedCat,delim);
		std::getline(ss, str);
		nu.coefAlloc=stod(str);

		nus.push_back (nu);
		j += 1; //increment number of lines


	}
	cout << "Table allocation coefficient for livestock" << "\n";
	for (unsigned int i=0; i< nus.size(); i++)
		cout << nus[i].cropName << ' '
		<< nus[i].feedCat << ' '
		<< nus[i].coefAlloc << ' '
		<< "\n";


	fin.close(); //closing the file
	std::cout << "Number of entries: " << j << endl;

	return nus;

}

std::vector<feedCatList> importFeedCatListData (char stringName[]){


	string file = stringName;
	ifstream fin(file); //opening the file.

	const char delim = ',';
	std::vector<feedCatList> FCLs;// what will be returned

	//ignore first line
	string line;
	getline(fin, line);
	if (fin) {
		cout << "" << "\n";
		cout << "File opened!" << endl;
	}
	int j=0;

	while(getline(fin,line)){
		istringstream ss(line);
		feedCatList FCL;

		string str;
		std::getline(ss, FCL.feedCat,delim);

		std::getline(ss, str,delim);
		istringstream(str) >> std::boolalpha >> FCL.importableOrNot;

		std::getline(ss, str,delim); 
		istringstream(str) >> std::boolalpha >> FCL.feedFoodCompOrNot;


		FCLs.push_back (FCL);
		j += 1; //increment number of lines


	}
	cout << "Table feed categories" << "\n";
	for (unsigned int i=0; i< FCLs.size(); i++)
		cout << FCLs[i].feedCat << ' '
		<< FCLs[i].importableOrNot<< ' '
		<< FCLs[i].feedFoodCompOrNot<< ' '
		<< "\n";


	fin.close(); //closing the file
	std::cout << "Number of entries: " << j << endl;

	return FCLs;

}


std::vector<feedCatList> feedCatListData=importFeedCatListData("feedCat_feedLevel1.txt");
std::vector<Nu> nuData=importNuData("allocationCoefLivestock_feedLevel1_linear2.txt");

std::vector<string> collectFeedCatList(string feedCatList::*f,  std::vector<feedCatList> const& v) {
	std::vector<string> output;
	for (auto const& elem : v) {
		output.push_back(elem.*f);
	}
	return output;
}

std::vector<double> collectNu(double Nu::*f,  std::vector<Nu> const& v) {
	std::vector<double> output;
	for (auto const& elem : v) {
		output.push_back(elem.*f);
	}
	return output;
}

std::vector<bool> collectBoolFeed(bool feedCatList::*f,  std::vector<feedCatList> const& v) {
	std::vector<bool> output;
	for (auto const& elem : v) {
		output.push_back(elem.*f);
	}
	return output;
}


std::vector<double> coefAlloc=collectNu(&Nu::coefAlloc,nuData);
std::vector<string> feedCat=collectFeedCatList(&feedCatList::feedCat,feedCatListData);
std::vector<bool> feedFoodCompOrNot=collectBoolFeed(&feedCatList::feedFoodCompOrNot,feedCatListData);
std::vector<bool> importableOrNot=collectBoolFeed(&feedCatList::importableOrNot,feedCatListData);


//double pvs1=printVectorString(feedCat);
//double pvs2=printVectorString(cropName);
//double pv1=printVector(coefAlloc);

std::vector<vector<double>> unmelt_nu(std::vector<string> cropName, std::vector<string> feedCat, std::vector<double> nu){
	std::vector<std::vector<double> > output(cropName.size(),std::vector<double>(feedCat.size(),0));
	int k=0;

	for (int j = 0; j < feedCat.size(); ++j) {
		for (int i = 0; i < cropName.size(); ++i) {

			output[i][j]=nu[k];
			k=k+1;
		}
	}
	return output;

}

std::vector<std::vector<double>> nu_lvstck=unmelt_nu(cropName,feedCat,coefAlloc);
double pvv1=printVectorofVectors(nu_lvstck);



///Allocation coefficients humans
struct NuH {

	std::string cropName;
	double nu;

};

std::vector<NuH> importNuHData (char stringName[]){


	string file = stringName;
	ifstream fin(file); //opening the file.

	const char delim = ',';
	std::vector<NuH> nuHs;// what will be returned

	//ignore first line
	string line;
	getline(fin, line);
	if (fin) {
		cout << "" << "\n";
		cout << "File opened!" << endl;
	}
	int j=0;

	while(getline(fin,line)){
		istringstream ss(line);
		NuH nuH;

		string str;
		std::getline(ss, nuH.cropName,delim);
		std::getline(ss, str);
		nuH.nu=stod(str);

		nuHs.push_back (nuH);
		j += 1; //increment number of lines


	}
	/*
  		for (unsigned int i=0; i< nuHs.size(); i++)
        		cout << nuHs[i].cropName << ' '
             			<< nuHs[i].nu << ' ' 
             			<< "\n";


  		fin.close(); //closing the file
  		std::cout << "Number of entries: " << j << endl;
	 */
	return nuHs;

}

std::vector<NuH> nuHData=importNuHData("allocationCoefHumans.txt");

std::vector<double> collectNuH(double NuH::*f,  std::vector<NuH> const& v) {
	std::vector<double> output;
	for (auto const& elem : v) {
		output.push_back(elem.*f);
	}
	return output;
}

std::vector<double> nu_H=collectNuH(&NuH::nu,nuHData);



std::vector<double> computeDelta(){
	std::vector<double> output;
	double tmp;
	int j=0;
	//cout << eta_max.size()<<endl;
	for (long unsigned int i = 0; i < eta_max.size(); i++) {
		j+=1;
		tmp=1/(sigma_H[i]+(sigma_R_A[i]*(1-HI[i])/HI[i])+(sigma_R_R[i]*1.65/(SR[i]*HI[i])));
		output.push_back(tmp);
		//cout << j << endl;
	}
	return output;
}
//std::vector<double> delta=1/(sigma_H+(sigma_R_A*(1-HI)/HI)+(sigma_R_R*1.65/(SR*HI)));
std::vector<double> delta=computeDelta();

//double pv1=printVector(beta_rum);
//double pvv1=printVectorofVectors(nu_lvstck);
//double pv1=printVector(nu_H);
//double pv2=printVector(delta);





//Functions to debug - check ODM




double beta_mono_tot=sumOf(beta_mono,"apport d'azote des monogastriques en kg/UGB");
double beta_rum_tot=sumOf(beta_rum,"apport d'azote des ruminants en kg/UGB");
double A_tot=sumOf(A,"surface totale en hectares");
double E_H=5*2*A_tot*0.1;

int sizeCrop =  eta_max.size();
int sizeFeed = beta_rum.size();

double maxOf(double a, double b){
	double max=0.0;
	if (a>b) {
		max=a;
	} else {
		max=b;
	}
	return max;

}

double minOf(double a, double b){
	double min=0.0;
	if (a>b) {
		min=b;
	} else {
		min=a;
	}
	return min;

}



//Perturbations(inkgN);


//double nbMonoMax=400000.0;
//double nbRumMax=200000.0;

double nbMonoMax=1.5*A_tot;
double nbRumMax=1.0*A_tot;

vector<double> initialisationPertFeedImp(double coefperturb){
	vector<double> I;
	double tmp;

	for (long unsigned int i =0;i<beta_rum.size();i++){
		tmp=(beta_mono[i]*nbMonoMax+beta_rum[i]*nbRumMax)*coefperturb;
		I.push_back(tmp);
		//cout << tmp << endl;
	}
	return I;
}

vector<double> I=initialisationPertFeedImp(0.0);
double S_a=0.0*A_F;//maximum and initial perturbation

double pv_I=printVector(I,"Feed import per category (kg nitrogen)");





#ifndef CORENTIN4D_DATA_H_include
#define CORENTIN4D_DATA_H_
/*! \var dim
 *  \brief State  dimension
 */
const int dim=4;

/*!
 * \var dynType  defines the type of dynamics in the model
 *      1 or CC : continuous in time and space
 *      2 or DC : discrete time continuous space
 *      3 or DD : discrete time discrete space
 *      4 or HD : hybrid \todo
 */
const int dynType=DC;

/*! \var dicret_type
 *  \brief the  type of discretization scheme used  for the dynamics
 *  EE or 1 = Euler Explicit  scheme
 *  EI or 2 = Euler Implicit
 *  RK2I or 3 = RK2 Implicit (RK2 for -F)
 *  RK2E or 4 = RK2 Explicit (RK2 for F)
 *  RK4I or 5 = RK2 Implicit (RK4 for -F)
 *  RK4E or 6 = RK2 Explicit (RK4 for F)
 */
const int discret_type=0;

/*! \var STATE_MIN[dim]
 *  \brief min values  for  state vector components
 *   \see gridParams
 */

double STATE_MIN[dim]={500.0,0.0,0.0,0.0};

/*! \var STATE_MAX[dim]
 *  \brief max values  for  state vector components
 *   \see gridParams
 */
double STATE_MAX[dim]={2000.0,200.0*A_tot,nbMonoMax,nbRumMax};  

/*!
 * \var nbPointsState[dim]
 * number of  discretization points for the state variable
 * \see gridParams
 */
//unsigned long long int nbPointsState[dim] = {20,21,81,61,tfinal+2}; 
unsigned long long int nbPointsState[dim] = {16,21,51,76}; //4 919 376 points ou 1 250 256 points
unsigned long long int dirTramage =1;

/*!
 * \var periodic[dim]
 * \brief indicator of periodicity
 * 1= periodic variable
 * 0=non periodic variable
 * \see gridParams
 */
int periodic[dim]={0,0,0,0};

int saveProjection=0;
unsigned long long int projection[dim]={0,0,0,0};
int saveBoundary=1;


/**************************
 *        Controls         *
 ***************************/

bool globalControlOrNot=false;


/*! \var dimc
 *  \brief Control variable   dimension
 */

const int dimc=0;




/*! \var CONTROL_MIN[dimc]
 *  \brief minimum values  for  control vector components
 *   \see controlParams
 */


double CONTROL_MIN[dimc]={};




/*! \var CONTROL_MAX[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */

double CONTROL_MAX[dimc]={};


/*! 
 * \var nbPointsControl[dimc]
 * number of  discretization points for the control variable
 * \see controlParams
 */

unsigned long long int nbPointsControl[dimc] =  {};


std::list<std::string> cereals = { "Barley", "GrainMaize","OtherCereals","SoftWheat"};
std::list<std::string> oilProt = { "ProtCrops", "Oilseeds"};
std::list<std::string> legMeadows = { "LegumeMeadows", "MixMeadows"};
std::list<std::string> grassFod = { "MaizeFodder", "GrassMeadows"};
std::list<std::string> others = { "PotatoConsumption", "RootVegetables","LegumesVegetables"};

double A_control=computeWithBool(A,controlOrNot,"Control area");

double getUInit(std::list<std::string> list) {

	double sum=0.0;
	int i;

	for (i =0;i<sizeCrop;++i){
		if (cultivatedOrNot[i]==true){
			if (controlOrNot[i]==true){

				if(std::find(list.begin(), list.end(), cropName[i]) != list.end()) {
					sum+=A[i];

				}

			}
		}
	};

	sum=sum/A_control;
	cout << "Valeur du controle:" << sum << endl;
	return sum;

}
double u_cereals=getUInit(cereals);
double u_oilProt=getUInit(oilProt);
double u_legMeadows=getUInit(legMeadows);
double u_grassFod=getUInit(grassFod);
double u_others=getUInit(others);


double residues_tot_init() {
	double sum=0.0;
	int i;


	for (i =0;i<sizeCrop;++i){
		if (cultivatedOrNot[i]==true){

			sum+=eta_max[i]*A[i]*(sigma_R_A[i]*((1-HI[i])/HI[i])+1.65*sigma_R_R[i]*(1/(SR[i]*HI[i])));

		}
	};
	cout << sum << endl;
	return sum;
}
double res_init=residues_tot_init();

//const int nbTrajs=0;
//double initPoints[dim*nbTrajs]={};
//double initControls[dimc*nbTrajs]={};

const int nbTrajs=5;
double initPoints[dim*nbTrajs]={870.0,res_init,1.5*A_tot,0.25*A_tot,//initial state
		870.0,res_init,1.5*A_tot,1*A_tot,//carnivore
		870.0,res_init,0*A_tot,1*A_tot,//carnivore-only ruminants
		870.0,res_init,1.5*A_tot,0*A_tot,//carnivore-only monogastrics
		870.0,res_init,0*A_tot,0.02*A_tot};//vegan
double initControls[dimc*nbTrajs]={};

/*!
 * sélection de type de reconstruction de trajectoire
 *  VD= 1, viable par défaut, on sélectionne le premier contrôle viable trouvé
 *  VL= 2, viable lourd: le contrôle reste constant tant qu'il viable;
 *  cette méthode nécessite une initialisation de contrôle
 */
int typeTraj=VD;



/*!
 * target = 1
 * departure =0;
 * Ce parametre determine le sens des trajectoires
 */
int target_or_departure_problem=0;







/*!
 * \var sortieOK[dim]
 * \brief Tableau qui définit les modalités d'interprétation de
 * la sortie du domaine de calcul lors du calcul de noyau de viabilité. Permet d'indiquer
 * si le domaine réel est infini le long d'une direction , alors que le domaine de calcul est forcément borné.
 * Si la valeur correspondante est 1 alors  si le successeur d'un point x a  la composante correspondante au-delà
 * des bornes du domaine de calcul il sera considéré comme viable.
 */
int sortieOK[dim]={1,1,1,1};

string getPrefix(bool switchVal){
	string str;

	if (switchVal==true){
		str="corentin_control_crop_comp_feed_synth_"+REG_ID;
	} else {
		str= "corentin_no_control_feed_synth_"+REG_ID;
	}

	return str;
}
string prefix=getPrefix(globalControlOrNot);




/*!
 * \var globalDeltaT
 *  bool�en indique si le pas de temps  doit �re choisi globalement
 *  ou localement pour les algorithmes de viabilit�
 */
bool globalDeltaT=false;


/*!
 * \var T maximum time horizon for the study
 */
double T=0.0004;


/*!
 * Sélection de la méthode de représentation de l'ensemble
 * Ceparamètre détermine quelle classe sera utilisée pour les calculs
 *
 *    - BS = BitSet, représentation par fonction caractéristique , pour
 *                   tout type de calculs
 *    - MM = MicroMacro, représentation par valeurs réelles, pour les calculs
 *          d'ensembles épigraphiques , associés aux système micro-macro
 */
int gridMethod=BS;

/*!
 * Sélection de l'ensemble à calculer
 */
int setType=VIAB;




double n_OM_mineral_PSU(double *x){
	double n_OM=0.0;

	double sewage_sludge_min=E_H*(1-(phi_H*CN_H/CN_sol));
	double residues_min=x[1]*(1-(phi_R*CN_R/CN_sol));
	double lvstck_effluent_min=xi_A*(1-0.5)*(beta_mono_tot*x[2]+beta_rum_tot*percTimeGraz*x[3])*(1-(phi_A_B*CN_A_B/CN_sol));


	n_OM=(sewage_sludge_min+lvstck_effluent_min+residues_min)*(1-epsilon_v_org)/A_cro;
	//cout << n_OM << endl;
	return n_OM;
}

double yield_i(int i,double *x){
	double y=0.0;
	double n_OM=0.0;
	n_OM=n_OM_mineral_PSU(x);
	//cout << n_OM << endl;

	y=minOf(eta_max[i],delta[i]*((mu*x[0]+dep+n_OM+(S_a/A_F*(1-epsilon_v_syn)))*(1-epsilon_Lapp-epsilon_N2Oapp-epsilon_N2OSM)+psi[i]));
	//cout << "eta_max " << eta_max[i] << endl;
	//cout << "actual yield " << y << endl;
	return y;
}



double residues_tot(double *x) {
	double sum=0.0;
	int i;


	for (i =0;i<sizeCrop;++i){
		if (cultivatedOrNot[i]==true){

			sum+=yield_i(i,x)*A[i]*(sigma_R_A[i]*((1-HI[i])/HI[i])+1.65*sigma_R_R[i]*(1/(SR[i]*HI[i])));

		}
	};

	return sum;
}

double harvToHumans_PSU(double *x) {//good ODM

	double sum=0.0;
	int i;




	for (i =0;i<sizeCrop;++i){
		if (cultivatedOrNot[i]==true){

			sum+=yield_i(i,x)*A[i]*sigma_H[i]*(nu_H[i]);

		}
	};



	sum=sum/A_tot;
	return sum;

}





double feed_av_k(int k,double *x, vector<double> I,std::string str) {

	double sum=0.0;
	int i;
	int j;
	int l;
	std::vector<double> beta;
	std::vector<double> beta_other;
	double ratio=1.0;
	if (str == "rum")
	{
		beta=beta_rum;
		beta_other=beta_mono;
		j=3;
		l=2;
	};
	if (str == "mono")
	{
		beta=beta_mono;
		beta_other=beta_rum;
		j=2;
		l=3;
	};




	for (i =0;i<sizeCrop;++i){

		if (cultivatedOrNot[i]==true){

			sum=sum+(yield_i(i,x)*A[i]*sigma_H[i]*nu_lvstck[i][k]);

		}
	}




	if (importableOrNot[k]==true){
		sum=sum+I[k];
	}

	ratio=((beta[k]*x[j])/((beta_other[k]*x[l])+(beta[k]*x[j])));
	if (isnan(ratio)){
		ratio=0.0;
	};

	sum=sum*ratio;

	return sum;
}




double feed_av(double *x, vector<double> I, std::string str) {

	double sum=0.0;
	int k;
	std::vector<double> beta;
	if (str == "rum")
		beta=beta_rum;
	if (str == "mono")
		beta=beta_mono;


	for (k=0;k<sizeFeed;++k){
		if (beta[k]>0){
			sum=sum+feed_av_k(k,x,I,str);
		}
	}

	return sum;
}

double tau_M(double *x,vector<double> I, std::string str){
	double tau_M=0.0;
	double shortage=0.0;
	int i=0;
	double beta_tot=0.0;
	if (str == "mono")
	{
		i=2;
		beta_tot=beta_mono_tot;
	};
	if (str == "rum")
	{
		i=3;
		beta_tot=beta_rum_tot;
	};

	shortage=(-1.0*maxOf(0.0,((beta_tot*x[i])-feed_av(x,I,str))/(beta_tot*x[i])));
	if (isnan(shortage)){
		shortage=0.0;

	}
	if (shortage>=0.0){
		tau_M=tau_M_init;
	} else {
		tau_M=shortage;
	};
	return tau_M;
}



double feed_surplus(double *x,vector<double> I){
	double sum=0.0;
	int k;
	for (k =0;k<I.size();++k){
		if (feedFoodCompOrNot[k]==true){
			sum=sum+maxOf(0.0,(feed_av_k(k,x,I,"rum")+feed_av_k(k,x,I,"mono"))-(beta_mono[k]*x[2]+beta_rum[k]*x[3]));
		}
	}
	return sum;
}




/*!
 * Definition of the dynamics  and associated functions and constants
 */

void dynamics(double * x, double *u, double * image)
{


	double sewage_sludge_org=E_H*phi_H*CN_H/CN_sol;
	double residues_org=x[1]*phi_R*CN_R/CN_sol;
	double lvstck_effluent_org=xi_A*(1.0-0.5)*(beta_mono_tot*x[2]+beta_rum_tot*percTimeGraz*x[3])*phi_A_B*CN_A_B/CN_sol;

	image[0]=x[0]*(1.0-mu)+(residues_org+lvstck_effluent_org+sewage_sludge_org)/A_cro;

	image[1]=residues_tot(x);

	image[2]=x[2]*(1.0+tau_M(x,I,"mono"));

	image[3]=x[3]*(1.0+tau_M(x,I,"rum"));

	//image[4]=x[4]+1.0;



}

/*!
 * Definition of the dynamics  and associated functions and constants
 */

double dynConstraintsForTraj(double * x, double * image)
{

	// Cette fonction doit renvoyer 1.0 si le point image respecte une contrainte par rapport au point de d�part, x
	// et PLUS_INF sinon
	double res = 1.0;

	res =  ((x[2]<=image[2]) && (x[3]<= image[3])) ? 1.0 : PLUS_INF;
	return res;


}


/*****************************************
 *  Definition of constraints and target *
 *****************************************/


/*!
 * \brief Function  defining the mixed  constraints
 *
 * This function defines the set U(x) for admissible controls as function of the state
 * @param x state variable
 * @param u control variable
 * @return  value that caraterise the constraints set
 */


inline double constraintsX( double * x)
{	

	double anim_prod_PSU=((C_mono_egg+C_mono_meat)*x[2]+(C_rum_meat+C_rum_milk)*x[3])/A_tot;


	//double feed_surplus_PSU=feed_surplus_no_control(x,I)/A_tot;
	//double con_1 =anim_prod_PSU+harvToHumans_PSU_no_control(x)+feed_surplus_PSU;


	double feed_surplus_PSU=feed_surplus(x,I)/A_tot;
	double con_1 =anim_prod_PSU+harvToHumans_PSU(x)+feed_surplus_PSU;


	//cout << "Anim_prod_PSU: " << anim_prod_PSU << endl;
	//cout << "Feed surplus PSU " << feed_surplus_PSU << endl;
	//cout << "Crop_prod_PSU: " << harvToHumans_PSU(x,u) << endl;
	//cout << "Contrainte: " << con_1 << endl;
	double res = 1.0;
	//if (x[4]<tfinal){
	if(con_1 >= 10.0) {res = 1.0;} else {res = PLUS_INF;};
	//};
	return res;
}

//#include "corentin_no_control_unused.h"
#include "corentin_no_control_traj_unused.h"
//#include "corentin_control_unused.h"
//#include "corentin_control_traj_unused.h"

#endif /* TESTDATA_H_ */
