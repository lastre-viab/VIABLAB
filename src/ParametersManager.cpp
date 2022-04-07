/*
 * ParametersManager.cpp
 *
 *  Created on: 27 juil. 2021
 *      Author: adesi
 */

#include "../include/ParametersManager.h"


ParametersManager::ParametersManager() {
	// TODO Auto-generated constructor stub

}

ParametersManager::ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams *cp,systemParams *sp)
{
	gridParameters =gp;
	algoParameters=avp;
	controlParameters=cp;
	systemParameters=sp;
}

ParametersManager::ParametersManager(gridParams *gp, algoViabiParams *avp, controlParams * cp, string controlParamsFile, systemParams *sp)
{
	gridParameters =gp;
	algoParameters=avp;
	controlParameters=  cp;
	systemParameters=sp;
	this->readControlParametersFromJson(controlParamsFile);
}

gridParams * ParametersManager::getGridParameters(){
	return gridParameters;
}
algoViabiParams * ParametersManager::getAlgoParameters(){
	return algoParameters;
}
controlParams * ParametersManager::getControlParameters(){
	return controlParameters;
}
systemParams * ParametersManager::getSystemParameters(){
	return systemParameters;
}

void ParametersManager::readControlParametersFromJson( string controlParamsFile)
{
	string line;
	ostringstream os;

	string tempStr, tmpStr;
	ptree dataRoot;
	string sedCommand;

	string input_tempfile="../INPUT/"+ controlParamsFile;
	//sedCommand="sed  -e 's/[ ]*#.*//g' -e '/^$/d' ../EXEC_DATA/"+ farmSimulationFile+" > "+input_tempfile;
	//unsigned long long int  zz=system(sedCommand.c_str()); zz++;

	//- chargement des donnees depuis le fichier json via le parseur boost
	read_json(input_tempfile, dataRoot);


	/*
	 * initialisation des paramètres des contrôles
	 */
	int dimc = dataRoot.get<int>("CONTROL_DIMENSION", 1);
	controlParameters->DIMC=(unsigned long long int) dimc;
	/*
	 * paramètres par défaut. Non utilisés ici
	 */
	controlParameters->LIMINFC = new double [dimc];
	controlParameters->LIMSUPC = new double [dimc];
	for(unsigned long long int tabIndice=0;tabIndice<dimc;tabIndice++)
	{
		controlParameters->LIMINFC[tabIndice]=0.0;//valeurs par default
		controlParameters->LIMSUPC[tabIndice]=1.0;//valeurs par default
	}
	this->readTabData(&dataRoot, controlParameters->LIMSUPC, "CONTROL_MAX_VALUES", dimc);
	this->readTabData(&dataRoot, controlParameters->LIMINFC, "CONTROL_MIN_VALUES", dimc);

	controlParameters->NBPOINTSC=new unsigned long long int[dimc];
	for(unsigned long long int tabIndice=0;tabIndice<dimc;tabIndice++)
	{
		controlParameters->NBPOINTSC[tabIndice]=1;//valeurs par default
	}
	this->readTabData(&dataRoot, controlParameters->NBPOINTSC, "CONTROL_GRID_POINTS", dimc);

	int dimc_ty;
	if(dataRoot.find("CONTROL_TYCHASTIC_DIMENSION")!=dataRoot.not_found())
	{
		dimc_ty= dataRoot.get<int>("CONTROL_TYCHASTIC_DIMENSION", 1);
		controlParameters->DIM_TY=(unsigned long long int) dimc_ty;
		/*
		 * paramètres par défaut. Non utilisés ici
		 */
		controlParameters->LIMINF_TY = new double [dimc_ty];
		controlParameters->LIMSUP_TY = new double [dimc_ty];
		for(unsigned long long int tabIndice=0;tabIndice<dimc_ty;tabIndice++)
		{
			controlParameters->LIMINF_TY[tabIndice]=0.0;//valeurs par default
			controlParameters->LIMSUP_TY[tabIndice]=1.0;//valeurs par default
		}
		this->readTabData(&dataRoot, controlParameters->LIMSUP_TY, "CONTROL_TY_MAX_VALUES", dimc_ty);
		this->readTabData(&dataRoot, controlParameters->LIMINF_TY, "CONTROL_TY_MIN_VALUES", dimc_ty);

		controlParameters->NBPOINTS_TY=new unsigned long long int[dimc_ty];
		for(unsigned long long int tabIndice=0;tabIndice<dimc_ty;tabIndice++)
		{
			controlParameters->NBPOINTS_TY[tabIndice]=1;//valeurs par default
		}
		this->readTabData(&dataRoot, controlParameters->NBPOINTS_TY, "CONTROL_TY_GRID_POINTS", dimc_ty);
	}
	else
	{
		dimc_ty= 0;
		controlParameters->DIM_TY=0;
		/*
		 * paramètres par défaut. Non utilisés ici
		 */
		controlParameters->LIMINF_TY = new double [dimc_ty];
		controlParameters->LIMSUP_TY = new double [dimc_ty];
		controlParameters->NBPOINTS_TY=new unsigned long long int[dimc_ty];
	}

}

void ParametersManager::readTabData (ptree *dataRoot, unsigned long long int * target, string label, int nbElements )
{
	if(dataRoot->find(label)!=dataRoot->not_found())
	{
		ptree  tabTree=dataRoot->find(label)->second;
		ptree::const_iterator it = tabTree.begin();
		for(unsigned long long int tabIndice=0;tabIndice<nbElements;tabIndice++)
		{
			target[tabIndice]=(unsigned long long int)(*it).second.get_value<int>();
			it++;

			cout<< " lu  depuis "<<  label << " "<<target[tabIndice]<<endl;
		}
	}
}

void ParametersManager::readTabData (ptree *dataRoot, double * target, string label, int nbElements )
{
	if(dataRoot->find(label)!=dataRoot->not_found())
	{
		ptree  tabTree=dataRoot->find(label)->second;
		ptree::const_iterator it = tabTree.begin();
		for(unsigned long long int tabIndice=0;tabIndice<nbElements;tabIndice++)
		{
			target[tabIndice]=(*it).second.get_value<double>();
			it++;

			cout<< " lu  depuis "<<  label << " "<<target[tabIndice]<<endl;
		}
	}
}
ParametersManager::~ParametersManager() {
	// TODO Auto-generated destructor stub
}

