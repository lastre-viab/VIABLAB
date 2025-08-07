 /*! \file  ModelDataInclusion.h
 *
 *  *
 *    VIABLAB : a numerical library for Mathematical Viability Computations
 *    Copyright (C) <2020>  <Anna DESILLES, LASTRE>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *   
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author: A. D�silles, LATSRE
 *
 *  \brief  Ce fichier  sert � d�clarer le mod�le utilisateur et le fichier de param�tres d'ex�cution associ�
 *
 *
 *
 *
 */


#ifndef MODELDATA_H_
#define MODELDATA_H_



//string paramsFile = "pareto_params.json";
//string paramsFile = "testPendule.json";
//string paramsFile = "Equilibres4D_params.json";
//string paramsFile = "ExempleMultiDim_params.json";
//string paramsFile = "JuliaSets_params.json";
//string paramsFile = "zermelo_tmin_params.json";
//string paramsFile = "zermelo_Lmin_params.json";
//string paramsFile = "ExempleViabi2D_params.json";
//string paramsFile = "Lac_params.json";
//string paramsFile = "LotkaVolterra_params.json";
string paramsFile = "VALIUM.json";
//string paramsFile = "Magique.json"
//string paramsFile = "allParams.json";

/************************************************************************************************
 * inclusion des déclarations du modèle
 */
//#include  "../data/pareto_data.h"   //--exemple 2D basique
//#include  "../data/Exemple_multiDim_data.h"   //--exemple 2D basique
//#include  "../data/zermelo_tmin.h" // Zermelo temps minimum
//#include  "../data/zermelo_Lmin.h"   //-- Zermelo crit�re int�gral
//#include  "../data/Julia2D_data.h" // Zermelo temps minimum
//#include "../data/data_Lac.h"
//#include "../data/LotkaVolterra.h"
//#include "../data/data_Magique.h"
#include "../data/VALIUM.h"
//#include  "../data/testPendule_data.h"   //-- Bassin de capture, pendule
//#include  "../data/resilience_data.h"   //-- Noyeu de viabilit� resilience
//#include "../data/equilibres4D_data.h"           // economie


#endif /* MODELDATA_H_ */
