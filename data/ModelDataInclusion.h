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
 *  \author: A. Dï¿½silles, LATSRE
 *
 *  \brief  Ce fichier  sert à déclarer le modèle utilisateur et le fichier de paramètres d'exécution associé
 *
 *
 *
 *
 */


#ifndef MODELDATA_H_
#define MODELDATA_H_



string paramsFile = "zermelo_tmin_params.json";

//string paramsFile = "zermelo_Lmin_params.json";
//string paramsFile = "testPSP.json";
//string paramsFile = "testPSPBis.json";
//string paramsFile = "resilience.json";
//string paramsFile = "testPendule.json";
//string paramsFile = "allParams.json";
//string paramsFile = "agroEcoDiv_D4.json";


/************************************************************************************************
 * inclusion des dÃ©clarations du modÃ¨le
 */
//#include  "../data/ex1Viabi2D_data.h"   //--exemple 2D basique
#include  "../data/zermelo_tmin_new.h"   //-- Zermelo temps minimum
//#include  "../data/zermelo_Lmin_new.h"   //-- Zermelo critère intégral

//#include  "../data/testPendule_data.h"   //-- Bassin de capture, pendule
//#include  "../data/resilience_data.h"   //-- Noyeu de viabilité resilience
//#include "../data/PSP_dataBis.h"           // economie
//#include "../data/dataAgroEcoDivMultiParcelsD4_new.h" // agro ecodiv bi parcelle


#endif /* MODELDATA_H_ */
