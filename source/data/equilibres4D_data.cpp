/*! \file  Equilibres4D_data.h
 *
 * *
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

#include "../include/utilities.h"

extern "C" {

std::string paramsFile = "equilibres4D_data";
    
/*
 * Definition of the dynamics  and associated functions and constants
 */

void dynamics(const double *x, const double *u, double *image)
{
  image[0]=  0.0;
  image[1]=  0.0;
  image[2]=  0.0;
  image[3]= ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) );
}

/*!
 * jacobian matrix of the dynamics
 */

void jacobian(const double *x, const double *u , double ** jacob)
{
  jacob[0][0]=0.0;
  jacob[0][1]=0.0;
  jacob[0][2]=0.0;
  jacob[0][3]=0.0;

  jacob[1][0]=0.0;
  jacob[1][1]=0.0;
  jacob[1][2]=0.0;
  jacob[1][3]=0.0;

  jacob[2][0]=0.0;
  jacob[2][1]=0.0;
  jacob[2][2]=0.0;
  jacob[2][3]=0.0;

  jacob[3][0]= sign(x[0]-sin(x[2]) )*( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) )+
             ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( sign(x[0]+sin(x[2]) )  );

  jacob[3][1]=( sign( x[1]-cos(x[2]) ) )* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) )+
              ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( sign( x[1]+cos(x[2]) ) );

  jacob[3][2]=( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* (sign(x[0]+sin(x[2]) )*cos(x[2])- sign( x[1]+cos(x[2]) )*sin(x[2]) )+
              ( -sign(x[0]-sin(x[2]) )*cos(x[2]) + sign( x[1]-cos(x[2]) ) * sin(x[2]))* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) );

  jacob[3][3]=0.0;

}
/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
void localDynBounds(const double *x, double * res)
{
  res[0]=0.;
  res[1]=0.;
  res[2]=0.;
  res[3]= ( abs(x[0]-sin(x[2]) ) + abs( x[1]-cos(x[2]) ) )* ( abs(x[0]+sin(x[2]) ) + abs( x[1]+cos(x[2]) ) );

}

}
