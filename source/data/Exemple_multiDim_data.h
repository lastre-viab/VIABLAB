/*! \file  ex3bis_Viabi2D_data.h
 *
 *
  *
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
 *  \author: A. Désilles, LATSRE
 *  \brief  Fichier contant la définition d'un modèle de viabilité pour les tests
 *
 *  Ce fichier  contient les définitions  de tous les éléments qui décrivent un
 *  un problème de viabilité concret et  qui sont considérés comme
 *  paramètres par le code; Ce fichier représente ainsi une interface d'entrée primitive pour le code.
 *
 *  Les éléments principaux définis ici  sont :
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
 */
#ifndef EX3BIS_VIABI2D_DATA_H_
#define EX3BIS_VIABI2D_DATA_H_
/*
 * paramètres du modèle
 */
#include "defs.h"
double r=1.5;      // rayon de la contrainte
double alpha=1.0;  // rayon de contrôles
int stateDim = 2;

void loadModelData(ParametersManager *PM)
{
	gridParams * gp = PM->getGridParameters();
	stateDim = gp->DIM;
}

void dynamics(double * x, double *u, double * image)
{
  double m=0.0;
  for(int i=2;i<stateDim;i++)
    {
    m+=x[i];
    }
  m=m*m/((stateDim-2)*(stateDim-2));
  image[0]=   x[0]+m-x[1];
  image[1]=  u[0];
  for(int k=2;k<stateDim;k++)
    {
    image[k]=  0.0;
    }


}
/*!
 * jacobian matrix of the dynamics
 */
inline void jacobian(double *x, double *u , double ** jacob)
{
  double m=0.0;
  for(int i=2;i<stateDim;i++)
    {
    m+=x[i];
    }
  jacob[0][0]=1.0;
  jacob[0][1]=-1.0;

  for(int j=2;j<stateDim;j++)
    {
    jacob[0][j]=2.0*m/((stateDim-2)*(stateDim-2));
    }

  for(int i=1;i<stateDim;i++)
    {
    for(int j=0;j<stateDim;j++)
      {
      jacob[i][j]=0.0;
      }
    }
}
/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
inline void localDynBounds(double * x, double * res)
{
  double m=0.0;
  for(int i=2;i<stateDim;i++)
    {
    m+=x[i];
    }
  m=m*m/((stateDim-2)*(stateDim-2));
  res[0]=  abs(  m+x[0]-x[1]);
  res[1]=  1.0;
  for(int i=2;i<stateDim;i++)
    res[i]=  0.0;
}
#endif /* EX3BIS_VIABI2D_DATA_H_ */
