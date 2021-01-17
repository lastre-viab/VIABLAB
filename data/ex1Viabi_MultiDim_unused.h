/*! \file  ex3bis_Viabi2D_unused.h
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
 */


#ifndef EX3BIS_VIABI2D_UNUSED_H_
#define EX3BIS_VIABI2D_UNUSED_H_

const int dimc_ty=1;
double CONTROL_MIN_ty[dimc_ty]={0.0};
/*! \var CONTROL_MAX_ty[dimc]
 *  \brief maximum values  for  control vector components
 *   \see controlParams
 */
double CONTROL_MAX_ty[dimc_ty]={0.0};
unsigned long long int nbPointsControl_ty[dimc_ty] =  {2};

/*
 * sélection de type de reconstruction de trajectoire
 *  VD= 1, viable par défaut, on sélectionne le premier contrôle viable trouvé
 *  VL= 2, viable lourd: le contrôle reste constant tant qu'il viable;
 *  cette méthode nécessite une initialisation de contrôle
 */
int typeTraj=VD;

const int nbTrajs=0;
double initPoints[]={};
double initControls[]={};

const int CaseCtrlDependsOnX=false; //true : u d�pend de x : lancement du programme KhiXU

int saveProjection=0;
/*!
 * Function  defining the  switch conditions  between the continueos and  discrete dynamics
 * @param x state variable
 * @param u control variable
 * @return  value that caraterise the restet  set
 */
inline double resetSet( double * x, double * u ){
  return 1.0;
}
/*!
 * \var dynType  defines the type of dynamics in the model
 *      1 or CC : continuous
 *      2 or DD : discrete
 *      3 or HD : hybrid
 */
const int dynType=CC;
int compteOnlyOptimalRetro=1;

/*!
 * Function that defines the  dynamical system
 */
void dynamics_hybrid(double * x, double *u, double * image){
}
/*!
 * Function that  characterise  the target set C, corresponds to c(x)
 * @param x state variable
 * @return value  to charaterise the target set
 */
inline  double target (double * x){
  return PLUS_INF;
}

/* *   *************************************************************************
 *      Definition of value function
 ***************************************************************************/
 /*!
 * Function  for optimisation criterion, corresponds to l(x,u)
 * @param x state variable
 * @param u control
 * @return value of l
 */
inline double l(double * x, double * u ){
  return 1.0;
}
 /*!
 * Function  for optimisation criterion, corresponds to m(x, u)
 * @param x state variable
 * @param u control
 * @return value of m
 */
inline double m(double * x, double * u ){
  return 0.0;
}
 void postProcess();
void postProcess(){
}
int compute_tmin=0;


#endif /* EX3BIS_VIABI2D_H_ */
