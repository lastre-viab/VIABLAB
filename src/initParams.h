/*
 * initParams.h
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
 *  Created on: 24 août 2018
 *      Author: anna_desilles
 */

#ifndef SRC_INITPARAMS_H_
#define SRC_INITPARAMS_H_

#include "../include/defs.h"

void initParams(gridParams &gp, algoViabiParams &avp, controlParams &cp,systemParams &sp);
void initParams(gridParams &gp, algoViabiParams &avp, controlParams &cp,systemParams &sp)
{
  gp.DIM=dim;
  gp.GRID_TYPE=0;
  gp.LIMINF=STATE_MIN;
  gp.LIMSUP=STATE_MAX;
if((computeSet==0)  & (refine >0))
  {
  for(int k=0;k<dim;k++)
    {
    for(int j=0;j<refine;j++)
      {
      nbPointsState[k]=2*nbPointsState[k]-1;
      }
    }
  }
  gp.NBPOINTS=nbPointsState;
  gp.PERIODIC=periodic;
  gp.FILE_PREFIX=prefix;
  gp.GRID_MAIN_DIR=dirTramage;
  gp.OMP_THREADS=ompThreads;

  /*
   * paramètres par défaut. Non utilisés ici
   */
  gp.SLICE_DIRS=sliceDirs;
  gp.SLICE_VALUES=sliceVals;
  gp.SORTIE_OK_INF=sortieOKinf;
  gp.SORTIE_OK_SUP=sortieOKsup;


  /*
   * initialisation des paramètres des contrôles
   */
  cp.DIMC=dimc;
  /*
   * paramètres par défaut. Non utilisés ici
   */
  cp.LIMINFC=CONTROL_MIN;
  cp.LIMSUPC=CONTROL_MAX;
  cp.NBPOINTSC=nbPointsControl;
  cp.DIM_TY=0;
  cp.LIMINF_TY=CONTROL_MIN_ty;
  cp.LIMSUP_TY=CONTROL_MAX_ty;
  cp.NBPOINTS_TY=nbPointsControl_ty;
  sp.COMPUTE_LC=computeLC;
  sp.COMPUTE_MF=computeM;
  sp.CONSTR_XU=&constraintsXU;
  sp.CONSTR_X=&constraintsX;
  /*
   * Initialisation des paramètres de systèmes dynamique
   */
  sp.DYN_TYPE=dynType;
  switch(dynType)
  {
    case 1:
      {
      sp.DYNAMICS=&dynamics;
      break;
      }
    case 2:
      {
      sp.DYNAMICS=&dynamics;
      break;
      }
    case 3:
      {
      sp.DYNAMICS=&dynamics_hybrid;
      break;
      }
  }
  /*
   * paramètres par défaut. Non utilisés ici
   */
  sp.globDeltat=globalDeltaT;
  sp.maxTime=T;
  sp.JACOBIAN=&jacobian;
  sp.LOCAL_DYN_BOUNDS=&localDynBounds;
  sp.M_FUNC=&m;
  sp.L_FUNC=&l;
  sp.LIP=LC;
  sp.MF=M;
  sp.SCHEME=discret_type;
  sp.SCALE_PARAM=false;
  for(int i=0;i<dim;i++)
    {
    sp.SCALE_PARAM|=scaling[i];
    }
  sp.TARGET=&target;


  /*
   * Initialisation des paramètres de systèmes dynamique
   * Ici toutes les valeurs sont par defaut, non utilisés
   */
  avp.NB_OMP_THREADS=ompThreads;

  avp.FILE_PREFIX=prefix;
  avp.TARGET_OR_DEPARTURE=target_or_departure_problem;
  avp.COMPUTE_TMIN=compute_tmin;

}


#endif /* SRC_INITPARAMS_H_ */
