/*
 * initParams.h
 *
 *  Created on: 24 août 2018
 *      Author: anna_desilles
 */

#ifndef SRC_INITPARAMS_H_
#define SRC_INITPARAMS_H_

#include "../include/defs.h"
#include "../include/ParametersManager.h"

ParametersManager * initParams(gridParams &gp, algoViabiParams &avp, controlParams &cp,systemParams &sp);
ParametersManager * initParams(gridParams &gp, algoViabiParams &avp, controlParams &cp,systemParams &sp)
{

	string input_tempfile="../INPUT/"+ controlParamsFile;
	ptree dataRoot;
	read_json(input_tempfile, dataRoot);
	int dimc = dataRoot.get<int>("CONTROL_DIMENSION", 1);

	 initControls = new double[dimc*nbTrajs];
	 initControls[0]=0.0;
	 initControls[1]=0.0;


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

	sp.COMPUTE_LC=computeLC;
	sp.COMPUTE_MF=computeM;

	sp.LIP=LC;
	sp.MF=M;

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
		sp.COMPUTE_LC=0.0;
		sp.COMPUTE_MF=0.0;
		sp.LIP=1.0;
		sp.MF=1.0;
		break;
	}
	case 3:
	{
		sp.DYNAMICS=&dynamics;

		sp.COMPUTE_LC=0.0;
		sp.COMPUTE_MF=0.0;
		sp.LIP=1.0;
		sp.MF=1.0;
		break;
	}
	case 4:
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
	if(l_fd)
	{
		sp.L_FUNC_FD=&l_fd;
	}
	else
	{
		sp.L_FUNC_FD=&l_fd_default;
	}

	sp.SCHEME=discret_type;
	sp.SCALE_PARAM=false;
	for(int i=0;i<dim;i++)
	{
		sp.SCALE_PARAM|=scaling[i];
	}
	sp.TARGET=&target;
	if(target_fd)
	{
		sp.TARGET_FD=&target_fd;
	}
	else
	{
		sp.TARGET_FD=&target_fd_default;
	}
	if(constraintsXU_fd)
		sp.CONSTR_XU_fd = &constraintsXU_fd;
	else
		sp.CONSTR_XU_fd  =NULL;
	if(dynamics_fd)
		sp.DYNAMICS_FD = &dynamics_fd;
	else
		sp.DYNAMICS_FD = NULL;
	if(dynamics_tych_fd)
		sp.DYNAMICS_TYCH_FD = &dynamics_tych_fd;
	else
		sp.DYNAMICS_TYCH_FD = NULL;
	if(constraintsX_fd)
		sp.CONSTR_X_fd = &constraintsX_fd;
	else
		sp.CONSTR_X_fd = NULL;
	if(dynConstraintsForTraj)
		sp.DYN_CONSTR_FOR_TRAJ = &dynConstraintsForTraj;
	else
		sp.DYN_CONSTR_FOR_TRAJ = &dynConstraintsForTraj_default;

	if(constraintsXUY_fd)
			sp.MU_FUNC_FD = &constraintsXUY_fd;
		else
			sp.MU_FUNC_FD = &constraintsXUY_fd_default;
	/*
	 * Initialisation des paramètres de systèmes dynamique
	 * Ici toutes les valeurs sont par defaut, non utilisés
	 */
	avp.NB_OMP_THREADS=ompThreads;

	avp.FILE_PREFIX=prefix;
	avp.TARGET_OR_DEPARTURE=target_or_departure_problem;
	avp.COMPUTE_TMIN=compute_tmin;
	avp.NB_TRAJS = nbTrajs;
	avp.TYPE_TRAJ = typeTraj;
	avp.INIT_POINTS = initPoints;
	avp.INIT_POINTS_FD = initPoints_fd;
	avp.INIT_VALUES = initValues;
	avp.INIT_VALUES_FD = initValues_fd;

	avp.INIT_CONTROLS = initControls;
	avp.GRID_REFINMENTS_NUMBER = refine;
	avp.INTERMEDIATE_SAVINGS = intermediate_savings;
	avp.SAVE_BOUNDARY = saveBoundary;
	avp.SAVE_PROJECTION = saveProjection;
	avp.SAVE_SLICE = saveCoupe;
	avp.SAVE_SLICE_BOUND = saveCoupeBound;
	avp.PROJECTION = projection;
	avp.SAVE_SUBLEVEL = saveSubLevel;
	avp.LEVEL=level;


	if(avp.COMPUTE_TMIN)
	{
		sp.L_LIP=1.0;
		sp.L_MF=sp.maxTime;
	}
	else{
		sp.L_LIP=l_Lip;
		sp.L_MF=l_max;
	}
	return new ParametersManager(&gp, &avp, &cp, controlParamsFile, &sp);

}


#endif /* SRC_INITPARAMS_H_ */
