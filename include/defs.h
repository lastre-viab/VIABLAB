/*! \file  defs.h
 *
 *
 *  \author: A. D�silles, LATSRE
 *  \brief  Fichier  contenant les d�finitions de structures globales et les ioncludes std
 *
 *  Ce fichier  contient les d�finition de quelques structures
 *   qui servent principalement  au passage de param�tres   aux constructeurs de classes
 *
 *   Il contient �galement en un seul endroit toutes les inclusions
 *   de classes std ou autres utilis�es dans le projet
 *
 *   Enfin  il d�finit �galement quelques fonctions g�n�riques, hors classes
 *    qui peuvent �tre utiles globalement dans le projet
 */

#ifndef DEFS_H_

#define DEFS_H_

//#include "CropData.h"
//#include "CycleData.h"

#include <stdlib.h>
#include <vector>
#include <strstream>
#include <valarray>
#include <algorithm>
#include <sstream>


//#include <boost/thread.hpp>
#include <string.h>

#include <math.h>
#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sys/types.h>
#include "boost/dynamic_bitset/dynamic_bitset.hpp"


#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/detail/info_parser_error.hpp"
#include "boost/property_tree/detail/info_parser_utils.hpp"
#include "boost/property_tree/json_parser.hpp"
using namespace boost::property_tree;
/*
#if defined(__APPLE__)

#include <sys/param.h>
#include <mach-o/dyld.h>
#endif
#include <inttypes.h>
 */

//#include <magick++.h>
#include <list>
#include <omp.h>
#include <ctime>
#include<sys/time.h>
#include<sys/types.h>

#include<valarray>
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif


#define PLUS_INF 1e+15
#define NB_MAX_TRAJ_ITER 100000000
#define numZERO  0.0000000001
#define  EE   1
#define  EI   2
#define  RK2I 3
#define  RK2E 4

#define  RK4I 5
#define  RK4E 6


#define  FUNC   1
#define  RETRO   2



#define DEV_PRINT 1


/*!
 * valeurs pour la constante  qui definit le type de dynamique
 */
#define CC  1  /*! continuous dynamics*/
#define DC  2  /*! discrete time dynamics*/
#define HD  4  /*! hybrid dynamics*/
#define DD  3  /*! full discrete ( time and space) dynamics*/
/*
 * Valeurs pour la constante qui définit la méthode de
 * représentation de l'ensemble
 */
#define BS 1
#define MM 2

#define pi  3.14159265359
/*!
 * valeurs pour la constante  qui definit le type de probleme
 */


/*!
 * valeurs pour la constante  de  type  de probleme
 */

#define VIAB 1
#define CAPT 2

/*!
 * valeurs pour la constante  de  type  de trajectoire
 */
#define VD 1 // viable par defaut
#define VL 2 // viable lourd
#define OP 3 // optimale

using namespace std;
//using namespace Magick;





#define GRID_DB_FILE "gridDB.db"

#define RETRO_OPTI_FILE "retro_opti.db"

#define  RETRO_VIAB_FILE "retro_viabi.db"
#define  RETRO_DYN_FILE "retroDynDB.db"
#define  RETRO_DYN_TEXT_FILE "retroDynFile.dat"

/*!
 *  \struct gridParams
 *  \brief Grid parameters container
 */
struct gridParams
{
	unsigned long long int  DIM;       //!< dimension de la  grille
	unsigned long long int  *NBPOINTS;  //!<   tableau indicant le nombre de points par axe
	double    *LIMINF;    //!<   tableau  de limites inf pour chaque variable
	double    *LIMSUP;    //!<   tableau  de limites sup pour chaque variable
	int       *PERIODIC;  //!<   tableau  de boolens pour indiquer   les variables  qui sont p�riodiques (1=p�riodique, 0=non)
	int        GRID_TYPE; //!<  0=les points de calcul sont les noeuds de la grille, 1= les points de calcul sont les milieux des mailles
	string FILE_PREFIX;
	int GRID_MAIN_DIR;        //!< for bitSet grids only :  la grille sera stock�e comme tableau de segments paralleles  a cette dimension
	int * SLICE_DIRS;
	double *SLICE_VALUES;
	int * SORTIE_OK_INF;
	int * SORTIE_OK_SUP;
	int OMP_THREADS;
	int GRID_METHOD;

};

/*!
 *  \struct controlParams
 *  \brief controle parameters container
 *  \see SysDyn
 */

struct controlParams
{
	unsigned long long int DIMC;        //!<  dimension  de controles
	unsigned long long int *NBPOINTSC;  //!<  nombre de points par axe en cas de discr�tisation
	double    *LIMINFC;      //!<  limites inf de chaque variable de controle
	double    *LIMSUPC;       //!<  limites sup de chaque variable de controle

	unsigned long long int DIM_TY;        //!<  dimension  de controles
	unsigned long long int *NBPOINTS_TY;  //!<  nombre de points par axe en cas de discr�tisation
	double    *LIMINF_TY;      //!<  limites inf de chaque variable de controle
	double    *LIMSUP_TY;       //!<  limites sup de chaque variable de controle

};

/*!
 * \struct systemParams
 * \brief  dynamical system functions cntainer
 * \see SysDyn
 */

struct systemParams
{
	void (*DYNAMICS)(double *, double *, double *);
	void (*DYNAMICS_FD)(unsigned long long int *, unsigned long long int *, unsigned long long int *);
	void (*DYNAMICS_TYCH_FD)(unsigned long long int *, unsigned long long int *, unsigned long long int *, unsigned long long int *);

	double (*CONSTR_XU)(double *, double *u);
	double (*CONSTR_XU_fd)(unsigned long long int *, unsigned long long int*);
	double (*CONSTR_X)(double *);
	double (*CONSTR_X_fd) (unsigned long long int *);
	double (*DYN_CONSTR_FOR_TRAJ)(double *, double *);
	double (*TARGET) (double *);
	double (*TARGET_FD) (unsigned long long int *);
	double (*L_FUNC_FD)(unsigned long long int *, unsigned long long int*);
	double (*L_FUNC_TYCH_FD)(unsigned long long int *, unsigned long long int*, unsigned long long int*);

	double (*L_FUNC)(double *, double *);
	double (*MU_FUNC_FD)(unsigned long long int *, unsigned long long int*);
	double (*M_FUNC)(double *, double *);
	void (*JACOBIAN)(double *x, double *u , double ** jacob);
	void (*LOCAL_DYN_BOUNDS)(double * x, double * res);
	int COMPUTE_LC;
	double LIP;
	double L_LIP;
	bool globDeltat;
	int COMPUTE_MF;
	double MF;
	double L_MF;
	double maxTime;
	int SCHEME;
	bool SCALE_PARAM;
	int DYN_TYPE;     // type of dynamics : FD, CD,
	int FD_DYN_TYPE;  // pour indiquer le type de repr�sentation de la dynamique : fonction ou fihciers r�tro
	std::string RETRO_FILE_NAME;
};

struct algoViabiParams
{

	int SET_TYPE;
	int COMPUTE_SET;
	int NB_OMP_THREADS;
	int COMPUTE_TMIN;

	string FILE_PREFIX;
	int TARGET_OR_DEPARTURE;
	int NB_TRAJS;
	int TYPE_TRAJ;
	double * INIT_POINTS;
	double * INIT_VALUES;
	double * INIT_VALUES_FD;
	unsigned long long int * INIT_POINTS_FD;
	double * INIT_CONTROLS;
	int GRID_REFINMENTS_NUMBER;
	int INTERMEDIATE_SAVINGS;
	int SAVE_SLICE;
	int SAVE_SLICE_BOUND;
	int SAVE_BOUNDARY;
	int SAVE_PROJECTION;
	int SAVE_SUBLEVEL;
	int LEVEL;
	unsigned long long int* PROJECTION;


};

struct algoViabiBoolParams
{
	unsigned long long int NB_GO_RAM;

};


/*!
 * \struct discretImageSet
 * \brief  Structure pour le stockage d'une image discrete d'un point;
 *
 * Decrit sous forme compacte l'ensemble image d'un point de grille par la dynamique discrete
 *  \var nbImageCells repr�sente le nombre de mailles de la grille, hors doublons, contenues dans l'image
 *  \var tabImageCells est une tableau d'entiers  de taille nbImageCells  repr�sentant les num�ros de mailles contenues
 *  dans l'image
 *  \var tabImageControls repr�sente une permutation de l'ensemble
 *  d'entiers num�ros des controles {0,... nbTotalControls-1} dans l'ordre d'association avec les mailles images
 *  \var tabCellsEntrees donne le debut dans tabImageControls de la  liste des controles
 *  associes  aux  cellules
 * \see viabiHJB
 */
struct discretImageSet
{
	unsigned long long int nbImageCells;
	unsigned long long int * tabImageCells;
	unsigned long long int * tabImageControls;
	unsigned long long int * tabCellEntrees;
};

struct discretImageSet_DD
{
	unsigned long long int nbImagePoints;
	unsigned long long int * tabImagePoints;
	unsigned long long int * tabImageControls;
	unsigned long long int * tabPointsEntrees;
};

struct discretImageSet_simple
{
	unsigned long long  int nbImageCells;
	unsigned long long  int * tabImageCells;

};



/*!
 * \struct triple
 *
 * \brief  d�finit un triplet pour un ant�c�dant de cellule
 * \see ViabiHJB
 */
struct triple
{
	unsigned long long int posX; /*!< num�ro de l'�tat ant�c�dant */
	unsigned long long int posU; /*!< num�roe de contr�le tel que x-rho*f(x,u) est dans la cellule*/
	double rho;        /*!< le pas de temps */
};


inline bool equal_triple(triple t1, triple t2)
{
	return ((t1.posX==t2.posX) & (t1.posU==t2.posU) & (t1.rho==t2.rho));
}
inline bool compare_triple_first(triple t1, triple t2)
{
	return ((t1.posX<t2.posX));
}
inline bool compare_triple_second(triple t1, triple t2)
{
	return ((t1.posU<t2.posU));
}
/*!
 * \struct imageCell
 *
 * \brief  d�finit une structure pour une cellule image
 */
struct imageCell
{
	unsigned long long int cellNum;    /*!< numero de cellule */
	double minVal;           /*!< la valeur optimale de la fonction valeur sur cette cellule*/
	list<triple> R_optimal;  /*!< la liste des ant�c�dents minimisant la fonction valeur*/
	list<triple> R_viable;   /*!< la liste des ant�c�dents viables arrivant dans cette  cellule */
};

/*!
 * \struct imageCellsList structure utilis�e pour stocker de fa�on compacte
 *  l'image  par la dynamique discrete  d'une couche en cours.
 *
 *  Les valeurs minNum et maxNum servent � gagner du temps pour les algorithmes d'insetion.
 */
struct imageCellsList
{
	list<imageCell> cellsList;
	int minNum;
	int maxNum;

};

struct imageCellsIntList
{
	list<unsigned long long int> cellsList;
	unsigned long long int minNum;
	unsigned long long int maxNum;
};


/*!
 * \struct imageCell
 *
 * \brief  d�finit une structure pour une cellule image
 */
struct imagePoint
{
	unsigned long long int PointNum;    /*!< numero de cellule */
	double minVal;           /*!< la valeur optimale de la fonction valeur sur cette cellule*/
	list<triple> R_optimal;  /*!< la liste des ant�c�dents minimisant la fonction valeur*/
	list<triple> R_viable;   /*!< la liste des ant�c�dents viables arrivant dans cette  cellule */
};
struct imageTempPoint
{
	unsigned long long int PointNum;    /*!< numero de cellule */
	double minVal;           /*!< la valeur optimale de la fonction valeur sur cette cellule*/
	list<triple>* R_optimal;  /*!< la liste des ant�c�dents minimisant la fonction valeur*/
	list<triple>* R_viable;   /*!< la liste des ant�c�dents viables arrivant dans cette  cellule */
};

/*!
 * \struct imageCellsList structure utilis�e pour stocker de fa�on compacte
 *  l'image  par la dynamique discrete  d'une couche en cours.
 *
 *  Les valeurs minNum et maxNum servent � gagner du temps pour les algorithmes d'insetion.
 */
struct imagePointsList
{
	std::list<imagePoint> *pointsList;
	unsigned long long int minNum;
	unsigned long long int maxNum;
};

struct imagePointsIntList
{
	list<unsigned long long int> pointsList;
	unsigned long long int minNum;
	unsigned long long int maxNum;
};


/*!
 * op�rateur de comparaison des ellules images  utilis� pour
 * trier les listes  de cellules  images par le num�ro de cellule (le premier entier)
 *
 */
inline bool imageCellCompare(imageCell c1, imageCell c2)
{
	return (c1.cellNum<=c2.cellNum);
}

/***********************************************************************************************
 *  Fonctions abstraites g�n�riques utiles � toutes les classes
 ***********************************************************************************************/


/*! \brief Transformation  de num�rotation inverse
 *  � partir du num�ro du point dans la num�rotation alpha-num�rique
 * on obient le vecteur des coordonn�s enti�res
 */
inline void numToIntCoords_gen(unsigned long long int num, unsigned long long int dim, unsigned long long int * nbPoints,unsigned long long int *res)
{
	int temp=num;

	for(  int d=dim-1;d>=0;d--)
	{
		res[d]=temp%nbPoints[d];                       // coordonn�es enti�res du point
		temp=(temp-res[d])/nbPoints[d];

	}

}

inline void numToIntDoubleCoords_gen(unsigned long long int num, unsigned long long int dim, unsigned long long int * nbPoints,unsigned long long int *res)
{
	int temp=num;

	for(  int d=dim-1;d>=0;d--)
	{
		res[d]=temp%nbPoints[d];                       // coordonn�es enti�res du point
		temp=(temp-res[d])/nbPoints[d];
	}
}

/*! \brief Transformation  de num�rotation alpha-num�rique
 *  *  � partir du num�ro des coordonn�es enti�res du points dans la grille
 *  sont num�ro
 */
inline void intCoordsToNum_gen(unsigned long long int dim, unsigned long long int * nbPoints, unsigned long long int * coords, unsigned long long int * res)
{
	(*res)=coords[0];

	for(unsigned long long int i=0;i<dim-1;i++)
	{
		(*res)=(*res)*nbPoints[i+1]+coords[i+1];

	}
}


/*!
 *  le type intPair est d�fini pour le calcul des images
 *   discretes d'un point
 */
typedef  pair<int,int> intPair;

/*!
 * op�rateur de comparaison des paires d'entiers utilis� pour
 * trier les listes  de cellules par le num�ro de cellule (le premier entier)
 * le second entier est le num�ro de contr�le  correspondant
 */
inline bool pairCompare(intPair p1, intPair p2)
{
	return (p1.first<=p2.first);
}


typedef  pair<unsigned long long int,unsigned long long int> uintPair;

/*!
 * op�rateur de comparaison des paires d'entiers utilis� pour
 * trier les listes  de cellules par le num�ro de cellule (le premier entier)
 * le second entier est le num�ro de contr�le  correspondant
 */
inline bool upairCompare(uintPair p1, uintPair p2)
{
	return (p1.first<=p2.first);
}

inline void printVector(double *vect, unsigned long long int dim)
{
	for(unsigned long long int k=0;k<dim;k++)
	{
		cout<< " "<<vect[k];
	}
	cout<<endl;
}


inline void printVector(int *vect, unsigned long long int dim)
{
	for(unsigned long long int k=0;k<dim;k++)
	{
		cout<< " "<<vect[k];
	}
	cout<<endl;
}
inline void printVector(unsigned long long int *vect, unsigned long long int dim)
{
	for(unsigned long long int k=0;k<dim;k++)
	{
		cout<< " "<<vect[k];
	}
	cout<<endl;
}

inline double sign(double x)
{
	return (double) ((x < 0) ? -1 : (x > 0));
}

struct imageDDPoint
{
	unsigned long long int PointNum;    /*!< numero de cellule */
	double val;           /*!< la valeur optimale de la fonction valeur sur cette cellule*/
	list<intPair> dynFD;  /*!< la liste des ant�c�dents minimisant la fonction valeur*/
	list<intPair> R_viable;   /*!< la liste des ant�c�dents viables arrivant dans cette  cellule */
	list<intPair> R_opti;  /*!< la liste des ant�c�dents viables arrivant dans cette  cellule */
};




/*
 * Fonctions de changement d'échelle universel
 */
inline void  fScale( const double *x, double *res, double *STATE_MINreel, double *STATE_MAXreel, int *scaling, int dim)
{
	/*
	 * le changement de variables est efectif seulement si scaling[i]=1
	 */
	for(int i=0;i<dim;i++)
	{
		res[i]=scaling[i]*((x[i]-STATE_MINreel[i])/(STATE_MAXreel[i]-STATE_MINreel[i]))+(1-scaling[i])*x[i];
	}
}

inline void  fScaleInv( const double *x, double *res, double *STATE_MINreel, double *STATE_MAXreel, int *scaling, int dim)
{
	/*
	 * le changement de variables est efectif seulement si scaling[i]=1
	 */
	for(int i=0;i<dim;i++)
	{
		res[i]=scaling[i]*(x[i]*(STATE_MAXreel[i]-STATE_MINreel[i])+STATE_MINreel[i])+(1-scaling[i])*x[i];
	}
}

inline double fTimeScale( const double  x, double scaleL, int timeScaling )
{
	/*
	 * le changement de variables est efectif seulement si scaling[i]=1
	 */

	return (timeScaling *((x)*scaleL)+(1-timeScaling)*x);

}

inline double  fTimeScaleInv( const double  x, double scaleL, int timeScaling)
{
	return (timeScaling *((x)/scaleL)+(1-timeScaling)*x);
}

inline void  fScalePrime(  double *res, double *STATE_MINreel, double *STATE_MAXreel, int *scaling, int dim)
{
	/*
	 * le changement de variables est efectif seulement si scaling[i]=1
	 */
	for(int i=0;i<dim;i++)
	{
		res[i]=scaling[i]/(STATE_MAXreel[i]-STATE_MINreel[i])+(1-scaling[i]);
	}
}

inline double  fScaleComponent( const double x, int i, double *STATE_MINreel, double *STATE_MAXreel, int *scaling)
{
	/*
	 * le changement de variables est efectif seulement si scaling[i]=1
	 */
	return   scaling[i]*(x -STATE_MINreel[i])/(STATE_MAXreel[i]-STATE_MINreel[i])+(1-scaling[i])*x ;

}

inline double  fScaleInvComponent( const double x, int i, double *STATE_MINreel, double *STATE_MAXreel, int *scaling)
{
	/*
	 * le changement de variables est efectif seulement si scaling[i]=1
	 */

	return  scaling[i]*((x*(STATE_MAXreel[i]-STATE_MINreel[i])+STATE_MINreel[i]))+(1-scaling[i])*x;

}
/*
 * le changement de variables est efectif seulement si scaling[i]=1
 */
inline double  fScalePrimeComponent(  int i, double *STATE_MINreel, double *STATE_MAXreel, int *scaling)
{
	return scaling[i]/(STATE_MAXreel[i]-STATE_MINreel[i])+(1-scaling[i]);
}





#endif /* DEFS_H_ */
