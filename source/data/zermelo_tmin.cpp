/*! \file  zermelo_tmin_data.h
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

#include "../include/utilities.h"

extern "C" {
  std::string paramsFile = "zermelo_tmin_params.json";
/*
 * Definition of the dynamics  and associated functions and constants
 */
double c=0.0, a=0.25;

void dynamics(const double *x, const double *u, double *image)
{
  image[0]= cos(u[0]) + c - a*tanh(x[1]);
  image[1]= sin(u[0]);
}

/*!
 * jacobian matrix of the dynamics
 */

void jacobian(const double *x, const double *u , double ** jacob)
{
  jacob[0][0]=0.0;
  jacob[0][1]=-a*(1.0-tanh(x[1])*tanh(x[1]));


  jacob[1][0]=0.0;
  jacob[1][1]=0.0;
}


/*!
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
void localDynBounds(const double *x, double * res)
{
  res[0]=1.0+c+a;
  res[1]=1.0;

}



/*      *****************************************
 *  Definition of constraints and target
 *************************************************** */

/*!
 * \brief Function  defining the state   constraints, corresponds  to k(x)
 *
 * This function defines the set K for admissible  states
 * @param x state variable
 * @return  value that characterize the constraints set
 */
int nbobs = 2; double obstacle[2*4]={-3.5, 1.0, 0.5, 1.0, -1.5, 0.0, 0.5, 0.75};

 double constraintsX( const double *x )
{
  double x1=x[0], x2=x[1];
   double res=-PLUS_INF, current;
   int i;
   for( i=0;i<nbobs;i++){
     current = -max(abs(x1-(obstacle[i*4]))-obstacle[i*4+2], abs(x2-(obstacle[i*4+1]))-obstacle[i*4+3]);
     res=max(res, current);
   }
   double resFinal=(res>0)?PLUS_INF:0.0;
   return resFinal;
}
/*!
 * Function that  characterise  the target set C, corresponds to c(x)
 * @param x state variable
 * @return value  to charaterise the target set
 */
  double target (const double *x)
{
        double res;
        if(sqrt(x[0]*x[0]+x[1]*x[1])<=0.2)
        {
                res=0.0;
        }
        else
        {
                res=PLUS_INF;
        }
        return res;
}

}
