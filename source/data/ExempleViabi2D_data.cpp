/*! \file  testZermelo.h
 *
 *
 *  \author: A. D�silles, LATSRE
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
 *      -# fonction l(t,x,u)
 *      -# fonction m(t,x,u)
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
/*
 * Some model specific parameters can be defined here
 */
double alpha=1.0;

void dynamics(const double *x, const double *u, double *image)
{

  image[0]=  alpha*x[0]*x[1];
  image[1]=  u[0];
  //cout<< " dynamique renvoie "<<image[0]<< " "<<image[1]<<endl;
}



/*!
 * jacobian matrix of the dynamics
 */
void jacobian(const double *x, const double *, double ** jacob)
{
  jacob[0][0]=alpha*x[1];
  jacob[0][1]=alpha*x[0] ;
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
  res[0]=  abs( alpha*x[0]*x[1]);
  res[1]=  0.8;
}

 double constraintsX( const double *x )
{
    double res1= std::max(0.2-x[0], x[0]-2.0);

    double res=(res1<=0) ? 1.0 : PLUS_INF;

  return res;
}

}
