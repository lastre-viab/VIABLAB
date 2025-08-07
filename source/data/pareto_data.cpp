/*! \file  testZermelo.h
 *
 *
 *  \author: A. D�silles, LATSRE
 *  \brief  Fichier contant la d�finition d'un mod�le de viabilit� pour les tests
 *
 *
 *
 *
 */

#include "../include/utilities.h"

extern "C" {


double a=0.06;
double rho=0.05;
double d=4.0;
double b=0.08;

 double prodFunc(double K)
{
  double x=max(K,0.0);
  return  sqrt(K);
}
 double prodFuncPrime(double K)
{
  double x=max(K,0.0001);

  return  1.0/(sqrt(x)*2.0);
}
 double utility(double C)
{
  return   (sqrt(C) );
  // return  log(C);
}
 double utilityPrime(double C)
{
  return   (0.5/sqrt(max(C, 0.0001)) );
  // return  log(C);
}




void dynamics(const double *x, const double *u, double *image)
{
  double alpha=u[0];
  double beta=u[1];

  double K=x[0];
  double P=x[1];
  double z1=x[2];

  image[0]=  (1.0-alpha-beta)*prodFunc(K)-a*K;
  image[1]=  (1-beta*d)*prodFunc(K)-b*P;
  image[2]=  rho*x[2];//+utility(alpha*prodFunc(K));
  //cout<< " dynamique renvoie "<<image[0]<< " "<<image[1]<<endl;
}

 double l(const double *x, const double *u )
{
   return pow(x[1], 1.1)/1.1;
}


 double m(const double *x, const double *u )
{
   return rho;
}


/*!
 * jacobian matrix of the dynamics
 */
void jacobian(const double *x, const double *u , double ** jacob)
{
  double alpha=u[0];
  double beta=u[1];

  double K=x[0];
  double P=x[1];
  double z1=x[2];

  jacob[0][0]=prodFuncPrime(K)-a;
  jacob[0][1]=0.0;
  jacob[0][2]=0.0;

  jacob[1][0]=(1-beta*d)*prodFuncPrime(K);
  jacob[1][1]=-b;
  jacob[1][2]=0.0;

  jacob[2][0]=0.0;
  jacob[2][1]=0.0;
  jacob[2][2]=0.0;

}
/*
 * Function  that  defines  analyticaly computed local bounds of all componennts of the dynamic function
 * @param[in] x  the state variable
 * @param[out] res  the result
 */
void localDynBounds(const double *x, double * image)
{
  double K=x[0];
  double P=x[1];
  double z1=x[2];

  image[0]=  abs(prodFunc(K)-a*K);
  image[1]= max(abs( (1- d)*prodFunc(K)-b*P),abs( prodFunc(K)-b*P));
  image[2]= 0.0;
}


/*!
 * \brief Function  defining the mixed  constraints
 *
 * This function defines the set U(x) for admissible controls as function of the state
 * @param x state variable
 * @param u control variable
 * @return  value that caraterise the constraints set
 */
 double constraintsXU( const double *x, const double *u )
{
  double alpha=u[0];
  double K=x[0];
double res=max( -utility(alpha*prodFunc(K))-x[2],u[0]+u[1]-1.0);
  return (res<=0)?1.0:PLUS_INF;
}



 double constraintsX( const double *x )
{
	return 0.0;
}

}
