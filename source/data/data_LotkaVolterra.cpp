#include "../include/utilities.h"

extern "C" {

std::string paramsFile = "LotkaVolterra_params.json";

double _r = 1.0;

double _m = 1.0;

double _umax = 0.5;

double seuil = 1.2;

void dynamics(const double *x, const double *u, double *image)

{

// programm here your dynamics function:

// x : the input state vector

// u : the input control vector

// image : the resulting state vector

image[0] = _r*x[0]-x[0]*x[1];

image[1] = -_m*x[1]+x[0]*x[1]-u[0]*_umax*x[1];

}

//------------------------------------------------------------------------------------------------------

// Jacobian matrix of the dynamics. It will be used to estimate locally the Lipschitz constant

// Please fill the code of teh below function without changing it's signature

//------------------------------------------------------------------------------------------------------

void jacobian(const double *x, const double *u , double ** jacob){

// programm here the jacobian matrix of your dynamics function:

// x : the input state vector

// u : the input control vector

// jacob : the resulting matrix

jacob[0][0]=_r-x[1];

jacob[0][1]=-x[0];

jacob[1][0]=x[1];

jacob[1][1]=-_m+x[0]-u[0]*_umax;

}

//------------------------------------------------------------------------------------------------------

// Local bound of the dynamics. Each element of the result vector is defined as max(over u) F_i(x,u)

// Please fill the code of teh below function without changing it's signature

//------------------------------------------------------------------------------------------------------

void localDynBounds(const double *x, double * bound){

// programm here the bound of your dynamics function, maximized over the control:

// x : the input state vector

// bound : the resulting vector

bound[0]=5.0;

bound[1]=5.0;

}

 double constraintsX( const double *x )
{
   double res=(x[0]>seuil)?1.0:PLUS_INF;

  return res;
}
}
