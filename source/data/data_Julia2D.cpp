//------------------------------------------------------------------------------------------------------ 
//  Model description header  for the problem Julia2D
//  Edit this file to define your model parameters 
//------------------------------------------------------------------------------------------------------
#include "../include/utilities.h"

extern "C" {

std::string paramsFile = "Julia2D_params.json";

double R=2.0;      // Ball radius for computations
double a=-0.8;  // Dynamics parameters = u=(a,b)
double b= 0.156;
//------------------------------------------------------------------------------------------------------ 
// Dynamics definition  
// Please fill the code of the below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void dynamics(const double *x, const double *u, double *image)
{
  image[0]=     x[0]*x[0] - x[1]*x[1]+a;
  image[1]=   2*x[0]*x[1]+ b;
}

 double constraintsX( const double *x )
{
  double res1=  x[0]*x[0] + x[1]*x[1]-R*R;
  double res=(res1<=0)?1.0:PLUS_INF;
  return res;
}


}
