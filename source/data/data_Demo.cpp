//------------------------------------------------------------------------------------------------------ 
//  Model description header  for the problem Demo
//  Edit this file to define your model parameters 
//------------------------------------------------------------------------------------------------------

#include "../include/utilities.h"

extern "C" {

std::string paramsFile = "Demo_params.json";
//------------------------------------------------------------------------------------------------------ 
// Dynamics definition  
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void dynamics(const double *x, const double *u, double *image)
{ 
// programm here your dynamics function: 
// x :     the input  state vector
// u :     the input control vector 
// image : the resulting state vector 
image[0]=x[0] - x[1];
image[1]=u[0];
} 
//------------------------------------------------------------------------------------------------------ 
// Jacobian matrix of the dynamics. It will be used to estimate locally the Lipschitz constant  
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void jacobian(const double *x, const double *u , double ** jacob){
// programm here the jacobian matrix of your dynamics function: 
// x :     the input  state vector
// u :     the input control vector 
// jacob : the resulting matrix 
jacob[0][0]=-1.0;
jacob[0][1]=1.0;
jacob[1][0]=0.0;
jacob[1][1]=0.0;
} 
//------------------------------------------------------------------------------------------------------ 
// Local bound of the dynamics. Each element of the result vector is defined as max(over u) F_i(x,u)  
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void localDynBounds(const double *x, double * bound){
// programm here the bound of your dynamics function, maximized over the control: 
// x :     the input  state vector
// bound : the resulting vector 
bound[0]=0.0;
bound[1]=0.0;
} 

}
