#include "../include/utilities.h"
#include "../include/ParametersManager.h"
#include "../include/defs.h"
#include <cstdlib>

extern "C" {

//------------------------------------------------------------------------------------------------------
//  Model description header  for the problem timeMintest
//  Edit this file to define your model parameters 
//------------------------------------------------------------------------------------------------------ 
string paramsFile = "timeMintest_params.json";

void loadModelData(ParametersManager *PM)
{

}
//------------------------------------------------------------------------------------------------------ 
// Dynamics definition  
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void dynamics(double * x, double * u, double * image)
{ 
// programm here your dynamics function: 
// x :     the input  state vector
// u :     the input control vector 
// image : the resulting state vector 
image[0]=u[0];
image[1]= (x[0] <=1.0) ? x[0] : 2.0 - x[0];
} 
//------------------------------------------------------------------------------------------------------ 
// Jacobian matrix of the dynamics. It will be used to estimate locally the Lipschitz constant  
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void jacobian(double * x, double * u , double ** jacob){
// programm here the jacobian matrix of your dynamics function: 
// x :     the input  state vector
// u :     the input control vector 
// jacob : the resulting matrix 
jacob[0][0]=0.0;
jacob[0][1]=0.0;
jacob[1][0]=(x[0] <=1.0) ? 1.0 : -1.0;
jacob[1][1]=0.0;
} 
//------------------------------------------------------------------------------------------------------ 
// Local bound of the dynamics. Each element of the result vector is defined as max(over u) F_i(x,u)  
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void localDynBounds(double * x, double * bound){
// programm here the bound of your dynamics function, maximized over the control: 
// x :     the input  state vector
// bound : the resulting vector 
bound[0]=1.0;
bound[1]=0.0;
} 
//------------------------------------------------------------------------------------------------------ 
// Function defining the target set. It should return a finite double value for any point  
// inside the target set, and the infinity for any point outside the target set  
// IMPORTANT : VIABLAB defines infinity as constant PLUS_INF. Please use it in this function  
// IMPORTANT : for target problems with epigraphic property, the values inside the target set are used as 
// initial state of teh value function 
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
double target (double * x){
// programm here the target function: 
// x :      the input  state vector
// return : PLUS_INF if x is not in the target set, otherwise, return a double value
return sqrt( (x[0] - 2.0)* (x[0] - 2.0) + (x[1] - 3.0) * (x[1] - 3.0) ) <= 0.025 ? 0.0 : PLUS_INF;
} 
//------------------------------------------------------------------------------------------------------ 
// Function defining the constraints set K. It should return a finite double value for any point  
// inside K, and the infinity for any point outside K  
// IMPORTANT : VIABLAB defines infinity as constant PLUS_INF. Please use it in this function  
// IMPORTANT : for viability problems with epigraphic property, this function is used as 
// initial state of the value function 
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
double constraintsX( double * x )
{
// programm here your constraints function:
// x : the input  state vector
// The function returns PLUS_INF for points outside the contrants set
// and a double value for points inside the constraints set

double res = PLUS_INF;
double test= max(x[0] - 2.0, 0.0 - x[0]);
test = max(test, max(x[1] - 3.0, 0.0 - x[1]));
res = test <= 0.0 ? 0.0 : PLUS_INF;
return res;
}
}
