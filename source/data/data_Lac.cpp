#include "../include/utilities.h"
#include "../include/ParametersManager.h"

extern "C" {

//------------------------------------------------------------------------------------------------------ 
//  Model description header  for the problem Lac
//  Edit this file to define your model parameters 
//------------------------------------------------------------------------------------------------------ 
std::string paramsFile = "Lac_params.json";
//------------------------------------------------------------------------------------------------------ 
// Dynamics definition  
// Please fill the code of the below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 

double b = 0.8;
double q = 8.0;
double mm = 1.0;
double r = 1.0;
double mq = pow(mm,q);
double Pmax;
double Lmax;
double umax;

void dynamics(const double *x, const double *u, double *image)
{ 
// programm here your dynamics function: 
// x :     the input  state vector
// u :     the input control vector 
// image : the resulting state vector 
double pq = x[1] > 0 ? pow(x[1], q) : 0.0;
image[0]= u[0];
image[1]= - b * x[1] + x[0] + r * pq / (mq + pq);
} 
//------------------------------------------------------------------------------------------------------ 
// Jacobian matrix of the dynamics. It will be used to estimate locally the Lipschitz constant  
// Please fill the code of the below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void jacobian(const double *x, const double *u , double ** jacob){
// programm here the jacobian matrix of your dynamics function: 
// x :     the input  state vector
// u :     the input control vector 
// jacob : the resulting matrix 
double pq = x[1] > 0 ? pow(x[1], q) : 0.0;
double pq1 = x[1] > 0 ? pow(x[1], q - 1.0) : 0.0;
jacob[0][0]=0.0;
jacob[0][1]=0.0;
jacob[1][0]=1.0;
jacob[1][1]=-b + mq * pq1 / ((mq + pq) * (mq + pq));
} 
//------------------------------------------------------------------------------------------------------ 
// Local bound of the dynamics. Each element of the result vector is defined as max(over u) F_i(x,u)  
// Please fill the code of the below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void localDynBounds(const double *x, double * bound){
// programm here the bound of your dynamics function, maximized over the control: 
// x :     the input  state vector
// bound : the resulting vector 
bound[0]=umax;
bound[1]=b * Pmax + Lmax + 1.0;
} 
//------------------------------------------------------------------------------------------------------ 
// This function can be used to initialize some model specific parameter in runtime mode   
// For example, the state dimention can be a parameter of the model, defined n the Json file.   
// In that case, it can be retrieved from the paramerManager object PM that is the argumentof this function  
// Please fill the code of the below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void loadModelData(const ParametersManager *PM)
{
  const gridParams * gp = PM->getGridParameters();
  const controlParams *cp = PM->getControlParameters();
  Lmax = gp->LIMSUP[0];
  Pmax = gp->LIMSUP[1];
  umax = cp->LIMSUPC[0];
}

}
