#include "../include/utilities.h"
#include "../include/ParametersManager.h"
#include "../include/ControlPickStrategy.h"

extern "C" {
//------------------------------------------------------------------------------------------------------ 
//  Model description header  for the problem SimplePopGrowth
//  Edit this file to define your model parameters 
//------------------------------------------------------------------------------------------------------ 
std::string paramsFile = "SimplePopGrowth_params.json";
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
image[0]=x[0]*x[1];
image[1]=u[0];
} 
//------------------------------------------------------------------------------------------------------ 
// Jacobian matrix of the dynamics. It will be used to estimate locally the Lipschitz constant  
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void jacobian(const double *x, const double *u, double ** jacob){
// programm here the jacobian matrix of your dynamics function: 
// x :     the input  state vector
// u :     the input control vector 
// jacob : the resulting matrix 
jacob[0][0]=0.0;
jacob[0][1]=0.0;
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
    
double controlWeight(const double *x, const double *u, double t, int, int) {
    if (x[0] > 1.5) {
        return -u[0];
    }
    else {
        return u[0];
    }
}

// Exemple : bulle de norme 3
bool isValidNeighbor(const unsigned long long int *coord1, const unsigned long long int *coord2, const double *bubbleRadius, int, int) {

    long long int diff1 = coord1[0] - coord2[0];
    double cube1 = std::fabs(diff1 * diff1 * diff1);

    long long int diff2 = coord1[1] - coord2[1];
    double cube2 = std::fabs(diff2 * diff2 * diff2);

    double d3 = bubbleRadius[0]*bubbleRadius[0]*bubbleRadius[0];
    
    return cube1 + cube2 <= d3 ;
}


class CustomStrategy final : public UserPickStrategy {
public:    
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
    ControlPickCriteria &criteria) override {

        double rho = criteria.getGridTimeStep();
        pickedControl p {
            rho,
            42                
        };
        
        return OptionalCu(p);
    }        
};

UserPickStrategy *newCustom(int strategyIndex, const TrajectoryParametersManager *) {
    return new CustomStrategy();
}

    
void temporalControl(double time, int, int, double *control) {
    if (time < 25.0) {
        control[0] = -0.5;
    }
    else if (time < 50.0) {
        control[0] = -0.17;
    }
    else if (time < 75.0) {
        control[0] = 0.17;
    }
    else  {
        control[0] = 0.5;
    }
}
}

