//------------------------------------------------------------------------------------------------------ 
//  Model description header  for the problem LacGaranti
//  Edit this file to define your model parameters 
//------------------------------------------------------------------------------------------------------ 
string paramsFile = "LacGaranti_params.json";

double mm = 26.90;
double r = 101.96;

double Pmax;
double Lmax;
double umax;
//------------------------------------------------------------------------------------------------------ 
// Dynamics definition  
// Please fill the code of the below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
void dynamics_tych(double  * x, double * u, double * v, double * image){ 
// programm here your dynamics function: 
// x :      the input  state vector
// u :      the input control vector 
// v :      the input tyche vector 
// image :  the resulting state vector 

double b = v[0];
double alpha = v[1];
double q = v[2];
double lam = v[3];
double pq = x[1] > 0 ? pow(x[1], q) : 0.0;
double mq = pow(mm, q);
double p = x[1];
double e = exp(-lam*(p-mm));
image[0]= u[0];
image[1]= - b * p + x[0] + r * ((1- alpha) * pq / (mq + pq) + alpha * p / (p + mm * e));
} 
//------------------------------------------------------------------------------------------------------ 
// Jacobian matrix of the dynamics. It will be used to estimate locally the Lipschitz constant  
// Please fill the code of the below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
inline void jacobian_tych(double * x, double * u , double * v , double ** jacob){
// programm here the jacobian matrix of your dynamics function: 
// x :     the input  state vector
// u :     the input control vector 
// v :     the input tyche vector 
// jacob : the resulting matrix 
double q = v[2];
double b = v[0];
double alpha = v[1];
double lam = v[3];
double pq = x[1] > 0 ? pow(x[1], q) : 0.0;
double mq = pow(mm, q);
double p = x[1];
double e = exp(-lam*(p-mm));
double pqPrime = x[1] > 0 ? pow(x[1], q - 1.0) : 0.0;
jacob[0][0]=0.0;
jacob[0][1]=0.0;
jacob[1][0]=1.0;
jacob[1][1]= -b + r * ( (1- alpha) * mq * pqPrime / ((mq + pq) * (mq + pq)) + alpha * e *  mm * (1 + p * lam) / ((p + mm * e) * (p + mm * e) ) );
} 
//------------------------------------------------------------------------------------------------------ 
// Local bound of the dynamics. Each element of the result vector is defined as max(over u) F_i(x,u)  
// Please fill the code of the below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
inline void localDynBounds(double * x, double * bound){
// programm here the bound of your dynamics function, maximized over the control: 
// x :     the input  state vector
// bound : the resulting vector 
bound[0]=umax;
bound[1]=2.5 * Pmax + Lmax + r;
} 
//------------------------------------------------------------------------------------------------------ 
// Function defining the constraints set K. It should return a finite double value for any point  
// inside K, and the infinity for any point outside K  
// IMPORTANT : VIABLAB defines infinity as constant PLUS_INF. Please use it in this function  
// IMPORTANT : for viability problems with epigraphic property, this function is used as 
// initial state of the value function 
// Please fill the code of the below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
inline double constraintsX( double * x )
{
// programm here your constraints function: 
// x : the input  state vector
// The function returns PLUS_INF for points outside the contrants set
// and a double value for points inside the constraints set
     return 1.0;
}

inline double constraintsXV_tych( double * x , double *v)
{
// programm here your constraints function:
// x : the input  state vector
// The function returns PLUS_INF for points outside the contrants set
// and a double value for points inside the constraints set
     return 1.0;
}
  void loadModelData(ParametersManager *PM)
{
  gridParams * gp = PM->getGridParameters();
  controlParams *cp = PM->getControlParameters();
  Lmax = gp->LIMSUP[0];
  Pmax = gp->LIMSUP[1];
  umax = cp->LIMSUPC[0];
 }


