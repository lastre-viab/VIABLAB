
string paramsFile = "labyrinthe_params.json";

//------------------------------------------------------------------------------------------------------ 
//  Model description header  for the problem labyrinthe
//  Edit this file to define your model parameters 
//------------------------------------------------------------------------------------------------------ 
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
image[0]=u[0] * sin(u[1]);
image[1]=u[0] * cos(u[1]);
} 

double r = 8.0;
double eps2 = 0.25;
//------------------------------------------------------------------------------------------------------ 
// Function defining the target set. It should return a finite double value for any point  
// inside the target set, and the infinity for any point outside the target set  
// IMPORTANT : VIABLAB defines infinity as constant PLUS_INF. Please use it in this function  
// IMPORTANT : for target problems with epigraphic property, the values inside the target set are used as 
// initial state of teh value function 
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 
inline  double target (double * x){ 
// programm here the target function: 
// x :      the input  state vector
// return : PLUS_INF if x is not in the target set, otherwise, return a double value
double rho = x[0]*x[0] + x[1] * x[1];
if(sqrt(rho) >= r)
{
return 0.0;
}
else
{
return PLUS_INF;
}
} 

//------------------------------------------------------------------------------------------------------ 
// Function defining the constraints set K. It should return a finite double value for any point  
// inside K, and the infinity for any point outside K  
// IMPORTANT : VIABLAB defines infinity as constant PLUS_INF. Please use it in this function  
// IMPORTANT : for viability problems with epigraphic property, this function is used as 
// initial state of the value function 
// Please fill the code of teh below function without changing it's signature  
//------------------------------------------------------------------------------------------------------ 


inline double constraintsX( double * x )
{
// programm here your constraints function: 
// x : the input  state vector
// The function returns PLUS_INF for points outside the contrants set
// and a double value for points inside the constraints set

double rho = sqrt(x[0]*x[0] + x[1] * x[1]);
int a = 0;

bool outOfK = false;
while (!outOfK &&  2*a +2 <= r)
{
outOfK =  ((double)(2*a+1) < rho ) && (rho < (double)(2*a+2));
if(a%2 == 0)
{
outOfK &= x[0]*x[0] > eps2;
}
else
{
outOfK &= x[1]*x[1] > eps2;
}
a++;
}
 
 return outOfK ? PLUS_INF : 1.0;
}

