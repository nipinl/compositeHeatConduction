#ifndef SHELL_H_
#define SHEL_H_
#include<iostream>
using namespace std;
class Shell{
protected:
//double x;
//double y;
bool axi;//axi symmetric case; like tube
int M;//Number of divisions in y direction
int N;//Number of divisions in x direction
double ri;//Inner radius of the tube
double Width;//thickness of the tube
double Length;//Length of the cylinder
//Material properties
double tCond; 
double spHeat;
double density;
double simTime;//in seconds

//Initial condition
double initTemp;//Initial uniform temperature of the tube

// temp bc
double TLeft,TRight,TBottom,TTop;

// heat flux bc
double qLeft,qRight,qBottom,qTop;

//convection
double TfLeft,TfRight,TfBottom,TfTop{;
double hfL,hfR,hfB,hfT;

int maxiter;
double re;//relaxation factor
double dt;//* 1e15;//dt >1e10 for steady state

public:
Shell();
Shell(double real, double img);
Shell(const Shell &C);
const Shell & operator=(const Shell &other);
printDetail();
double getReal()const {return M;}
double getImg()const {return N;}
};
ostream &operator<<(ostream& out,const Shell& C);

#endif
