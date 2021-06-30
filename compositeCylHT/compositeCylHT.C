#include<iostream>
#include "Shell.H"
using namespace std;
int main(){                                                                          
	Shell c1,c2,c3,c4;
	Shell Connection [2][2] = { {c1,c2} ,
								{c3,c4} };
	c1.setGeometry(0.5,0.05,true,0.15);
	c1.printDetail();
	c1.solveTransient();
	c2.solveSteady();
	return 0;
}
