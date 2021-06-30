#include<iostream>
#include "Shell.H"
using namespace std;
int main(){
	double simTime = 20;
	double dt = 0.2;//* 1e15;//dt >1e10 for steady state
	double initialTemp = 300;                                                                          
	Shell c1,c2,c3,c4;
	Shell Connection [2][2] = { {c1,c2} ,
								{c3,c4} };
	c2.setGeometry(0.5,0.05,true,0.15);


	c1.printDetail();
	
	c1.preprocessShell();	
	c2.preprocessShell();

	double t=0;//time
	while(t<simTime){
		c1.advanceOneTimeStep();
		c2.advanceOneTimeStep();
		t=t+dt;//increment time step
	} 
	return 0;
}
