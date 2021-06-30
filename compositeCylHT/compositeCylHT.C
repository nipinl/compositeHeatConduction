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
	connectShells(c1,c2);

	c1.populateNodes();
	c1.populateMaterialProperties();
	c1.initialiseField();
	c1.printDetail();

	c2.setConstantTempBC("Left",100);
	c2.setConstantTempBC("Right",300);
	c2.setConstantTempBC("Bottom",400);
	c2.setConstantTempBC("Top",200);
	c2.populateNodes();
	c2.populateMaterialProperties();
	c2.initialiseField();
	double t=0;//time
	while(t<simTime){
		c1.solveIt();
		c2.solveIt();
		t=t+dt;//increment time step
	} 
	return 0;
}
