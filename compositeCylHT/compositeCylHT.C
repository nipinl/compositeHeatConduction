#include<iostream>
#include "Shell.H"
using namespace std;
int main(){
	double simTime = 20;
	double dt = 0.2;//* 1e15;//dt >1e10 for steady state
	double initialTemp = 300;
	Shell c1;
	c1.setTimes(simTime);
	c1.setMaterialProperties(10);
	c1.setConvectionBC("Top",500,100);
	c1.setConvectionBC("Bottom",400,90);
	c1.setConvectionBC("Left",300,80);
	c1.setConvectionBC("Right",200,70);
	c1.printDetail();
	 c1.populateNodes();
	c1.populateMaterialProperties();
	c1.initialiseField(800);
	double t=0;//time
	while(t<simTime){
		c1.solveIt();
		t=t+dt;//increment time step
	} 
	return 0;
}
