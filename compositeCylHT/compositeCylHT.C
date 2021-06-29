#include<iostream>
#include "Shell.H"
using namespace std;
int main(){
	double simTime = 20;
	double dt = 0.1;//* 1e15;//dt >1e10 for steady state
	Shell c1;
	c1.setTimes(simTime,dt);
	c1.printDetail();
	c1.populateNodes();
	c1.populateMaterialProperties();
	c1.initialiseField();
	double t=0;//time
	while(t<simTime){
		c1.solveIt();
		t=t+dt;//increment time step
	}
	return 0;
}
