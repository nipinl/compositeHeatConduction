#include<iostream>
#include "Shell.H"
using namespace std;
int main()
{                                                                          
	Shell s1,s2,s3,s4;
	s1.setSimulationTime(50);
	s1.setConstantTempBC("Top",1500);
	s1.setConstantTempBC("Left",1500);
	s1.setConstantTempBC("Bottom",1500);
	s1.setConstantTempBC("Right",1500);
	s1.setInitialTemp(300);

	s1.setThermalConductivity(16,0.1);
	//s3.setConstantTempBC("Bottom",500);
	//s3.setConstantTempBC("Left",500);
	//s4.setConstantTempBC("Bottom",500);
	//s4.setConstantTempBC("Right",500);

	/* s1.setBottomInterShellRadiationWithConvection(0.8,0.8,0.1,100,1000);
	s2.setBottomInterShellRadiation(0.8,0.8,0.1);
	//s2.setBottomInterShellContactResistance(1000);
	s3.setRightInterShellContactResistance(100);

	vector<vector<Shell>> s
    {
        {s1,s2},
		{s3,s4}
    }; */
	solveTransient(s1);
	//solveSystem(s);
	//solveTransient(s1);
	//s1.solveTransient();

	return 0;
}