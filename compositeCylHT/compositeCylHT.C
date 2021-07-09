#include<iostream>
#include "Shell.H"
using namespace std;
int main()
{                                                                          
	Shell s1,s2,s3,s4;
	s1.setSimulationTime(100);
	s1.setConstantTempBC("Top",1500);
	s1.setConstantTempBC("Left",1500);
	s2.setConstantTempBC("Top",1500);
	s2.setConstantTempBC("Right",1500);
	//s3.setConstantTempBC("Bottom",500);
	//s3.setConstantTempBC("Left",500);
	//s4.setConstantTempBC("Bottom",500);
	//s4.setConstantTempBC("Right",500);

	s1.setBottomInterShellRadiationWithConvection(0.8,0.8,0.1,100,1000);
	s2.setBottomInterShellRadiation(0.8,0.8,0.1);
	//s2.setBottomInterShellContactResistance(1000);
	s3.setRightInterShellContactResistance(100);
	
	vector<vector<Shell>> s
    {
        {s1,s2},
		{s3,s4}
    };
	
	solveSystem(s);

	return 0;
}