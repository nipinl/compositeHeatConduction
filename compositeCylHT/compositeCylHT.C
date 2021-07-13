#include<iostream>
#include "Shell.H"
using namespace std;
int main()
{                                                                          
	Shell s1,s2,s3,s4;
	s1.setSimulationTime(10.15);
	s1.setConstantTempBC("Top",1500);
	s1.setConstantTempBC("Left",1500);
	s1.setConstantTempBC("Bottom",1500);
	s1.setConstantTempBC("Right",1500);
	s1.setInitialTemp(300);

	s1.setThermalConductivity(16,0.1);
	vector <double> temp ={100,400,300,200};
	vector <double> tkx ={10,40,30,20}; 
	vector <double> tky ={100,400,300,200}; 
	s1.setVariableThermalConductivity(0,temp,tkx);
	s1.setVariableThermalConductivity(1,temp,tky);

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
	s1.printDetail();
	//solveSystem(s);
	//solveTransient(s1);
	//s1.solveTransient();

	return 0;
}