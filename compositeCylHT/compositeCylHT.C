#include<iostream>
#include "Shell.H"
using namespace std;
//A single shell case with variable thermal conductivity
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
	//the following will override the above definition of thermal conductivity
	
	vector <double> temp ={100,3000,400,300,200};
	vector <double> tkx ={2, 20, 8, 6, 4}; 
	vector <double> tky ={1, 10, 4, 3, 2}; 
	s1.setVariableThermalConductivity(0,temp,tkx);
	s1.setVariableThermalConductivity(1,temp,tky);

	solveTransient(s1);
	s1.printDetail();

	return 0;
}

//A system solution case
/* int main()
{                                                                          
	Shell s1,s2,s3,s4;
	
	s3.setConstantTempBC("Bottom",500);
	s3.setConstantTempBC("Left",500);
	s4.setConstantTempBC("Bottom",500);
	s4.setConstantTempBC("Right",500);

	s1.setBottomInterShellRadiationWithConvection(0.8,0.8,0.1,100,1000);
	s2.setBottomInterShellRadiation(0.8,0.8,0.1);

	vector<vector<Shell>> s
    {
        {s1,s2},
		{s3,s4}
    };

	solveSystem(s);
	return 0;
} */