#include<iostream>
#include "Shell.H"
using namespace std;
//A single shell case with variable thermal conductivity
int main()
{                                                                          
	Shell s1,s2,s3,s4;
	s1.cylindrical(false);//if not set, shell will be cylindrical type
	s1.setSimulationTime(10);
	s1.setConstantTempBC("Top",1500);
	s1.setConstantTempBC("Left",1500);
	s1.setConstantTempBC("Bottom",1500);
	s1.setConstantTempBC("Right",1500);
	s1.setInitialTemp(300);

	s1.setThermalConductivity(16,0.1);
	//the following will override the above definition of thermal conductivity
	
	vector <double> T ={100,3000,400,300,200};// need not be in ascending or descending order
	vector <double> TCx ={2, 20, 8, 6, 4}; 
	vector <double> TCy ={1, 10, 4, 3, 2}; 
	s1.setVariableThermalConductivity(0,T,TCx);//0: direction along length (x)  1: direction along width (y or r)
	s1.setVariableThermalConductivity(1,T,TCy);

	solveTransient(s1);
	s1.printDetail();

	return 0;
}

//--------------------------------------------------------------------------------------
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
//-----------------------------------END---------------------------------------------------

//comments
// 1. InterShellBoundaryCondition  Contact Resistance + Radiation/ Radiation with Conv can be thought of, but not implemented.