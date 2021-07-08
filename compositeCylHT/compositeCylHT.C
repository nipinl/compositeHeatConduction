#include<iostream>
#include "Shell.H"
using namespace std;
int main()
{                                                                          
	Shell s1,s2,s3,s4,s5,s6,s7,s8,s9;

	s5.setInitialTemp(1000);
	s3.setMaterialProperties(0.1);
	s4.setMaterialProperties(0.1);
	s6.setMaterialProperties(0.1);
	s7.setMaterialProperties(0.1);

	s7.setConstantTempBC("Bottom",500);
	s8.setConstantTempBC("Bottom",500);
	s9.setConstantTempBC("Bottom",500);
	vector<vector<Shell>> s
    {
        {s1,s2,s3},
		{s4,s5,s6},
		{s7,s8,s9}
    };
	
	solveSystem(s);
	/* s2.setSimulationTime(0.08);
	s2.setConstantTempBC("Top",1000);
	s2.setConstantTempBC("Right",1000);
	s2.solveTransient();
	s2.printDetail(); */

	return 0;
}