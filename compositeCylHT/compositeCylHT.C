#include<iostream>
#include "Shell.H"
using namespace std;
int main()
{                                                                          
	Shell s1,s2,s3,s4,s5,s6,s7,s8,s9, s10, s11, s12;
    
	vector<vector<Shell>> Shells
    {
        {s1, s2, s3, s4},
        {s5, s6, s7, s8},
		{s9, s10, s11, s12}
    };


	//s1.setGeometry(0.5,0.05,true,0.15);
	s1.setRightInterShellContactResistance(100);
	s1.setBottomInterShellRadiation();
	s2.setBottomInterShellRadiationWithConvection();
	//s1.solveTransient();
	//s2.solveSteady();
	s2.setInitialTemp(5000);
	s2.setConvectionBC();
	s3.setInitialTemp(5000);
	s3.setConvectionRadiationBC();
	vector<vector<Shell>> s
    {
        {s1,s2},
		{s3,s4}
    };
	solveSystem(s);
	//s2.solveTransient();
	//s3.solveTransient();
	//s3.printTe();
    	
	return 0;
}