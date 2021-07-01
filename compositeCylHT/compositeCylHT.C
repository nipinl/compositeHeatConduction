#include<iostream>
#include "Shell.H"
using namespace std;
int main(){                                                                          
	Shell s1,s2,s3,s4,s5,s6;
	vector<vector<Shell>> Shells
    {
        {s1, s2},
        {s3, s4},
		{s5, s6}
    };
/* 
	s1.setGeometry(0.5,0.05,true,0.15);
	s1.printDetail();
	s1.solveTransient();
	s2.solveSteady(); */
	solveSystem(Shells);
	return 0;
}