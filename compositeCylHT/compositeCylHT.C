#include<iostream>
#include "Shell.H"
using namespace std;
int main()
{                                                                          
	Shell s1,s2,s3,s4;
	/* s1.setConstantTempBC("Top",500);
	s1.setConstantTempBC("Left",500);
	s2.setConstantTempBC("Top",500);
	s2.setConstantTempBC("Right",500);
	s3.setConstantTempBC("Bottom",500);
	s3.setConstantTempBC("Left",500);
	s4.setConstantTempBC("Bottom",500);
	s4.setConstantTempBC("Right",500); */
	
	vector<vector<Shell>> s
    {
        {s1,s2},
		{s3,s4}
    };
	
	solveSystem(s,"s",5);

	return 0;
}