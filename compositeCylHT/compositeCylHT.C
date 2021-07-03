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
    int shellNo=0;
    for (int i = 0; i < Shells.size(); i++){
        for (int j = 0; j < Shells[i].size(); j++){
			shellNo++;
			Shells[i][j].setLength(shellNo);
            Shells[i][j].setWidth(shellNo);
        }  
    }

	//s1.setGeometry(0.5,0.05,true,0.15);
	s1.setLength(2*s1.getLength());//doubling length
	s1.printDetail();
	//s1.solveTransient();
	//s2.solveSteady();
	s2.setInitialTemp(5000);
	s2.setConvectionBC();
	s3.setInitialTemp(5000);
	s3.setConvectionRadiationBC();
	vector<vector<Shell>> s
    {
        {s2}
    };
	//solveSystem(s);
	s2.solveTransient();
	s3.solveTransient();
	cout<<"-----s2-------"<<endl;
	s2.printTe();
	cout<<"-----s3-------"<<endl;
	s3.printTe();
    	
	return 0;
}