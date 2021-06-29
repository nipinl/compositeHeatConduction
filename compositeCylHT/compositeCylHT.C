#include<iostream>
#include "Shell.H"
using namespace std;

int main(){
	Shell c1;
	Shell c2(2.0,3.0);
	Shell c3(c1);
	Shell c4=c2;
	c1 = c4;
	cout<<c4<<endl<<c3<<endl;

	return 0;
}
