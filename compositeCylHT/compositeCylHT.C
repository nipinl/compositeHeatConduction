#include<iostream>
#include "Shell.H"
using namespace std;

int main(){
	Shell c1;
	Shell c2;
	c1.printDetail();
	c1.populateNodes();
	c1.populateMaterialProperties();

	return 0;
}
