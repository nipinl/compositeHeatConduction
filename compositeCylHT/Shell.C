#include "Shell.H"

Shell::Shell():
				axi(true),//axisymmetric, i.e., for tube
				M(10),//divisions along y(r)
				N(10),//divisions along x, i.e., along length
				r
				{}
Shell::Shell(double real,double img):x(real),y(img){
}
 Shell::Shell(const Shell &C){
    cout<<"invoking copy constructor !!"<<endl;
    x=C.x;y=C.y;
} 
const Shell & Shell::operator=(const Shell &other){
	this->x = other.x;
	this->y = other.y;
	return *this;
}

ostream &operator<<(ostream& out,const Shell& C) {
	out<<C.getReal()<<"i +"<<C.getImg()<< " j ";
	return out;
}
