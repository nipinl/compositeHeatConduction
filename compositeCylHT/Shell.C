#include "Shell.H"

Shell::Shell():
				axi(true),//axisymmetric, i.e., for tube
				M(10),//divisions along y(r)
				N(10),//divisions along x, i.e., along length
				ri(0.1),//inner radius
				Length(0.04),
				Width(0.04),
				tCond(2),spHeat(10),density(8000),
				simTime(20),
				initTemp(300),
				lbc(constTemp),rbc(constTemp),bbc(constTemp),tbc(constTemp),
				TLeft(300),TRight(300),TBottom(300),TTop(500),
				qLeft(0),qRight(0),qBottom(100),qTop(-50),
				TfLeft(300),TfRight(300),TfBottom(300),TfTop(500),
				hfL(10),hfR(10),hfB(10),hfT(50),
				maxiter(1000),
				re(1),
				dt(0.1)

				{}

//Shell::Shell(double real,double img):x(real),y(img){}
/*  Shell::Shell(const Shell &C){
    cout<<"invoking copy constructor !!"<<endl;
   // x=C.x;y=C.y;
}  */
//const Shell & Shell::operator=(const Shell &other){	this->x = other.x;this->y = other.y;return *this;}

ostream &operator<<(ostream& out,const Shell& C) {
	out<<C.getReal()<<"i +"<<C.getImg()<< " j ";
	return out;
}


void Shell::printDetail(){
				if (axi) cout<<"Cylindrical shell"<<endl;
				if (!axi) cout<<"2D rectangular shell"<<endl;

				if (axi) cout<<"Inner radius(m) :"<<ri<<endl;
				cout<<"Shell Length(m) :"<<Length<<endl;
				if (!axi) cout<<"Shell Width(m) :"<<Width<<endl;
				if (axi) cout<<"Shell Thickness(m) :"<<Width<<endl;
				cout<<"Material properties"<<endl;
				cout<<"-----------------------------------"<<endl;
				cout<<"Thermal Conductivity :"<<tCond<<endl;
				cout<<"Specific Heat :"<<spHeat<<endl;
				cout<<"Density :"<<density<<endl<<endl;

				cout<<"Simulation time(s) :"<<simTime<<endl;
				cout<<"Number of divisions along y :"<<M<<endl;
				cout<<"Number of divisions along x :"<<N<<endl;
				cout<<"Initial temperature(K) :"<<initTemp<<endl;
				cout<<"Boundary conditions :"<<endl;
				switch (lbc)
				{
				case constTemp:
					cout<<"Left boundary condition is constTemp of "<<TLeft<<" K"<<endl;
					break;
				case constHeatFlux:
					cout<<"Left boundary condition is constHeatFlux of "<<qLeft<<" watts"<<endl;
					break;
				case convection:
					cout<<"Left boundary condition is connvection with Tf = "<<TfLeft<<" K and hf = "<<hfL<<"W/m-k "<<endl;
					break;
				default:
					cout<<"Invalid boundary condition "<<endl;
					break;
				}
				switch (rbc)
				{
				case constTemp:
					cout<<"Right boundary condition is constTemp of "<<TRight<<" K"<<endl;
					break;
				case constHeatFlux:
					cout<<"Right boundary condition is constHeatFlux of "<<qRight<<" watts"<<endl;
					break;
				case convection:
					cout<<"Right boundary condition is connvection with Tf = "<<TfRight<<" K and hf = "<<hfR<<"W/m-k "<<endl;
					break;
				default:
					cout<<"Invalid boundary condition "<<endl;
					break;
				}
				switch (bbc)
				{
				case constTemp:
					cout<<"Bottom boundary condition is constTemp of "<<TBottom<<" K"<<endl;
					break;
				case constHeatFlux:
					cout<<"Bottom boundary condition is constHeatFlux of "<<qBottom<<" watts"<<endl;
					break;
				case convection:
					cout<<"Bottom boundary condition is connvection with Tf = "<<TfBottom<<" K and hf = "<<hfL<<"W/m-k "<<endl;
					break;
				default:
					cout<<"Invalid boundary condition "<<endl;
					break;
				}
				switch (tbc)
				{
				case constTemp:
					cout<<"Top boundary condition is constTemp of "<<TTop<<" K"<<endl;
					break;
				case constHeatFlux:
					cout<<"Top boundary condition is constHeatFlux of "<<qTop<<" watts"<<endl;
					break;
				case convection:
					cout<<"Top boundary condition is connvection with Tf = "<<TfTop<<" K and hf = "<<hfT<<"W/m-k "<<endl;
					break;
				default:
					cout<<"Invalid boundary condition "<<endl;
					break;
				}

				cout<<"Maximum iteration :"<<maxiter<<endl;
				cout<<"Relaxation coefficient :"<<re<<endl;
				cout<<"Time step(s) :"<<dt<<endl;
}
