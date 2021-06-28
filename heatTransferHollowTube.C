#include<iostream>
#include<vector>
using namespace std;
bool debug = true;
const int M = 100;//Number of divisions in y direction
const int N = 100;//Number of divisions in x direction
const double r1 = 0.03;//Inner radius of the tube
const double t1 = 0.01;//thickness of the tube
const double Length = 0.1;//Length of the cylinder

//Material properties
const double tCond = 1000.0; 
const double spHeat = 1000;
const double density = 1000.0;

const double simTime = 100.0;//in seconds
const double initTemp = 300.0;//Initial uniform temperature of the tube

const double TLeft{300},TRight{300},TBottom{300},TTop{300};




int main(){
	//Preprocessing the data

	//Define arrays for x, y, dx and dy
	double x [N+3]{0}, y[M+3]{0}, T[M+2][N+2]{0};
	double dx[N+2]{Length/N}, dy[M+2]{t1/M};

	//gost cells, hence giving almost zero thickness
	dx[0] = 1e-10;dx[N+1] = 1e-10;dy[0] = 1e-10;dy[M+1] = 1e-10;
	if (debug) cout<<"dx = "<<dx[0]<<" and dy = "<<dy[0]<<endl;//Just a check point
	
	//populate x and y
	for (int i=0;i<N+2;i++){x[i+1] = x[i]+dx[i];}
	for (int i=0;i<M+2;i++){y[i+1] = y[i]+dy[i];}
	if (debug) cout<<"x[0] = "<<x[0]<<", x[1] = "<<x[1]<<" and y[0] = "<<y[0]<<", y[1] = "<<y[1]<<endl;//Just another check point

	// keep the option of variable thermal conductivity, Cp and density
	double tk[M+2][N+2]{tCond},  cp[M+2][N+2]{spHeat},  rho[M+2][N+2]{density};

	// temperature variables for storing at different points of the algorithm
	double te0[M+2][N+2]{initTemp},  te[M+2][N+2]{initTemp},  tep[M+2][N+2]{initTemp};
	
	//source term modelling
	double sc[M+2][N+2]{0},sp[M+2][N+2]{0};

	//for TDMA
	double ta[M+1][N+1]{0},  tb[M+1][N+1]{0},  tc[M+1][N+1]{0},  td[M+1][N+1]{0};
	double alpha[N+2]{0}, beta[N+2]{0},  dum[N+2]{0};

	//for convergence checking
	double errorTe[M+2][N+2]{0};

	//defining variables for equation
	double ke{0},  kw{0},  ks{0},  kn{0};
	double de{0},  dw{0},  ds{0},  dn{0};
	double ae{0},  aw{0},  as{0},  an{0};
	double a0{0},  ap{0},  b{0},  vol{0};

	//variables required for loop
	int iter{0}, iflag{1};
	double t=0;//time

	//--------------------Outer Loop ---------------------------------------------
	while(t<simTime){
		iter = 0;
		iflag =1;

		//......................Inner Loop ........................................
		while(iflag==1){
			//update properties if are dependent on temperature. Currently it is not

			//Enforce boundary conditions
			
			//tube top, which is along r(or y) direction
            /* Const temperature*/ for (int j=0;j<M+2;j++){te[j][0]=TLeft;}
            /* Const heat flux*/   for (int j=0;j<M+2;j++){te[j][0]=te[j][1]+0.5*dx[1]*qLeft/tk[j][1];}
			                            


			//tube bottom, which is along r(or y) direction
            /* Const temperature*/ for (int j=0;j<M+2;j++){te[j][N+1]=TRight;}
            /* Const heat flux*/   for (int j=0;j<M+2;j++){te[j][N+1]=te[j][N]+0.5*dx[N]*qRight/tk[j][N];}


			//inside tube , which is along z(or x) direction
            /* Const temperature*/ for (int i=0;i<N+2;i++){te[0][i]=TBottom;}

			//outside tube, which is along z(or x) direction
            /* Const temperature*/ for (int i=0;i<N+2;i++){te[M+1][i]=TTop;}

	
	






		}
		//......................Inner Loop ........................................


	}
	//--------------------Outer Loop ---------------------------------------------

	
	return 0;
}
