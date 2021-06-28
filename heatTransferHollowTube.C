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

//const temp bc
const double TLeft{300},TRight{300},TBottom{300},TTop{300};

//const heat flux bc
const double qLeft{0},qRight{0},qBottom{0},qTop{0};





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
			
			//tube inlet, which is along r(or y) direction            
			for (int j=0;j<M+2;j++){
				/* Const temperature*/ 
				te[j][0]=TLeft;
            	
				/* Const heat flux*/   
				te[j][0]=te[j][1]+0.5*dx[1]*qLeft/tk[j][1];//heat flux into the system is +ve
				
				/* Convection */ 
				te[j][0]=(hfL*TfLeft+2*te[j][1]*tk[j][1]/dx[1])/(hfL+2*tk[j][1]/dx[1]);
			}                           

			//tube outlet, which is along r(or y) direction
			for (int j=0;j<M+2;j++){
            	/* Const temperature*/ 
				te[j][N+1]=TRight;
            	/* Const heat flux*/
				te[j][N+1]=te[j][N]+0.5*dx[N]*qRight/tk[j][N];
				/* Convection */
				te[j][N+1]=(hfR*TfRight+2*te[j][N]*tk[j][N]/dx[N])/(hfR+2*tk[j][N]/dx[N]);
			}

			//inside tube, which is along z(or x) direction
			for (int i=0;i<N+2;i++){
            	/* Const temperature*/ 
				te[0][i]=TBottom;
            	/* Const heat flux*/ 
				te[0][i]=te[1][i]+0.5*dy[1]*qBottom/tk[1][i];
				/* Convection */ 
				te[0][i]=(hfB*TfBottom+2*te[1][i]*tk[1][i]/dy[1])/(hfB+2*tk[1][i]/dy[1]);
			}

			//outside tube, which is along z(or x) direction
			for (int i=0;i<N+2;i++){
            	/* Const temperature*/ 
				te[M+1][i]=TTop;
            	/* Const heat flux*/ 
				te[M+1][i]=te[M][i]+0.5*dy[M]*qTop/tk[M][i];
				/* Convection */
				te[M+1][i]=(hfT*TfTop+2*te[M][i]*tk[M][i]/dy[M])/(hfT+2*tk[M][i]/dy[M]);
			}

			//add source terms if any. Now nothing
			for (j=0;j<M+2;j++) {
                for (i=0;i<N+2;i++) {
                    sp[j][i] = 0;
                    sc[j][i] = 0;
                }
            }
// solve -------------------------------------------
			//marching in y
            for (j=1;j<M+1;j++) {      
				//marching in x
                for (i=1;i<N+1;i++) {  
                    if (geo==1){sae = dy[j];saw = sae;}
                    if (geo==2){sae = dy[j]*(y[j]+y[j+1])/2.0;saw = sae;}
                    if (geo==3){sae = dy[j]*x[i+1];saw = dy[j]*x[i];}
                    ke=tk[j][i]*tk[j][i+1]*(dx[i]+dx[i+1])/(dx[i]*tk[j][i+1] + dx[i+1]*tk[j][i]);
                    de = 2.0*ke*sae/(dx[i]+dx[i+1]);
                    ae=de;
                    kw=tk[j][i]*tk[j][i-1]*(dx[i]+dx[i-1])/(dx[i]*tk[j][i-1] + dx[i-1]*tk[j][i]);
                    dw = 2.0*kw*saw/(dx[i]+dx[i-1]);
                    aw=dw;
                    if (geo==1){san = dx[i];sas = san;}
                    if (geo==2){san = dx[i]*y[j+1];sas = dx[i]*y[j];}
                    if (geo==3){san = dx[i];sas = san;}

                    kn = tk[j][i] *tk[j+1][i]*(dy[j]+dy[j+1])/(dy[j]*tk[j+1][i]+dy[j+1]*tk[j][i]);
                    dn = 2.0*kn*san/(dy[j]+dy[j+1]);
                    an=dn;
                    ks = tk[j][i] *tk[j-1][i]*(dy[j]+dy[j-1])/(dy[j]*tk[j-1][i]+dy[j-1]*tk[j][i]);
                    ds = 2.0*ks*sas/(dy[j]+dy[j-1]);
                    as=ds;
                    if (geo==1){vol=dx[i]*dy[j];}
                    if (geo==2){vol=dx[i]*dy[j]*(y[j]+y[j+1])/2.0;}
                    if (geo==3){vol=dx[i]*dy[j]*(x[i]+x[i+1])/2.0;}
                    a0 = rho[j][i]*cp[j][i]*vol/dt;// =0 for steady state
                    ap = ae+aw+an+as+a0 - sp[j][i]*vol;
                    b = sc[j][i]*vol+a0*te0[j][i];
                    ta[0][i]=ap/re;
                    tb[0][i]=ae;
                    tc[0][i]=aw;
                    td[0][i]=b+ap/re*(1-re)*te[j][i]+an*te[j+1][i]+as*te[j-1][i];
                    if(i==1) td[0][i]=td[0][i]+aw*te[j][0];
                    if(i==n) td[0][i]=td[0][i]+ae*te[j][N+1];
                }//marching in x
                //start of tdma
                beta[1]=tb[0][1]/ta[0][1];
                alpha[1]=td[0][1]/ta[0][1];
                //forward substitution
                for (ii=2;ii<N+1;ii++){
                    beta[ii]=tb[0][ii]/(ta[0][ii] - tc[0][ii]*beta[ii-1]);
                    alpha[ii]=(td[0][ii]+tc[0][ii]*alpha[ii-1])/(ta[0][ii] - tc[0][ii]*beta[ii-1]);
                }
                //backward substitution
                dum[N]=alpha[N];
                for (jj=0;jj<N-1;jj++){
                    ii=N-1-jj;
                    dum[ii]=beta[ii]*dum[ii+1]+alpha[ii];
                }
                //end of tdma
                for (i=1;i<N+1;i++){
                    te[j][i] = dum[i];
                }
            }//marching in y
            //solve ********************************************
            //start convergence checking ---------------------
            maxErr=1e-8;
            for (j=0;j<M+2;j++) {
                for (i=0;i<N+2;i++) {
                    errorTe[j][i] = Math.abs(te[j][i]-tep[j][i])/te[j][i];
                    if (errorTe[j][i]>maxErr) maxErr =errorTe[j][i];
                }
            }
            if(maxErr>error){
                iter++;
                for (i=0;i<N+2;i++) {
                    for (j=0;j<M+2;j++) {
                        tep[j][i]=te[j][i];
                    }
                }
                iflag=1;
            }
            if(maxErr<=error){
                iflag=0;
                cout<<" Converged. MaxErr = "<< maxErr<<" iter = "<<iter<<" Time = "<<t<<endl;
            }
            //end convergence checking ***********************

            if (iter>maxiter) {
                cout<<" Iterations need to be inreased. Error in Temp is "<<maxErr*100<<" %"<<endl;
                break;
            }

        }//end of inner while loop checking iflag
		//......................Inner Loop ........................................


        t=t+dt;//increment time step
        for (i=0;i<N+2;i++) {
            for (j=0;j<M+2;j++) {
                te0[j][i]=te[j][i];
                tep[j][i]=te0[j][i];
            }
        }
        if (iwrite>mwrite) {
            //Toast.makeText(getApplicationContext(),"Iteration no:" +iter+" and time : "+t,Toast.LENGTH_SHORT).show();
            iwrite=0;
        }
        iwrite++;
	}//end of while loop checking t<simTime
	//--------------------Outer Loop ---------------------------------------------

	
	return 0;
}
