#include "Shell.H"

//Constructor
Shell::Shell():
				axi(true),//axisymmetric, i.e., for tube
				M(5),//divisions along y(r)
				N(3),//divisions along x, i.e., along length
				ri(0.1),//inner radius
				Length(0.04),
				Width(0.04),
				tCond(16),spHeat(500),density(8000),
				initTemp(300),
				lbc(constTemp),rbc(constTemp),bbc(constTemp),tbc(constTemp),
				TLeft(300),TRight(300),TBottom(300),TTop(500),
				qLeft(0),qRight(0),qBottom(100),qTop(-50),
				TfLeft(300),TfRight(300),TfBottom(300),TfTop(500),
				hfL(10),hfR(10),hfB(10),hfT(50),
				maxiter(1000),
				re(1),
				simTime(20),
				dt(0.1)
				{}

void Shell::setTimes(double simtime, double delt){
	simTime = simtime;
	dt = delt;
}
void Shell::setMaterialProperties(double Thermal_cond,double Cp,double Density){
	tCond = Thermal_cond;
	spHeat = Cp;
	density = Density;
}
void Shell::setConstantTempBC(string boundary, double Temp){
	if (boundary=="Left"){lbc = constTemp; TLeft = Temp;}
	else if (boundary=="Right"){rbc = constTemp; TRight = Temp;}
	else if (boundary=="Bottom"){bbc = constTemp; TBottom = Temp;}
	else if (boundary=="Top"){tbc = constTemp; TTop = Temp;}
	else{ cout<<"Wrong boundary specified in constantTemp boundary condition"<<endl;exit(1);}
	
}
void Shell::setConstantHeatfluxBC(string boundary, double q){
	if (boundary=="Left"){lbc = constHeatFlux; qLeft = q;}
	else if (boundary=="Right"){rbc = constHeatFlux; qRight = q;}
	else if (boundary=="Bottom"){bbc = constHeatFlux; qBottom = q;}
	else if (boundary=="Top"){tbc = constHeatFlux; qTop = q;}
	else{ cout<<"Wrong boundary specified in constantHeatflux boundary condition"<<endl;exit(1);}
}
void Shell::setConvectionBC(string boundary, double Tinf, double h){
	if (boundary=="Left"){lbc = convection; TfLeft = Tinf, hfL = h;}
	else if (boundary=="Right"){rbc = convection; TfRight = Tinf, hfR = h;}
	else if (boundary=="Bottom"){bbc = convection; TfBottom = Tinf, hfB = h;}
	else if (boundary=="Top"){tbc = convection; TfTop = Tinf, hfT = h;}
	else{ cout<<"Wrong boundary specified in Convection boundary condition"<<endl;exit(1);}
	
}
void Shell::populateNodes(){
	//populate dx, dy
	for (int i=0;i<N+2;i++){dx.push_back(Length/double(N));}
	for (int j=0;j<M+2;j++){dy.push_back(Width/double(M));}
	dx[0] = 1e-10;dx[N+1] = 1e-10;dy[0] = 1e-10;dy[M+1] = 1e-10;
	
	//populate x & y
	for (int i=0;i<N+3;i++){
		if (i==0) x.push_back(0);
		if (i>0) x.push_back(x[i-1]+dx[i-1]);
	}
	for (int j=0;j<M+3;j++){
		if (j==0){
			if (!axi) y.push_back(0);
			if (axi) y.push_back(ri);
		}
		if (j>0) y.push_back(y[j-1]+dy[j-1]);
	}
	/* cout<<  "dx:" <<endl<<endl;
	for (int i=0;i<N+2;i++){
		cout<<dx[i]<<endl;
	}
	
	cout<<  "dy:" <<endl<<endl;
	for (int j=0;j<M+2;j++){
		cout<<dy[j]<<endl;
	}
	
	cout<<  "x:" <<endl<<endl;
	for (int i=0;i<N+3;i++){
		cout<<x[i]<<endl;
	}

	cout<<  "y:" <<endl<<endl;
	for (int j=0;j<M+3;j++){
		cout<<y[j]<<endl;
	} */
}

void Shell::populateMaterialProperties(){
	//std::vector<std::vector<double>> tk(M, std::vector<double>(N, tCond));
	vector<double> tk1d(N+2, tCond);
	vector<double> cp1d(N+2, spHeat);
	vector<double> rho1d(N+2, density);
	for (int i = 0 ; i < M+2 ; i++) {
        tk.push_back(tk1d);
		cp.push_back(cp1d);
		rho.push_back(rho1d);		
    }
	Shell::print2dVector(tk);
	Shell::print2dVector(cp);
	Shell::print2dVector(rho);
}
void Shell::initialiseField(){
	vector<double> te1d(N+2, initTemp);
	vector<double> scsp1d(N+2, 0);
	for (int i = 0; i < N+1; i++)
	{
		ta.push_back(0);
		tb.push_back(0);
		tc.push_back(0);
		td.push_back(0);
	}

	for (int i = 0 ; i < M+2 ; i++) {
        te.push_back(te1d);
		te0.push_back(te1d);
		tep.push_back(te1d);
		sc.push_back(scsp1d);	
		sp.push_back(scsp1d);
    }
	Shell::print2dVector(te);
	Shell::print2dVector(te0);
	Shell::print2dVector(tep);
}
void::Shell::solveIt(){
	//for convergence checking
	double maxErr = 1e-10;
    double error = 1.0e-9;

	//defining variables for equation
	double ke{0},  kw{0},  ks{0},  kn{0};
	double de{0},  dw{0},  ds{0},  dn{0};
	double ae{0},  aw{0},  as{0},  an{0};
	double a0{0},  ap{0},  b{0},  vol{0};
	double sae{0},  saw{0}, san{0}, sas{0};
	//variables required for loop
	int iter{0}, iflag = 1;
	

	//--------------------Outer Loop ---------------------------------------------
	
		
		//......................Inner Loop ........................................
		while(iflag==1){
			//update properties if are dependent on temperature. Currently it is not

			//Enforce boundary conditions
			
			 


			//add source terms if any. Now nothing
			for (int j=0;j<M+2;j++) {
                for (int i=0;i<N+2;i++) {
                    sp[j][i] = 0;
                    sc[j][i] = 0;
                }
            }
			applyBoundaryConditions();
// solve -------------------------------------------
			//START marching in y
            for (int j=1;j<M+1;j++) {      
				//START marching in x
                for (int i=1;i<N+1;i++) {  
                    if (!axi){sae = dy[j];saw = sae;}
                    if (axi){sae = dy[j]*(y[j]+y[j+1])/2.0;saw = sae;}
                    ke=tk[j][i]*tk[j][i+1]*(dx[i]+dx[i+1])/(dx[i]*tk[j][i+1] + dx[i+1]*tk[j][i]);
                    de = 2.0*ke*sae/(dx[i]+dx[i+1]);
                    ae=de;
                    kw=tk[j][i]*tk[j][i-1]*(dx[i]+dx[i-1])/(dx[i]*tk[j][i-1] + dx[i-1]*tk[j][i]);
                    dw = 2.0*kw*saw/(dx[i]+dx[i-1]);
                    aw=dw;
                    if (!axi){san = dx[i];sas = san;}
                    if (axi){san = dx[i]*y[j+1];sas = dx[i]*y[j];}

                    kn = tk[j][i] *tk[j+1][i]*(dy[j]+dy[j+1])/(dy[j]*tk[j+1][i]+dy[j+1]*tk[j][i]);
                    dn = 2.0*kn*san/(dy[j]+dy[j+1]);
                    an=dn;
                    ks = tk[j][i] *tk[j-1][i]*(dy[j]+dy[j-1])/(dy[j]*tk[j-1][i]+dy[j-1]*tk[j][i]);
                    ds = 2.0*ks*sas/(dy[j]+dy[j-1]);
                    as=ds;
                    if (!axi){vol=dx[i]*dy[j];}
                    if (axi){vol=dx[i]*dy[j]*(y[j]+y[j+1])/2.0;}
                    a0 = rho[j][i]*cp[j][i]*vol/dt;// =0 for steady state
                    ap = ae+aw+an+as+a0 - sp[j][i]*vol;
                    b = sc[j][i]*vol+a0*te0[j][i];
                    ta[i]=ap/re;	
                    tb[i]=ae;
                    tc[i]=aw;
                    td[i]=b+ap/re*(1-re)*te[j][i]+an*te[j+1][i]+as*te[j-1][i];
                    if(i==1) td[i]=td[i]+aw*te[j][0];
                    if(i==N) td[i]=td[i]+ae*te[j][N+1];
                }//marching in x ends here
				//END marching in x
                //solve in x direction using tdma
				tdma(j);
                
            }//marching in y ends here
			//END marching in y
//solve ----------------------------------------------------------------------------------
            //start convergence checking ---------------------
            iflag = checkConvergence(error);
			if(iflag==1) iter++;
			if (iflag==0) cout<<"  Solution Converged for one time step. dt = "<<dt<<" in "<<iter<<" iterations"<<endl;
            //end convergence checking ***********************

            if (iter>maxiter) {
                cout<<" Iterations need to be inreased. Error in Temp is "<<maxErr*100<<" %"<<endl;
                break;
            }

        }//end of inner while loop checking iflag
		//......................Inner Loop ........................................


        
        for (int i=0;i<N+2;i++) {
            for (int j=0;j<M+2;j++) {
                te0[j][i]=te[j][i];
                tep[j][i]=te0[j][i];
            }
        }
        
		for (int i=0;i<N+2;i++) {
         	for (int j=0;j<M+2;j++) {
                cout<<te[j][i]<<"  ";
           	}
			cout<<endl;
       	}
            
        
	//end of while loop checking t<simTime
	//--------------------Outer Loop ---------------------------------------------

}
void Shell::applyBoundaryConditions(){
	//tube inlet, which is along r(or y) direction            
			for (int j=0;j<M+2;j++){
				switch (lbc)
				{
				case constTemp:
					te[j][0]=TLeft;
					break;
				case constHeatFlux:
					te[j][0]=te[j][1]+0.5*dx[1]*qLeft/tk[j][1];//heat flux into the system is +ve
					break;
				case convection:
					te[j][0]=(hfL*TfLeft+2*te[j][1]*tk[j][1]/dx[1])/(hfL+2*tk[j][1]/dx[1]);
					break;
				default:
					te[j][0]=TLeft;
					break;
				}
			}                           

			//tube outlet, which is along r(or y) direction
			for (int j=0;j<M+2;j++){
				switch (rbc)
				{
				case constTemp:
					te[j][N+1]=TRight;
					break;
				case constHeatFlux:
					te[j][N+1]=te[j][N]+0.5*dx[N]*qRight/tk[j][N];//heat flux into the system is +ve
					break;
				case convection:
					te[j][N+1]=(hfR*TfRight+2*te[j][N]*tk[j][N]/dx[N])/(hfR+2*tk[j][N]/dx[N]);
					break;
				default:
					te[j][0]=TRight;
					break;
				}
			} 


			//inside tube, which is along z(or x) direction
			for (int i=0;i<N+2;i++){
				switch (bbc)
				{
				case constTemp:
					te[0][i]=TBottom;
					break;
				case constHeatFlux:
					te[0][i]=te[1][i]+0.5*dy[1]*qBottom/tk[1][i];//heat flux into the system is +ve
					break;
				case convection:
					te[0][i]=(hfB*TfBottom+2*te[1][i]*tk[1][i]/dy[1])/(hfB+2*tk[1][i]/dy[1]);
					break;
				default:
					te[0][i]=TBottom;
					break;
				}
			} 
			

			//outside tube, which is along z(or x) direction
			for (int i=0;i<N+2;i++){
				switch (tbc)
				{
				case constTemp:
					te[M+1][i]=TTop;
					break;
				case constHeatFlux:
					te[M+1][i]=te[M][i]+0.5*dy[M]*qTop/tk[M][i];//heat flux into the system is +ve
					break;
				case convection:
					te[M+1][i]=(hfT*TfTop+2*te[M][i]*tk[M][i]/dy[M])/(hfT+2*tk[M][i]/dy[M]);
					break;
				default:
					te[M+1][i]=TTop;
					break;
				}
			}
}
void Shell::tdma(int j){
	double alpha[N+2]{0}, beta[N+2]{0},  dum[N+2]{0};
	for(int i = 0;i<N+2;i++){
		alpha[i] = 1;
		beta[i] = 1;
		dum[i] = 1;
	}
	beta[1]=tb[1]/ta[1];
    alpha[1]=td[1]/ta[1];
    //forward substitution
    for (int ii=2;ii<N+1;ii++){
        beta[ii]=tb[ii]/(ta[ii] - tc[ii]*beta[ii-1]);
        alpha[ii]=(td[ii]+tc[ii]*alpha[ii-1])/(ta[ii] - tc[ii]*beta[ii-1]);
    }
    //backward substitution
    dum[N]=alpha[N];
    for (int jj=0;jj<N-1;jj++){
        int ii=N-1-jj;
        dum[ii]=beta[ii]*dum[ii+1]+alpha[ii];
    }
    //solved value
    for (int i=1;i<N+1;i++){te[j][i] = dum[i];}
}
int Shell::checkConvergence(double error){
	double maxErr{1e-10},errorTe{0};
	int iflag =1;

    for (int j=0;j<M+2;j++) {
        for (int i=0;i<N+2;i++) {
            errorTe = abs(te[j][i]-tep[j][i])/te[j][i];
            if (errorTe>maxErr) maxErr = errorTe;
        }
    }
    if(maxErr>error){
        for (int i=0;i<N+2;i++) {
            for (int j=0;j<M+2;j++) {
                tep[j][i]=te[j][i];
            }
        }
        iflag=1;
    }
    if(maxErr<=error){
        iflag=0;
		cout<<"MaxErr = "<< maxErr;
        
    }
	return iflag;
}
void Shell::print2dVector(vector<vector<double>> const &v){
	 for (auto i : v) {
     //i is now an 1D vector
     	for (auto j : i) {
            cout << j << " ";
        }
        cout << endl;
    }
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
