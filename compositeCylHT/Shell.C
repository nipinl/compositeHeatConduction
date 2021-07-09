#include "Shell.H"

//Constructor
Shell::Shell() : axi(true),		   //axisymmetric, i.e., for tube
				 connected(false), //this shell is not connected to other
				 M(10),			   //divisions along y(r)
				 N(8),			   //divisions along x, i.e., along length
				 ri(0.1),		   //inner radius
				 Length(0.04),
				 Width(0.02),
				 tCond(16), spHeat(500), density(8000),
				 initTemp(300),
				 lbc(constTemp), rbc(constTemp), bbc(constTemp), tbc(constTemp),
				 TLeft(300), TRight(300), TBottom(300), TTop(300),
				 qLeft(0), qRight(0), qBottom(100), qTop(-50),
				 //TLeft(300), TRight(300), TBottom(300), TTop(500),
				 hfL(10), hfR(10), hfB(10), hfT(50),
				 eps(0.8),
				 maxiter(1000),
				 re(1),
				 simTime(20),
				 dt(0.1)
{}

//setters
void Shell::setConnected() { connected = true; }
void Shell::setSimulationTime(double simtime)
{
	simTime = simtime;
}
void Shell::setTimeStep(double delt)
{
	dt = delt;
}
void Shell::setGeometry(double length, double width, bool cylindrical, double innerRadius)
{
	bool err = false;
	if (length <= 0)
	{
		cout << "Length shall be realistic" << endl;
		err = true;
	}
	if (width <= 0)
	{
		cout << "Width shall be realistic" << endl;
		err = true;
	}
	if (err)
		exit(1);
	Length = length;
	Width = width;
	axi = cylindrical;
	ri = innerRadius;
	if (!cylindrical)
		ri = -1;
}
void Shell::setLength(double length){
	if (length <= 0){
		cout << "Length shall be realistic" << endl;
		exit(1);
	}
	Length = length;
}
void Shell::setWidth(double width){
	if (width <= 0){
		cout << "Width shall be realistic" << endl;
		exit(1);
	}
	Width = width;
}
void Shell::setInnerRadius(double innerRadius){
	if (innerRadius <= 0){
		cout << "InnerRadius value shall be realistic" << endl;
		exit(1);
	}
	ri = innerRadius;
}
void Shell::setMaterialProperties(double Thermal_cond, double Cp, double Density)
{
	tCond = Thermal_cond;
	spHeat = Cp;
	density = Density;
}
void Shell::setConstantTempBC(string boundary, double Temp)
{
	if (boundary == "Left")
	{
		lbc = constTemp;
		TLeft = Temp;
	}
	else if (boundary == "Right")
	{
		rbc = constTemp;
		TRight = Temp;
	}
	else if (boundary == "Bottom")
	{
		bbc = constTemp;
		TBottom = Temp;
	}
	else if (boundary == "Top")
	{
		tbc = constTemp;
		TTop = Temp;
	}
	else
	{
		cout << "Wrong boundary specified in constantTemp boundary condition" << endl;
		exit(1);
	}
}
void Shell::setConstantHeatfluxBC(string boundary, double q)
{
	if (boundary == "Left")
	{
		lbc = constHeatFlux;
		qLeft = q;
	}
	else if (boundary == "Right")
	{
		rbc = constHeatFlux;
		qRight = q;
	}
	else if (boundary == "Bottom")
	{
		bbc = constHeatFlux;
		qBottom = q;
	}
	else if (boundary == "Top")
	{
		tbc = constHeatFlux;
		qTop = q;
	}
	else
	{
		cout << "Wrong boundary specified in constantHeatflux boundary condition" << endl;
		exit(1);
	}
}
void Shell::setConvectionBC(string boundary, double Tinf, double h)
{
	if (boundary == "Left")
	{
		lbc = convection;
		TLeft = Tinf, hfL = h;
	}
	else if (boundary == "Right")
	{
		rbc = convection;
		TRight = Tinf, hfR = h;
	}
	else if (boundary == "Bottom")
	{
		bbc = convection;
		TBottom = Tinf, hfB = h;
	}
	else if (boundary == "Top")
	{
		tbc = convection;
		TTop = Tinf, hfT = h;
	}
	else
	{
		cout << "Wrong boundary specified in Convection boundary condition" << endl;
		exit(1);
	}
}
void Shell::setConvectionRadiationBC(string boundary, double Tamb, double h, double emissivity)
{
	if (boundary == "Left")
	{
		lbc = convectionRadiation;
		TLeft = Tamb, hfL = h, eps = emissivity;
	}
	else if (boundary == "Right")
	{
		rbc = convectionRadiation;
		TRight = Tamb, hfR = h, eps = emissivity;
	}
	else if (boundary == "Bottom")
	{
		bbc = convectionRadiation;
		TBottom = Tamb, hfB = h, eps = emissivity;
	}
	else if (boundary == "Top")
	{
		tbc = convectionRadiation;
		TTop = Tamb, hfT = h, eps = emissivity;
	}
	else
	{
		cout << "Wrong boundary specified in Convection - Radiation boundary condition" << endl;
		exit(1);
	}
}
void Shell::setRadiationBC(string boundary, double Tamb, double emissivity)
{
	if (boundary == "Left")
	{
		lbc = ambientRadiation;
		TLeft = Tamb, eps = emissivity;
	}
	else if (boundary == "Right")
	{
		rbc = ambientRadiation;
		TRight = Tamb, eps = emissivity;
	}
	else if (boundary == "Bottom")
	{
		bbc = ambientRadiation;
		TBottom = Tamb, eps = emissivity;
	}
	else if (boundary == "Top")
	{
		tbc = ambientRadiation;
		TTop = Tamb, eps = emissivity;
	}
	else
	{
		cout << "Wrong boundary specified in ambientRadiation boundary condition" << endl;
		exit(1);
	}
}

//interShellBC
void Shell::setRightInterShellContactResistance(double contactResistance){
	if (rightISBC.empty())
	{
		rightISBC.push_back(contactResistance);
	}
	else
	{
		rightISBC[0] = contactResistance;
	}	
}
void Shell::setBottomInterShellContactResistance(double contactResistance)
{
	if (bottomISBC.empty())
	{
		bottomISBC.push_back(contactResistance);
	}
	else
	{
		bottomISBC[0] = contactResistance;
	}	
}
void Shell::setBottomInterShellRadiation(double thisShellEmissivity, double otherShellEmissivity, double radialGap)
{
	if (bottomISBC.empty())
	{
		bottomISBC.push_back(thisShellEmissivity);
		bottomISBC.push_back(otherShellEmissivity);
		bottomISBC.push_back(radialGap);
	}
	else
	{
		cout<<"Bottom Inter-Shell boundary condition given twice"<<endl;
		exit(1);
	}	

}
void Shell::setBottomInterShellRadiationWithConvection(double thisShellEmissivity,
             double otherShellEmissivity, double radialGap, double hf, double Tf)
{
	if (bottomISBC.empty())
	{
		bottomISBC.push_back(thisShellEmissivity);
		bottomISBC.push_back(otherShellEmissivity);
		bottomISBC.push_back(radialGap);
		bottomISBC.push_back(hf);
		bottomISBC.push_back(Tf);
	}
	else
	{
		cout<<"Bottom Inter-Shell boundary condition given twice"<<endl;
		exit(1);
	}	

}


void Shell::setInitialTemp(double initialTemp){
	initTemp = initialTemp;
}
void Shell::setTe(int j, int i, double val){
	//if (val<0){cout<<"A negative value cannot be set to temperature: setTe method"<<endl;exit(1);}
	te[j][i] = val;
}
void Shell::setTe0(int j, int i, double val){
	//if (val<0){cout<<"A negative value cannot be set to temperature: setTe0 method"<<endl;exit(1);}
	te0[j][i] = val;
}
void Shell::setTep(int j, int i, double val)
{
	//if (val<0){cout<<"A negative value cannot be set to temperature: setTep method"<<endl;exit(1);}
	tep[j][i] = val;
}
void Shell::setSp(int j, int i, double val){sp[j][i] = val;}
void Shell::setSc(int j, int i, double val){sc[j][i] = val;}


//getters
bool Shell::isConnected() { return connected; }
bool Shell::getType() { return axi; }
double Shell::getLength() { return Length; }
double Shell::getWidth() { return Width; }
double Shell::getInnerRadius()
{
	if (!axi)
	{
		cout << " No inner radius for a rectangular shell" << endl;
		exit(1);
	}
	return ri;
}
double Shell::getTimeStep() { return dt; }
double Shell::getSimulationTime() { return simTime; }
double Shell::getTe(int j, int i){return te[j][i];}
double Shell::getTk(int j, int i){return tk[j][i];}
double Shell::getRho(int j, int i){return rho[j][i];}
double Shell::getCp(int j, int i){return cp[j][i];}
double Shell::getSp(int j, int i){return sp[j][i];}
double Shell::getSc(int j, int i){return sc[j][i];}
double Shell::getTe0(int j, int i){return te0[j][i];}
double Shell::getY(int j){
	if (j<0) {
		cout<<"Negative index in getY";
		exit(1);
	}
return y[j];
}
double Shell::getDy(int j){
	if (j<0) {
		cout<<"Negative index in getDy";
		exit(1);
	}
return dy[j];
}
double Shell::getDx(int i){
	if (i<0) {
		cout<<"Negative index in getDx";
		exit(1);
	}
return dx[i];
}
int Shell::getMaxIter(){return maxiter;}
double Shell::getRelaxationCoeff(){return re;}
int Shell::getM(){return M;}
int Shell::getN(){return N;}

//Transient solver
void Shell::solveTransient(string fileName, int writeInterval){
	preprocessShell();
	double t{0};
	ofstream outFile;
	//writing geometry data for plotting using python
	//no of shells Length Width N M no of time data
	string geomFile = fileName+".geom"; 
	string tempFile = fileName+".temp";
	outFile.open(geomFile.c_str());
	outFile<<1<<'\t'<<writeInterval<<'\t'<<Length<<'\t'<<Width<<'\t'<<N<<'\t'<<M<<'\t'<<floor(simTime / (dt*writeInterval) )<<endl;
	outFile.close();

	//open tempFile in re write mode to overwrite any existing default file
	outFile.open(tempFile.c_str());
	outFile.close();
	outFile.open(tempFile.c_str(),std::ios_base::app);
	while(t<simTime){
		advanceOneTimeStep();
		int count{0};
		if ((int)(std::round(t / dt)) % writeInterval == 0)
		{
			for (int j = 1; j < M+1; j++)
			{
				for (int i = 1; i < N+1; i++)
				{
					outFile<<te[j][i];
					if (i <N ) outFile<<",";
					count++;
				}
				outFile<<endl;
			}
		}
		t=t+dt;//increment time step
	}

	outFile.close();
}
//Steady solver
void Shell::solveSteady(int maxIter){
	maxiter = maxIter;
	setTimeStep(1e5);
	solveTransient();
}
	
	//preprocessors 
	void Shell::preprocessShell(){
		if (!dx.empty()) {cout<<"This Shell looks to be solved"<<endl;exit(1);}
		populateNodes();
		populateMaterialProperties();
		initialiseField();
	}
		void Shell::populateNodes()
	{
		//populate dx, dy
		for (int i = 0; i < N + 2; i++)
		{
			dx.push_back(Length / double(N));
		}
		for (int j = 0; j < M + 2; j++)
		{
			dy.push_back(Width / double(M));
		}
		dx[0] = 1e-10;
		dx[N + 1] = 1e-10;
		dy[0] = 1e-10;
		dy[M + 1] = 1e-10;

		//populate x & y
		for (int i = 0; i < N + 3; i++)
		{
			if (i == 0)
				x.push_back(0);
			if (i > 0)
				x.push_back(x[i - 1] + dx[i - 1]);
		}
		for (int j = 0; j < M + 3; j++)
		{
			if (j == 0)
			{
				if (!axi)
					y.push_back(0);
				if (axi)
					y.push_back(ri);
			}
			if (j > 0)
				y.push_back(y[j - 1] + dy[j - 1]);
		}
	}
		void Shell::populateMaterialProperties()
	{
		vector<double> tk1d(N + 2, tCond);
		vector<double> cp1d(N + 2, spHeat);
		vector<double> rho1d(N + 2, density);
		for (int i = 0; i < M + 2; i++)
		{
			tk.push_back(tk1d);
			cp.push_back(cp1d);
			rho.push_back(rho1d);
		}
		/* 	Shell::print2dVector(tk);
		Shell::print2dVector(cp);
		Shell::print2dVector(rho); */
	}
		void Shell::initialiseField(){
		vector<double> te1d(N + 2, initTemp);
		vector<double> scsp1d(N + 2, 0);
		for (int i = 0; i < N + 1; i++){
			ta.push_back(0);
			tb.push_back(0);
			tc.push_back(0);
			td.push_back(0);
		}
		for (int i = 0; i < M + 2; i++){
			te.push_back(te1d);
			te0.push_back(te1d);
			tep.push_back(te1d);
			sc.push_back(scsp1d);
			sp.push_back(scsp1d);
		}
	}

	//sub-solver methods
	void Shell::advanceOneTimeStep()
	{
		//for convergence checking
		double maxErr = 1e-10;
		double error = 1.0e-9;

		//defining variables for equation
		double ke{0}, kw{0}, ks{0}, kn{0};
		double de{0}, dw{0}, ds{0}, dn{0};
		double ae{0}, aw{0}, as{0}, an{0};
		double a0{0}, ap{0}, b{0}, vol{0};
		double sae{0}, saw{0}, san{0}, sas{0};
		//variables required for loop
		int iter{0}, iflag = 1;

		//--------------------Outer Loop ---------------------------------------------

		//......................Inner Loop ........................................
		while (iflag == 1)
		{
			//update properties if are dependent on temperature. Currently it is not

			//Enforce boundary conditions

			//add source terms if any. Now nothing
			for (int j = 0; j < M + 2; j++)
			{
				for (int i = 0; i < N + 2; i++)
				{
					sp[j][i] = 0;
					sc[j][i] = 0;
				}
			}
			applyBoundaryConditions();
			// solve -------------------------------------------
			//START marching in y
			for (int j = 1; j < M + 1; j++)
			{
				//START marching in x
				for (int i = 1; i < N + 1; i++)
				{
					if (!axi)
					{
						sae = dy[j];
						saw = sae;
					}
					if (axi)
					{
						sae = dy[j] * (y[j] + y[j + 1]) / 2.0;
						saw = sae;
					}
					ke = tk[j][i] * tk[j][i + 1] * (dx[i] + dx[i + 1]) / (dx[i] * tk[j][i + 1] + dx[i + 1] * tk[j][i]);
					de = 2.0 * ke * sae / (dx[i] + dx[i + 1]);
					ae = de;
					kw = tk[j][i] * tk[j][i - 1] * (dx[i] + dx[i - 1]) / (dx[i] * tk[j][i - 1] + dx[i - 1] * tk[j][i]);
					dw = 2.0 * kw * saw / (dx[i] + dx[i - 1]);
					aw = dw;
					if (!axi)
					{
						san = dx[i];
						sas = san;
					}
					if (axi)
					{
						san = dx[i] * y[j + 1];
						sas = dx[i] * y[j];
					}

					kn = tk[j][i] * tk[j + 1][i] * (dy[j] + dy[j + 1]) / (dy[j] * tk[j + 1][i] + dy[j + 1] * tk[j][i]);
					dn = 2.0 * kn * san / (dy[j] + dy[j + 1]);
					an = dn;
					ks = tk[j][i] * tk[j - 1][i] * (dy[j] + dy[j - 1]) / (dy[j] * tk[j - 1][i] + dy[j - 1] * tk[j][i]);
					ds = 2.0 * ks * sas / (dy[j] + dy[j - 1]);
					as = ds;
					if (!axi)
					{
						vol = dx[i] * dy[j];
					}
					if (axi)
					{
						vol = dx[i] * dy[j] * (y[j] + y[j + 1]) / 2.0;
					}
					a0 = rho[j][i] * cp[j][i] * vol / dt; // =0 for steady state
					ap = ae + aw + an + as + a0 - sp[j][i] * vol;
					b = sc[j][i] * vol + a0 * te0[j][i];
					ta[i] = ap / re;
					tb[i] = ae;
					tc[i] = aw;
					td[i] = b + ap / re * (1 - re) * te[j][i] + an * te[j + 1][i] + as * te[j - 1][i];
					if (i == 1)
						td[i] = td[i] + aw * te[j][0];
					if (i == N)
						td[i] = td[i] + ae * te[j][N + 1];
				} //marching in x ends here
				//END marching in x
				//solve in x direction using tdma
				
				tdma(j);

			} //marching in y ends here
			//END marching in y
			//solve ----------------------------------------------------------------------------------
			//start convergence checking ---------------------
			iflag = checkConvergence(error);
			if (iflag == 1)
				iter++;
			if (iflag == 0)
				cout << "  Solution Converged for one time step. dt = " << dt << " in " << iter << " iterations" << endl;
			//end convergence checking ***********************

			if (iter > maxiter)
			{
				cout << " Iterations need to be inreased. Error in Temp is " << maxErr * 100 << " %" << endl;
				break;
			}

		} //end of inner while loop checking iflag
		//......................Inner Loop ........................................

		for (int i = 0; i < N + 2; i++)
		{
			for (int j = 0; j < M + 2; j++)
			{
				te0[j][i] = te[j][i];
				tep[j][i] = te0[j][i];
			}
		}

		//printTe();

		//end of while loop checking t<simTime
		//--------------------Outer Loop ---------------------------------------------
	}
		void Shell::applyBoundaryConditions(){ //calls inside advanceOneTimeStep
			//tube inlet, which is along r(or y) direction
			double qr = 0;
			for (int j = 0; j < M + 2; j++)
			{
				switch (lbc)
				{
				case constTemp:
					te[j][0] = TLeft;
					break;
				case constHeatFlux:
					te[j][0] = te[j][1] + 0.5 * dx[1] * qLeft / tk[j][1]; //heat flux into the system is +ve
					break;
				case convection:
					te[j][0] = (hfL * TLeft + 2 * te[j][1] * tk[j][1] / dx[1]) / (hfL + 2 * tk[j][1] / dx[1]);
					break;
				case ambientRadiation:
					qr = eps*5.67e-8*(pow(TLeft,4) -pow(te[j][1],4));
					te[j][0] = te[j][1] + 0.5 * dx[1] * qr / tk[j][1]; //heat flux into the system is +ve
					break;
				case convectionRadiation:
					hfL = hfL + eps*5.67e-8*(pow(TLeft,2) + pow(te[j][1],2)* (te[j][1] + TLeft) );
					te[j][0] = (hfL * TLeft + 2 * te[j][1] * tk[j][1] / dx[1]) / (hfL + 2 * tk[j][1] / dx[1]);
					break;
				default:
					te[j][0] = TLeft;
					break;
				}
			}

			//tube outlet, which is along r(or y) direction
			for (int j = 0; j < M + 2; j++)
			{
				switch (rbc)
				{
				case constTemp:
					te[j][N + 1] = TRight;
					break;
				case constHeatFlux:
					te[j][N + 1] = te[j][N] + 0.5 * dx[N] * qRight / tk[j][N]; //heat flux into the system is +ve
					break;
				case convection:
					te[j][N + 1] = (hfR * TRight + 2 * te[j][N] * tk[j][N] / dx[N]) / (hfR + 2 * tk[j][N] / dx[N]);
					break;
				case ambientRadiation:
					qr = eps*5.67e-8*(pow(TRight,4) -pow(te[j][N],4));
					te[j][N + 1] = te[j][N] + 0.5 * dx[N] * qr / tk[j][N]; //heat flux into the system is +ve
					break;
				case convectionRadiation:
					hfR = hfR + eps*5.67e-8*(pow(TRight,2) + pow(te[j][N],2)* (te[j][N] + TRight) );
					te[j][N + 1] = (hfR * TRight + 2 * te[j][N] * tk[j][N] / dx[N]) / (hfR + 2 * tk[j][N] / dx[N]);
					break;
				default:
					te[j][0] = TRight;
					break;
				}
			}

			//inside tube, which is along z(or x) direction
			for (int i = 0; i < N + 2; i++)
			{
				switch (bbc)
				{
				case constTemp:
					te[0][i] = TBottom;
					break;
				case constHeatFlux:
					te[0][i] = te[1][i] + 0.5 * dy[1] * qBottom / tk[1][i]; //heat flux into the system is +ve
					break;
				case convection:
					te[0][i] = (hfB * TBottom + 2 * te[1][i] * tk[1][i] / dy[1]) / (hfB + 2 * tk[1][i] / dy[1]);
					break;
				case ambientRadiation:
					qr = eps*5.67e-8*(pow(TBottom,4) -pow(te[1][i],4));
					te[0][i] = te[1][i] + 0.5 * dy[1] * qr / tk[1][i]; //heat flux into the system is +ve
					break;
				case convectionRadiation:
					hfB = hfB + eps*5.67e-8*(pow(TBottom,2) + pow(te[1][i],2)* (te[1][i] + TBottom) );	
					te[0][i] = (hfB * TBottom + 2 * te[1][i] * tk[1][i] / dy[1]) / (hfB + 2 * tk[1][i] / dy[1]);
					break;
				default:
					te[0][i] = TBottom;
					break;
				}
			}

			//outside tube, which is along z(or x) direction
			for (int i = 0; i < N + 2; i++)
			{
				switch (tbc)
				{
				case constTemp:
					te[M + 1][i] = TTop;
					break;
				case constHeatFlux:
					te[M + 1][i] = te[M][i] + 0.5 * dy[M] * qTop / tk[M][i]; //heat flux into the system is +ve
					break;
				case convection:
					te[M + 1][i] = (hfT * TTop + 2 * te[M][i] * tk[M][i] / dy[M]) / (hfT + 2 * tk[M][i] / dy[M]);
					break;
				case ambientRadiation:
					qr = eps*5.67e-8*(pow(TTop,4) -pow(te[M][i],4));
					te[M + 1][i] = te[M][i] + 0.5 * dy[M] * qr / tk[M][i]; //heat flux into the system is +ve
					break;
				case convectionRadiation:
					hfT = hfT + eps*5.67e-8*(pow(TTop,2) + pow(te[M][i],2)* (te[M][i] + TTop) );		
					te[M + 1][i] = (hfT * TTop + 2 * te[M][i] * tk[M][i] / dy[M]) / (hfT + 2 * tk[M][i] / dy[M]);
					break;
				default:
					te[M + 1][i] = TTop;
					break;
				}
			}
		}
		void Shell::tdma(int j)
		{ //calls inside advanceOneTimeStep

			//print1dVector(ta);		
			double alpha[N + 2]{0}, beta[N + 2]{0}, dum[N + 2]{0};
			for (int i = 0; i < N + 2; i++)
			{
				alpha[i] = 1;
				beta[i] = 1;
				dum[i] = 1;
			}
			beta[1] = tb[1] / ta[1];
			alpha[1] = td[1] / ta[1];
			//forward substitution
			for (int ii = 2; ii < N + 1; ii++)
			{
				beta[ii] = tb[ii] / (ta[ii] - tc[ii] * beta[ii - 1]);
				alpha[ii] = (td[ii] + tc[ii] * alpha[ii - 1]) / (ta[ii] - tc[ii] * beta[ii - 1]);
			}
			//backward substitution
			dum[N] = alpha[N];
			for (int jj = 0; jj < N - 1; jj++)
			{
				int ii = N - 1 - jj;
				dum[ii] = beta[ii] * dum[ii + 1] + alpha[ii];
			}
			//solved value
			for (int i = 1; i < N + 1; i++)
			{
				te[j][i] = dum[i];
				//cout<<"dum["<<i<<"] = "<<dum[i]<<endl;
			}
		}
		int Shell::checkConvergence(double error)
		{ //calls inside advanceOneTimeStep
			double maxErr{1e-10}, errorTe{0};
			int iflag = 1;

			for (int j = 0; j < M + 2; j++)
			{
				for (int i = 0; i < N + 2; i++)
				{
					errorTe = abs(te[j][i] - tep[j][i]) / te[j][i];
					//cout<<"errorTe = "<<errorTe;
					if (errorTe > maxErr)
						maxErr = errorTe;
				}
				//cout<<endl;
			}
			if (maxErr > error)
			{
				for (int i = 0; i < N + 2; i++)
				{
					for (int j = 0; j < M + 2; j++)
					{
						tep[j][i] = te[j][i];
					}
				}
				iflag = 1;
			}
			if (maxErr <= error)
			{
				iflag = 0;
				cout << "MaxErr = " << maxErr;
			}
			return iflag;
		}

	
//printers
void Shell::printDetail()
{
	cout<<"-----------------------------------------------"<<endl;
	if (axi)
		cout << "Cylindrical shell" << endl;
	if (!axi)
		cout << "2D rectangular shell" << endl;

	if (axi)
		cout << "Inner radius(m) :" << ri << endl;
	cout << "Shell Length(m) :" << Length << endl;
	if (!axi)
		cout << "Shell Width(m) :" << Width << endl;
	if (axi)
		cout << "Shell Thickness(m) :" << Width << endl;
	cout << "Material properties  :" << endl<<endl;
	cout << "Thermal Conductivity :" << tCond << endl;
	cout << "Specific Heat :" << spHeat << endl;
	cout << "Density :" << density << endl;

	cout << "MSimulation details  :" << endl<<endl;
	cout << "Simulation time(s) :" << simTime << endl;
	cout << "Number of divisions along y :" << M << endl;
	cout << "Number of divisions along x :" << N << endl;
	cout << "Initial temperature(K) :" << initTemp << endl;
	cout << "Boundary conditions :" << endl;
	switch (lbc)
	{
	case constTemp:
		cout << "Left boundary condition is constTemp of " << TLeft << " K" << endl;
		break;
	case constHeatFlux:
		cout << "Left boundary condition is constHeatFlux of " << qLeft << " watts" << endl;
		break;
	case convection:
		cout << "Left boundary condition is connvection with Tf = " << TLeft << " K and hf = " << hfL << " W/m-k " << endl;
		break;
	case ambientRadiation:
		cout << "Left boundary condition is ambientRadiation with Tamb = " << TLeft << " K and emissivity = " << eps<< endl;
		break;
	case convectionRadiation:
		cout << "Left boundary condition is combined connvection and radiation with Tamb = " << TLeft << " K, hf = " << hfL << " W/m-k  and emissivity = " << eps << endl;
		break;
	default:
		cout << "Invalid boundary condition " << endl;
		break;
	}
	switch (rbc)
	{
	case constTemp:
		cout << "Right boundary condition is constTemp of " << TRight << " K" << endl;
		break;
	case constHeatFlux:
		cout << "Right boundary condition is constHeatFlux of " << qRight << " watts" << endl;
		break;
	case convection:
		cout << "Right boundary condition is connvection with Tf = " << TRight << " K and hf = " << hfR << "W/m-k " << endl;
		break;
	case ambientRadiation:
		cout << "Right boundary condition is ambientRadiation with Tamb = " << TRight << " K and emissivity = " << eps<< endl;
		break;
	case convectionRadiation:
		cout << "Right boundary condition is combined connvection and radiation with Tamb = " << TRight << " K, hf = " << hfR << " W/m-k  and emissivity = " << eps << endl;
	default:
		cout << "Invalid boundary condition " << endl;
		break;
	}
	switch (bbc)
	{
	case constTemp:
		cout << "Bottom boundary condition is constTemp of " << TBottom << " K" << endl;
		break;
	case constHeatFlux:
		cout << "Bottom boundary condition is constHeatFlux of " << qBottom << " watts" << endl;
		break;
	case convection:
		cout << "Bottom boundary condition is connvection with Tf = " << TBottom << " K and hf = " << hfL << "W/m-k " << endl;
		break;
	case ambientRadiation:
		cout << "Bottom boundary condition is ambientRadiation with Tamb = " << TBottom << " K and emissivity = " << eps<< endl;
		break;
	case convectionRadiation:
		cout << "Bottom boundary condition is combined connvection and radiation with Tamb = " << TBottom << " K, hf = " << hfB << " W/m-k  and emissivity = " << eps << endl;
	default:
		cout << "Invalid boundary condition " << endl;
		break;
	}
	switch (tbc)
	{
	case constTemp:
		cout << "Top boundary condition is constTemp of " << TTop << " K" << endl;
		break;
	case constHeatFlux:
		cout << "Top boundary condition is constHeatFlux of " << qTop << " watts" << endl;
		break;
	case convection:
		cout << "Top boundary condition is connvection with Tf = " << TTop << " K and hf = " << hfT << "W/m-k " << endl;
		break;
	case ambientRadiation:
		cout << "Top boundary condition is ambientRadiation with Tamb = " << TTop << " K and emissivity = " << eps<< endl;
		break;
	case convectionRadiation:
		cout << "Top boundary condition is combined connvection and radiation with Tamb = " << TTop << " K, hf = " << hfT << " W/m-k  and emissivity = " << eps << endl;
	default:
		cout << "Invalid boundary condition " << endl;
		break;
	}

	cout << "Maximum iteration :" << maxiter << endl;
	cout << "Relaxation coefficient :" << re << endl;
	cout << "Time step(s) :" << dt << endl;
	cout<<"-----------------------------------------------"<<endl;
}
void Shell::printTe()
{
	for (auto i : te)
	{
		//i is now an 1D vector
		for (auto j : i)
		{
			cout << j << " ";
		}
		cout << endl;
	}
}

//non member function methods

void solveSystem(vector<vector<Shell>> &v, string fileName, int writeInterval){
	//width for horizontal, length for vertical connection
	//Ri for axi
	int rows = v.size();
	int cols = v[0].size();
	int totalShells = rows*cols;
	double Length{0}, Width{0};
	int M{0},N{0};
	double simTime = v[0][0].getSimulationTime();
	double dt = v[0][0].getTimeStep();
	double t=0;

	for (int i = 0; i < v.size(); i++)
	{
		M += v[i][0].getM();
		Width += v[i][0].getWidth();
	}
	for (int j = 0; j < v[0].size(); j++)
	{
		N += v[0][j].getN();
		Length += v[0][j].getLength();
	}

	
	//Display the details
	int shellNo{0};
	for (int i = 0; i < v.size(); i++){
        for (int j = 0; j < v[i].size(); j++){
			shellNo++;
			cout<<"Shell No : "<<shellNo<<endl;
			v[i][j].setSimulationTime(v[0][0].getSimulationTime());
			v[i][j].setTimeStep(v[0][0].getTimeStep());
			if (j>0) {v[i][j].setWidth(v[i][0].getWidth()); }
			if (i>0) {v[i][j].setLength(v[0][j].getLength());}

			v[i][j].printDetail();
			v[i][j].preprocessShell();
			v[i][j].applyBoundaryConditions();
        }   
    }
	shellNo=0;
	applyInterShellBC(v);
	//for file printing
	ofstream outFile;
	string geomFile = fileName+".geom"; 
	string tempFile = fileName+".temp";
	outFile.open(geomFile.c_str());
	outFile<<rows<<'\t'<<writeInterval<<'\t'<<Length<<'\t'<<Width<<'\t'<<N<<'\t'<<M<<'\t'<<floor(simTime / (dt*writeInterval) )<<endl;
	outFile.close();
    outFile.open(tempFile.c_str());
	outFile.close();
	outFile.open(tempFile.c_str(),std::ios_base::app);

	
	while(t<simTime)
	{
		for (int i = 0; i < v.size(); i++)
		{
			for (int j = 0; j < v[i].size(); j++)
			{
				advanceOneTimeStep(v[i][j]);
				v[i][j].applyBoundaryConditions();
			}  
			cout<<endl; 
		}
		applyInterShellBC(v);//bcs set on the internal faces will be bye-passed here
		t = t+dt;
		//outputting file(csv)
		if ((int)(std::round(t / dt)) % writeInterval == 0)//write only at a user specified intervals
		{
			for (int i = v.size()-1; i >=0 ; i--)
			{
				for(int jj=1;jj<v[i][0].getM()+1;jj++)
				{
					for (int j = 0; j < v[i].size(); j++)
					{
						for(int ii=1;ii<v[i][j].getN()+1;ii++)
						{
							outFile<<v[i][j].getTe(jj,ii);
							if (ii == v[i][j].getN() && j == v[i].size()-1){}
							else outFile<<",";
						}
				
					}
					outFile<<endl;
				}
    		}
		}
		

	}//end of while loop checking t<simTime
	outFile.close();
	shellNo=0;
	for (int i = 0; i < v.size(); i++){
        for (int j = 0; j < v[i].size(); j++){
			shellNo++;
			cout<<shellNo<<endl;
			
			v[i][j].printTe();
			v[i][j].printDetail();
        }  
		cout<<endl; 
    }


}
//this will take care (only) the internal bcs 
void applyInterShellBC(vector<vector<Shell>> &v){
	for (int i = 0; i < v.size(); i++){
        for (int j = 0; j < v[i].size(); j++){
			//setting the connection bc
			int M1 = v[i][j].getM();
			int N1 = v[i][j].getN();
			//vertical connection
			if (j<(v[i].size()-1)) {
				int N2 = v[i][j+1].getN();//N of bottom shall can be different but M2=M1
				for (int jj = 0; jj < M1 + 2; jj++){
					double T1 = v[i][j].getTe(jj,N1);
					double T2 = v[i][j+1].getTe(jj,1);//left but 1 column of right shell
					double dx1 = v[i][j].getLength()/N1;
					double dx2 = v[i][j+1].getLength()/N2;
					double dy  = v[i][j].getWidth()/M1;
					double k1 = v[i][j].getTk(jj,N1);
					double k2 = v[i][j+1].getTk(jj,1);
					double Rc = 0;
					if(!v[i][j].rightISBC.empty()) Rc = v[i][j].rightISBC[0];
					double q = (T2 -T1)/(dy * Rc + 0.5 * (dx1/k1 + dx2/k2) );
					v[i][j].setTe(jj, N1+1 , v[i][j].getTe(jj,N1) + 0.5 * dx1 * q / v[i][j].getTk(jj,N1)); 
					v[i][j+1].setTe(jj, 0 , v[i][j+1].getTe(jj,1) - 0.5 * dx2 * q / v[i][j+1].getTk(jj,1));		
				}
			}
			//horizontal connection
			if (i<v.size()-1)
			{	
				int M2 = v[i+1][j].getM();//M of bottom shall can be different but N2=N1
				for (int ii = 0; ii < N1 + 2; ii++){
					double T1 = v[i][j].getTe(1,ii);
					double T2 = v[i+1][j].getTe(M2,ii);//top but 1 row of bottom shell
					double dy1 = v[i][j].getWidth()/M1;
					double dy2 = v[i+1][j].getWidth()/M2;
					double dx  = v[i][j].getLength()/N1;
					double k1 = v[i][j].getTk(1,ii);
					double k2 = v[i+1][j].getTk(M2,ii);
					double Rc = 0;
					double q{0},qconv1{0},qconv2{0};
					if(v[i][j].bottomISBC.size()<2)//non zero contact resistance
					{
						Rc = v[i][j].bottomISBC[0];//other bcs have atleast 2 arguments
						q = (T2 -T1)/(dx * Rc + 0.5 * (dy1/k1 + dy2/k2) );
						v[i][j].setTe(0, ii , v[i][j].getTe(1,ii) + 0.5 * dy1 * q / v[i][j].getTk(1,ii)); 
						v[i+1][j].setTe(M2+1, ii , v[i+1][j].getTe(M2,ii) - 0.5 * dy2 * q / v[i+1][j].getTk(M2,ii));
					}
					if(v[i][j].bottomISBC.size()>2)//radiation(3) & radiation with convection (5)
					{
						//get avg temp of both shells
						double thisShellTemp{0},otherShellTemp{0};
						for (int k = 0; k < N1; k++)
						{
							thisShellTemp += v[i][j].getTe(1,k);
							otherShellTemp += v[i+1][j].getTe(1,k);
						}
						thisShellTemp = thisShellTemp/N1;
						otherShellTemp = otherShellTemp/N1;
						//estimate q
						double eps1 = v[i][j].bottomISBC[0];
						double eps2 = v[i][j].bottomISBC[1];
						double r1 = v[i][j].getInnerRadius();
						double r2 = v[i+1][j].getInnerRadius() + v[i][j].bottomISBC[2];//radial gap added
						q = 5.67e-8*( pow(otherShellTemp,4) - pow(thisShellTemp,4) )
							/(1/eps1 + (1 - eps2)/eps2 * pow(r1/r2,2));
						if(v[i][j].bottomISBC.size()==5)//radiation with convection
						{
							double h = v[i][j].bottomISBC[3];
							double Tf = v[i][j].bottomISBC[4];
							qconv1 = h*(Tf - thisShellTemp);		
							qconv2 = h*(Tf - otherShellTemp);		
						}

						//setTe for both
						v[i][j].setTe(0, ii , v[i][j].getTe(1,ii) + 0.5 * dy1 * (q + qconv1) / v[i][j].getTk(1,ii)); 
						v[i+1][j].setTe(M2+1, ii , v[i+1][j].getTe(M2,ii) - 0.5 * dy2 * (q - qconv2) / v[i+1][j].getTk(M2,ii));
					}
					
				}			
			}

        }   
    }

}


void advanceOneTimeStep(Shell &s)
{
	int M = s.getM();
	int N = s.getN();
		//for convergence checking
		double maxErr = 1e-10;
		double error = 1.0e-9;

		//defining variables for equation
		double ke{0}, kw{0}, ks{0}, kn{0};
		double de{0}, dw{0}, ds{0}, dn{0};
		double ae{0}, aw{0}, as{0}, an{0};
		double a0{0}, ap{0}, b{0}, vol{0};
		double sae{0}, saw{0}, san{0}, sas{0};

		//for  tdma
		vector<double> ta, tb, tc, td;
		for (int i = 0; i < N + 1; i++)
		{
			ta.push_back(0);
			tb.push_back(0);
			tc.push_back(0);
			td.push_back(0);
		}
		//variables required for loop
		int iter{0}, iflag = 1;

		bool axi=s.getType();
		double re = s.getRelaxationCoeff();
		while (iflag == 1)
		{
			//add source terms if any. Now nothing
			for (int j = 0; j < M + 2; j++)
			{
				for (int i = 0; i < N + 2; i++){
					s.setSp(j,i,0);
					s.setSc(j,i,0);
				}
			}

			//START marching in y
			for (int j = 1; j < M + 1; j++)
			{
				double Yj = s.getY(j); 
				double Yjp1 = s.getY(j+1);
				double dyj = s.getDy(j); 
				double dyjp1 = s.getDy(j+1);
				double dyjm1 = s.getDy(j-1);
				//START marching in x
				for (int i = 1; i < N + 1; i++)
				{
					double dxi = s.getDx(i);
					double dxip1 = s.getDx(i+1);
					double dxim1 = s.getDx(i-1);
					double Tkji = s.getTk(j,i);
					double Tkjip1 = s.getTk(j,i+1);
					double Tkjim1 = s.getTk(j,i-1);
					double Tkjp1i = s.getTk(j+1,i);
					double Tkjm1i = s.getTk(j-1,i);
					
					if (!axi)
					{
						sae = dyj;
						saw = sae;
					}
					if (axi)
					{
						sae = dyj * (Yj + Yjp1) / 2.0;
						saw = sae;
					}
					ke = Tkji * Tkjip1 * (dxi + dxip1) / (dxi * Tkjip1 + dxip1 * Tkji);
					de = 2.0 * ke * sae / (dxi + dxip1);
					ae = de;
					kw = Tkji * Tkjim1 * (dxi + dxim1) / (dxi * Tkjim1 + dxim1 * Tkji);
					dw = 2.0 * kw * saw / (dxi + dxim1);
					aw = dw;
					if (!axi)
					{
						san = dxi;
						sas = san;
					}
					if (axi)
					{
						san = dxi * Yjp1;
						sas = dxi * Yj;
					}

					kn = Tkji * Tkjp1i * (dyj + dyjp1) / (dyj * Tkjp1i + dyjp1 * Tkji);
					dn = 2.0 * kn * san / (dyj + dyjp1);
					an = dn;
					ks = Tkji * Tkjm1i * (dyj + dyjm1) / (dyj * Tkjm1i + dyjm1 * Tkji);
					ds = 2.0 * ks * sas / (dyj + dyjm1);
					as = ds;
					if (!axi)
					{
						vol = dxi * dyj;
					}
					if (axi)
					{
						vol = dxi * dyj * (Yj + Yjp1) / 2.0;
					}
					a0 = s.getRho(j,i) * s.getCp(j,i) * vol / s.getTimeStep(); // =0 for steady state
					ap = ae + aw + an + as + a0 - s.getSp(j,i) * vol;
					b = s.getSc(j,i) * vol + a0 * s.getTe0(j,i);
					ta[i] = ap / re;
					tb[i] = ae;
					tc[i] = aw;
					td[i] = b + ap / re * (1 - re) * s.getTe(j,i) + an * s.getTe(j+1,i) + as * s.getTe(j-1,i);
					if (i == 1)
						td[i] = td[i] + aw * s.getTe(j,0);
					if (i == N)
						td[i] = td[i] + ae * s.getTe(j,N+1);
				} //marching in x ends here
				//END marching in x
				//solve in x direction using tdma
				
				tdma(s, j, N, ta, tb, tc, td);

			} //marching in y ends here
			//END marching in y
			//solve ----------------------------------------------------------------------------------
			//start convergence checking ---------------------
			iflag = s.checkConvergence(error);
			if (iflag == 1)
				iter++;
			if (iflag == 0)
				cout << "  Solution Converged for one time step. dt = " << s.getTimeStep() << " in " << iter << " iterations" << endl;
			//end convergence checking ***********************

			if (iter > s.getMaxIter())
			{
				cout << " Iterations need to be inreased. Error in Temp is " << maxErr * 100 << " %" << endl;
				break;
			}

		} //end of inner while loop checking iflag

		for (int i = 0; i < N + 2; i++)
		{
			for (int j = 0; j < M + 2; j++)
			{
				s.setTe0(j,i, s.getTe(j,i));//te0[j][i] = te[j][i];
				s.setTep(j,i, s.getTe0(j,i));//tep[j][i] = te0[j][i];
			}
		}
	}
	

void tdma(Shell &s, int &j, int &N, vector<double> &ta, vector<double> &tb, vector<double> &tc, vector<double> &td )
{ //calls inside advanceOneTimeStep
	double alpha[N + 2]{0}, beta[N + 2]{0}, dum[N + 2]{0};

	for (int i = 0; i < N + 2; i++)
	{
		alpha[i] = 1;
		beta[i] = 1;
		dum[i] = 1;
	}
	beta[1] = tb[1] / ta[1];
	alpha[1] = td[1] / ta[1];
	//forward substitution
	for (int ii = 2; ii < N + 1; ii++)
	{
		beta[ii] = tb[ii] / (ta[ii] - tc[ii] * beta[ii - 1]);
		alpha[ii] = (td[ii] + tc[ii] * alpha[ii - 1]) / (ta[ii] - tc[ii] * beta[ii - 1]);
	}
	//backward substitution
	dum[N] = alpha[N];
	for (int jj = 0; jj < N - 1; jj++)
	{
		int ii = N - 1 - jj;
		dum[ii] = beta[ii] * dum[ii + 1] + alpha[ii];
	}
	//solved value
	for (int i = 1; i < N + 1; i++)
	{
		s.setTe(j,i,dum[i]);
	}
}	
void print2dVector(vector<vector<double>> const &v)
{
	for (auto i : v)
	{
		//i is now an 1D vector
		for (auto j : i)
		{
			cout << j << " ";
		}
		cout << endl;
	}
}
void print1dVector(vector<double> const &v)
{
	for (auto i : v)
	{
		cout << i << " ";
	}
	cout << endl;
}
