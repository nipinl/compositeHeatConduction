#ifndef SHELL_H_
#define SHEL_H_
#include<iostream>
#include<vector>
using namespace std;
class Shell{
    protected:
        bool axi;//axi symmetric case; like tube
        bool connected;
        int M;//Number of divisions in y direction
        int N;//Number of divisions in x direction
        double ri;//Inner radius of the tube
        double Length;//Length of the cylinder
        double Width;//thickness of the tube

        //Material properties
        double tCond; 
        double spHeat;
        double density;

        //Initial condition
        double initTemp;//Initial uniform temperature of the tube

        //Boundary conditions
        enum bc {constTemp, constHeatFlux, convection};
        bc lbc,rbc,bbc,tbc;

        // variables for bc
        double TLeft,TRight,TBottom,TTop;
        // heat flux bc
        double qLeft,qRight,qBottom,qTop;
        //convection
        double TfLeft,TfRight,TfBottom,TfTop;
        double hfL,hfR,hfB,hfT;

        //variables for simulation control
        int maxiter;
        double re;//relaxation factor
        double simTime;
        double dt;
        
        //variables for solution
        vector<double> x,y,dx,dy;
        vector<vector<double>> tk,cp,rho;
        vector<vector<double>> te0,tep,te;
        vector<vector<double>> sc,sp;
        vector<double> ta,tb,tc,td;

    public:
        Shell();
        //setters
        void setConnected();//change the connection flag to true
        void setGeometry(double length=0.04, double width = 0.02, bool cylindrical=true, double innerRadius = 0.1);
        void setLength(double length=0.04);
        void setWidth(double width=0.04);
        void setInnerRadius(double innerRadius = 0.1);
        void setSimulationTime(double simtime=20);
        void setTimeStep( double delt=0.1);
        void setMaterialProperties(double Thermal_cond=16,double Cp = 500,double Density = 8000);
        void setInitialTemp(double initialTemp = 300);
        void setConstantTempBC(string boundary= "Left", double Temp = 300);
        void setConstantHeatfluxBC(string boundary= "Bottom", double q = 0);
        void setConvectionBC(string boundary= "Top", double Tinf=300, double h = 10);
        void setTe(int j, int i, double val);
        void setTe0(int j, int i, double val);
        void setTep(int j, int i, double val);
        void setSp(int j, int i, double val);
        void setSc(int j, int i, double val);
                
        //getters
        bool isConnected();//checks if shell is connected to other.
        int getM();
        int getN();
        bool getType();
        double getLength();
        double getWidth();
        double getInnerRadius();
        double getSimulationTime();
        double getTimeStep();
        double getTk(int j, int i);
        double getTe(int j, int i);
        double getTe0(int j, int i);
        double getRho(int j, int i);
        double getCp(int j, int i);
        double getSp(int j, int i);
        double getSc(int j, int i);
        double getY(int j);
        double getDy(int j);
        double getDx(int j);
        int getMaxIter();
        double getRelaxationCoeff();
        

        //Transient solver
        void solveTransient(double simulationTime = 20, double delt = 0.1);
        //Steady solver
        void solveSteady(int maxIter=1000);
        
       
        //preprocessors
        void preprocessShell();
            void populateNodes();
            void populateMaterialProperties();
            void initialiseField();

        //sub-solver methods

        void advanceOneTimeStep();
        void applyBoundaryConditions();
        void tdma(int);
        int checkConvergence(double);
        


        //printers
        void printDetail();
        void printTe();
};

void connectShells(Shell& s1,Shell& s2, double gap=0.01, double interfaceResistance=0);//assume outer(top) of first is connected to inner(bottom) of second
void tdma(Shell &s, int &j, int &N, vector<double> &ta, vector<double> &tb, vector<double> &tc, vector<double> &td );
void solveSystem(vector<vector<Shell>> &v);
void setConnectionBC(vector<vector<Shell>> &v);
void print2dVector(vector<vector<double>> const &v);


#endif
