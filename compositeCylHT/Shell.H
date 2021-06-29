#ifndef SHELL_H_
#define SHEL_H_
#include<iostream>
#include<vector>
using namespace std;
class Shell{
    protected:
        //double x;
        //double y;
        bool axi;//axi symmetric case; like tube
        const int M;//Number of divisions in y direction
        int N;//Number of divisions in x direction
        double ri;//Inner radius of the tube
        double Length;//Length of the cylinder
        double Width;//thickness of the tube
        //Material properties
        double tCond; 
        double spHeat;
        double density;
        double simTime;//in seconds

        //Initial condition
        double initTemp;//Initial uniform temperature of the tube

        //Boundary conditions
        enum bc {constTemp, constHeatFlux, convection};
        bc lbc,rbc,bbc,tbc;

        // temp bc
        double TLeft,TRight,TBottom,TTop;

        // heat flux bc
        double qLeft,qRight,qBottom,qTop;

        //convection
        double TfLeft,TfRight,TfBottom,TfTop;
        double hfL,hfR,hfB,hfT;

        int maxiter;
        double re;//relaxation factor
        double dt;//* 1e15;//dt >1e10 for steady state

        vector<double> x,y,dx,dy;
        vector<vector<double>> tk,cp,rho;
        vector<vector<double>> te0,tep,te;
        vector<vector<double>> sc,sp;
        vector<double> ta,tb,tc,td;
        //{0}, y[M+3]{0}, dx[N+2]{0}, dy[M+2]{0};




    public:
        Shell();
        //Shell(double real, double img);
        //Shell(const Shell &C);
        const Shell & operator=(const Shell &other);

        double getReal()const {return M;}
        double getImg()const {return N;}
        void populateNodes();
        void populateMaterialProperties();
        void initialiseField();
        void solveIt();
        void applyBoundaryConditions();
        void tdma(int);
        int checkConvergence(double);
        void print2dVector(vector<vector<double>> const &v);
        void printDetail();
};
ostream &operator<<(ostream& out,const Shell& C);


#endif
