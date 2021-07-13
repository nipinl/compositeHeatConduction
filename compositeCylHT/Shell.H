#ifndef SHELL_H_
#define SHEL_H_
#include <iostream>
#include<fstream>
#include <vector>
#include <math.h>
#include<algorithm>
using namespace std;
class Shell
{
protected:
    bool axi; //axi symmetric case; like tube
    bool connected;
    int M;         //Number of divisions in y direction
    int N;         //Number of divisions in x direction
    double ri;     //Inner radius of the tube
    double Length; //Length of the cylinder
    double Width;  //thickness of the tube

    //Material properties
    double tCondx,tCondy;
    double spHeat;
    double density;

    //Initial condition
    double initTemp; //Initial uniform temperature of the tube

    //Boundary conditions
    enum bc
    {
        constTemp,
        constHeatFlux,
        convection,
        ambientRadiation,
        convectionRadiation
    };
    bc lbc, rbc, bbc, tbc;

    // variables for bc
    double TLeft, TRight, TBottom, TTop;
    // heat flux bc
    double qLeft, qRight, qBottom, qTop;
    //convection
    //double TLeft,TRight,TBottom,TTop;
    double hfL, hfR, hfB, hfT;
    //ambientRadiation
    double eps;
    

    //variables for simulation control
    int maxiter;
    double re; //relaxation factor
    double simTime;
    double dt;

    //variables for solution
    vector<double> x, y, dx, dy;
    vector<vector<double>> tkx, tky, cp, rho;
    vector<pair<double, double> > tkxTable, tkyTable;
    vector<vector<double>> te0, tep, te;
    vector<vector<double>> sc, sp;
    vector<double> ta, tb, tc, td;
public:
    //connection BC
    vector<double> rightISBC,bottomISBC;//otherwise a lot of getters are needed :)

public:
    Shell();
    //setters
    void setConnected(); //change the connection flag to true
    void cylindrical( bool cylindrical = true);
    void setLength(double length = 0.04);
    void setWidth(double width = 0.04);
    void setInnerRadius(double innerRadius = 0.1);
    void setSimulationTime(double simtime = 20);
    void setTimeStep(double delt = 0.1);
    void setThermalConductivity(double Thermal_cond=16);
    void setThermalConductivity(double Thermal_cond_x=16,double Thermal_cond_y=16);
    void setVariableThermalConductivity(int direction, const vector<double> &temp, const vector<double> &tk);
    void setHeatCapacity(double Cp = 500);
    void setDensity(double Density = 8000);
    void setInitialTemp(double initialTemp = 300);
    void setConstantTempBC(string boundary = "Left", double Temp = 300);
    void setConstantHeatfluxBC(string boundary = "Bottom", double q = 0);
    void setConvectionBC(string boundary = "Top", double Tinf = 300, double h = 10);
    void setRadiationBC(string boundary = "Top", double Tamb = 300, double emissivity = 0.8);
    void setConvectionRadiationBC(string boundary = "Top", double Tamb = 300, double h = 10, double emissivity = 0.8);
    void setRightInterShellContactResistance(double contactResistance = 0);
    void setBottomInterShellContactResistance(double contactResistance = 0);
    void setBottomInterShellRadiation(double thisShellEmissivity = 0.8, double otherShellEmissivity = 0.8, double radialGap = 0.1);
    void setBottomInterShellRadiationWithConvection(double thisShellEmissivity = 0.8,
                                                    double otherShellEmissivity = 0.8, double radialGap = 0.1, double hf = 10, double Tf = 300);
    void setTe(int j, int i, double val);
    void setTe0(int j, int i, double val);
    void setTep(int j, int i, double val);
    void setSp(int j, int i, double val);
    void setSc(int j, int i, double val);
    

    //getters
    bool isConnected(); //checks if shell is connected to other.
    int getM();
    int getN();
    bool getType();
    double getLength();
    double getWidth();
    double getInnerRadius();
    double getSimulationTime();
    double getTimeStep();
    double getTkx(int j, int i);
    double getTky(int j, int i);
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
    void solveTransient(string fileName= "s", int writeInterval = 1);
    //Steady solver
    void solveSteady(int maxIter = 1000);
 

    //preprocessors
    void preprocessShell();
    void populateNodes();
    void populateMaterialProperties();
    void initialiseField();

    //sub-solver methods

    void advanceOneTimeStep();
    void applyBoundaryConditions();
    void updateThermalConductivity();
    void tdma(int);
    int checkConvergence(double);

    //printers
    void printDetail();
    void printTe();
};

void connectShells(Shell &s1, Shell &s2, double gap = 0.01, double interfaceResistance = 0); //assume outer(top) of first is connected to inner(bottom) of second
void tdma(Shell &s, int &j, int &N, vector<double> &ta, vector<double> &tb, vector<double> &tc, vector<double> &td);
void solveTransient(vector<vector<Shell>> &v, string fileName= "s", int writeInterval = 1);
void solveTransient(Shell & s, string fileName="s", int writeInterval=1);
void solveSteady(vector<vector<Shell>> &v, string fileName= "s");
void solveSteady(Shell & s, string fileName="s");
void advanceOneTimeStep(Shell &s);
void applyInterShellBC(vector<vector<Shell>> &v);
double interpolate(double x,vector<pair<double, double> > & table );
void print2dVector(vector<vector<double>> const &v);
void print1dVector(vector<double> const &v);

#endif
