/*
first draft of parallel code
*/

#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>

using namespace std;
double * FindPair_brute(double * x, int N);
// function to decompose particles to all processes
void decompose_particles(int N, int rank, int size, double * loc, 
                         int & N_proc, double ** substart);

class SPH_parallel{

private:
    int rank = 0;               // rank of the SPH process
    double k = 2000;            // Gas constant
    double rho_0 = 1000;        // Resting density
    double miu = 1.0;           // Viscosity
    double g = 9.81;            // Acceleration of gravity
    double h = 0.01;            // Radius of influence
    double e = 0.5;             // Coefficient of Restitution
    double m = 1.0;             // mass
    double domain[2] = {0, 1};  // domain
    int N;                      // total No. of particles
    int N_proc;                 // No. of particles for 1 process
    double dt = 1E-4;           // time step 


    double * x;                // cordinate of particles
    double * v;                // velocity of particles
    double * r;                // distance of particles
    double * q;                // norm of r_ijs
    double * p;                // presure of paritcles
    double * Fp;              // Presure Force
    double * Fv;              // Viscous Force
    double * Fg;              // Gravity Force
    double * rho;              // Density of particle
    double * phi_d;            // phi_d for calculating density
    double * phi_p;            // \nabla phi_p for calculating pressure
    double * phi_v;            // \nabla^2 phi_v for calculating velocity
    double * a;                // acceleration

public:
    // constructor
    SPH_parallel(int N, int rank, int N_proc);
    // destructor
    ~SPH_parallel();
    // input location
    void inputLocation(double * loc);
    // rho calculation
    void calRho(double * loc);
    // Pressure Force Calculation
    void calPre();
    // Viscous Force Calculation
    void calVis();
    // Gravity Force Calculation
    void calGra();
    // time integration;
    void timeInte();
    // write result;
    void writeResult();
    // scale m and recalculate rho
    void scaleRecal();
    // print matrix
    void printMatrix(double * Matrix, int m, int n){
        for(int i = 0; i < n; i++){
            for(int j = 0; j < m; j++){
                cout << Matrix[i * m + j] << " ";
            }
            cout << endl;
        }
    }
};




