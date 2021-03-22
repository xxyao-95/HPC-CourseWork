/*
Parallel SPH header file made by xiyao
this file contains all the function declaration for the 
SPH_parallel class
*/
// header guard
# pragma once
// includes
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "datastructure/DataStructure.h"
#include "mpi.h"

using namespace std;

// SPH_parallel class
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
    double Ek = 0;              // kinetic energy
    double Ep = 0;              // potential energy
    double E_total = 0;         // total energy


    double * x;                // cordinate of particles
    double * v;                // velocity of particles
    double * r;                // distance of particles
    double * q;                // norm of r_ijs
    double * p;                // presure of paritcles
    double * Fp;               // Presure Force
    double * Fv;               // Viscous Force
    double * Fg;               // Gravity Force
    double * rho;              // Density of particle
    double * phi_d;            // phi_d for calculating density
    double * phi_p;            // \nabla phi_p for calculating pressure
    double * phi_v;            // \nabla^2 phi_v for calculating velocity
    double * a;                // acceleration
    double * rho_global;       // rho of all processes
    double * p_global;         // p of all processes
    double * x_global;         // x of all processes
    double * v_global;         // v of all processes
    
    DecomposeInfo * info;       // hold information of world info
    vector<int> coordQ;         // coordinate of q that is < 1

public:
    // constructor
    SPH_parallel(int N, int size, int rank);
    // destructor
    ~SPH_parallel();
    // input location
    void inputLocation(double * loc);
    // rho calculation
    void calRho();
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
    // send/gather locations of all particle
    void sendRecvLoc();
    // calculate rij and q
    void calRijQ();
    // calculate energy
    void calEnergy();
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




