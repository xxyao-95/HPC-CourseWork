#include <iostream>
#include <cmath>
using namespace std;

#define F77NAME(x) x##_

extern "C"{
    // level 1
    double F77NAME(ddot) (const int & n,
                          const double *x, const int & incx,
                          const double *y, const int & incy);
    
    void F77NAME(dcopy) (const int & N,
                         double * x,
                         const int & incx,
                         double * y,
                         const int & incy);
    
    void F77NAME(daxpy)(const int & N,
                        const double & Alpha,
                        double * x,
                        const int & incx,
                        double * y,
                        const int & incy);

    double F77NAME(dnrm2)(const int & N,
                          double * x,
                          const int & incx);
    
    double F77NAME(dasum)(const int & N,
                          double * x,
                          const int & incx);

    double F77NAME(dscal)(const int & N,
                          const double & alpha,
                          double * x,
                          const int & incx);

    
        // level 2
    void F77NAME(dsymv) (const char & UPLO, // specify upper or lower triangular part to be referenced
                         const int & N, const double & Alpha,
                         double * A, const int & Lda,
                         double * x, const int & incx,
                         const double & Beta,
                         double * y, const int & incy);


}



class SPH_serial{

private:
    double k = 2000;            // Gas constant
    double rho_0 = 1000;        // Resting density
    double miu = 1.0;           // Viscosity
    double g = 9.81;            // Acceleration of gravity
    double h = 0.01;            // Radius of influence
    double e = 0.5;             // Coefficient of Restitution
    double m = 1.0;             // mass
    double domain[2] = {0, 1};  // domain
    int N;                      // No. of particles
    double dt = 1E-4;           // time step       

    double ** x;                // cordinate of particles
    double ** v;                // velocity of particles
    double ** r;                // distance of particles
    double * q;                 // norm of r_ijs
    double * p;                 // presure of paritcles
    double ** F_p;              // Presure Force
    double ** F_v;              // Viscous Force
    double ** F_g;              // Gravity Force
    double * rho_i;             // Density of particle
    double * phi_d;             // phi_d for calculating density
    double ** phi_p;            // \nabla phi_p for calculating pressure
    double * phi_v;             // \nabla^2 phi_v for calculating velocity
    double ** a;                // acceleration

public:
    // constructor
    SPH_serial(int N);
    // input starting location
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
    // destructor();
    ~SPH_serial();
    // print matrix
    void printMatrix(double * Matrix, int n){
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                cout << Matrix[i * n + j] << " ";
            }
            cout << endl;
        }
    }
    // scale m and recalculate rho
    void scaleRecal();


};

// constructor
SPH_serial::SPH_serial(int N){
    // N particles, initialise all parameters based on number of particles
    x = new double * [N];
    r = new double * [N * N];
    q = new double [N * N];
    phi_d = new double [N * N];
    rho_i = new double [N];
    p = new double [N];
    phi_p = new double * [N * N];
    phi_v = new double [N * N];

    F_p = new double * [N];
    F_v = new double * [N];
    F_g = new double * [N];

    a = new double * [N]; 
    v = new double * [N];
    
    for (int i = 0; i < N; i++){
        x[i] = new double[2];
        v[i] = new double[2];
        F_p[i] = new double[2];
        F_v[i] = new double[2];
        F_g[i] = new double [2];
        a[i] = new double[2];
        v[i][0] = 0;
        v[i][1] = 0;
    }

    for (int i = 0; i < N*N; i++){
        r[i] = new double[2];
        phi_p[i] = new double[2]; 
    }

}

// destructor
SPH_serial::~SPH_serial(){


    for(int i=0; i<N; i++){
        delete [] x[i];
    }
    delete [] x;

    for(int i=0; i<N; i++){
        delete [] v[i];
    }
    delete [] v;

    for(int i=0; i < N*N; i++){
        delete [] r[i];
    }
    delete [] r;
    delete [] q;
    delete [] p;

    for(int i=0; i<N; i++){
        delete [] F_p[i];
    }
    delete [] F_p;

    for(int i=0; i<N; i++){
        delete [] F_v[i];
    }
    delete [] F_v;

    for(int i=0; i<N; i++){
        delete [] F_g[i];
    }
    delete [] F_g;

    delete [] rho_i;

    delete [] phi_d;

    for(int i=0; i < N*N; i++){
        delete [] phi_p[i];
    }
    delete [] phi_p;
    delete [] phi_v;
    for(int i=0; i<N; i++){
        delete [] a[i];
    }
    delete [] a;

}



// input location
void SPH_serial::inputLocation(double * loc){
    // read locations from array loc(col major)
    for(int i = 0; i < N; i++){
        x[i][0] = loc[2*i];
        x[i][1] = loc[2*i + 1];
        cout << x[i][0] << endl;
        cout << x[i][0] << endl;
    }
}

// calculate density for each particle
void SPH_serial::calRho(){

    // fill r matrix
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            double * temp = x[j];
            F77NAME(daxpy) (2, -1.0, x[i], 1, temp, 1);
            F77NAME(dcopy) (2, temp, 1, r[i*N + j], 1);
        }
    }

    // calculate Norm to fill q matrix
    for(int i = 0; i < N; i++){
        for (int j = i; j < N; j++){
            q[i*N + j] = F77NAME(dnrm2)(2, r[i*N + j], 1);
        }
    }


    // fill phi_d matrix
    for(int i = 0; i < N; i++){
        for (int j = i; j < N; j++){
            if (q[i*N + j] < 1){
                phi_d[i*N + j] = pow(1.0 - pow(q[i*N + j],2), 3.0) * 4/(M_PI * h * h) ;
            }else{
                phi_d[i*N + j] = 0;
            }
        }
    }

    // calculate rho
    double * temp = new double[N];
    fill(temp, temp + N, 1.0);
    F77NAME(dsymv)('L', N, m, phi_d, N, temp, 1, 0.0, rho_i, 1);
    delete [] temp;
    for(int i=0; i < N; i++){
        cout << "..." << endl;
        cout<<rho_i[i]<<endl;
        cout << "........." << endl;
    }

}

// scale m and recalculate rho
void SPH_serial::scaleRecal(){
    m = N * rho_0 / F77NAME(dasum)(N, rho_i, 1);
    double * temp = new double[N];
    fill(temp, temp + N, 1.0);
    F77NAME(dsymv)('L', N, m, phi_d, N, temp, 1, 0.0, rho_i, 1);
    delete [] temp;
    for(int i=0; i<N; i++){
        cout<<rho_i[i]<<endl;
        cout << "........." << endl;
    }

}

// calculate pressure force
void SPH_serial::calPre(){
    // calculate p
    fill(p, p + N, -k * rho_0);
    F77NAME(daxpy)(N, k, rho_i, 1, p, 1);
    // calculate phi_p
    for (int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            F77NAME(dcopy)(2, r[i*N +j], 1, phi_p[i*N+j], 1);
            double scal_fac = 0.0;
            if (q[i*N + j] < 1 && i != j){
                scal_fac = -30 / M_PI / pow(h,3.0);
                scal_fac *= pow((1-q[i*N + j]),2)/q[i*N + j];
            }
            F77NAME(dscal)(N, scal_fac, phi_p[i*N + j], 1);
        }
    }
    // calculate F_p
    for (int i=0; i<N; i++){
        F_p[i][0] = 0;
        F_p[i][1] = 0;

        for (int j = 0; j<N; j++){
            double scale_fac = -(m/rho_i[j]) * (p[i] + p[j])/2; 
            F77NAME(dscal)(N, scale_fac, phi_p[i*N + j], 1);
            F77NAME(daxpy)(N, 1.0, phi_p[i*N + j], 1, F_p[i], 1);
        }
        cout << F_p[i][0] << endl;
        cout << F_p[i][1] << endl;
        cout << "......." << endl;
    }
}

// calculate Viscous Force







int main(int argc, char const *argv[])
{
    SPH_serial test = SPH_serial(1);
    double loc[2] = {0.5, 0.5};
    test.inputLocation(loc);
    test.calRho();
    test.scaleRecal();
    test.calPre();
    return 0;
}
