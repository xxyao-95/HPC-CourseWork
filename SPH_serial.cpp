#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>

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
    this->N = N;
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
    }
}

// calculate density for each particle
void SPH_serial::calRho(){

    // fill r matrix, r is N*N matrix of [delta_x, delta_y]
    for(int i = 0; i < N; i++){ // col 
        for (int j = 0; j < N; j++){ // row
            F77NAME(dcopy)(2, x[j], 1, r[i*N + j], 1); // r = xj
            F77NAME(daxpy) (2, -1.0, x[i], 1, r[i*N + j], 1); // r = xj - xi
        }
    }

    // calculate Norm to fill q matrix
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            q[i*N + j] = F77NAME(dnrm2)(2, r[i*N + j], 1)/h; // q = norm2(r)/h
        }
    }

    // printMatrix(q, N);
    // cout << "q matrix" << endl;
    // fill phi_d matrix
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (q[i*N + j] < 1){
                phi_d[i*N + j] = pow(1.0 - pow(q[i*N + j],2), 3.0) * 4/(M_PI * h * h) ; // phi_d = 4/(pi*h^2) * (1-q^2)^3
            }else{
                phi_d[i*N + j] = 0;
            }
        }
    }

    // calculate rho
    double * temp = new double[N]; 
    fill(temp, temp + N, 1.0); // temp = [1,1,...1]
    F77NAME(dsymv)('L', N, m, phi_d, N, temp, 1, 0.0, rho_i, 1); // m * phi_d * temp
    delete [] temp;
    // for(int i=0; i<N; i++){
    //     cout<< rho_i[i] <<endl;
    //     cout << "Density calculated" << endl;
    // }

}

// scale m and recalculate rho
void SPH_serial::scaleRecal(){
    m = N * rho_0 / F77NAME(dasum)(N, rho_i, 1);
    double * temp = new double[N];
    fill(temp, temp + N, 1.0);
    F77NAME(dsymv)('L', N, m, phi_d, N, temp, 1, 0.0, rho_i, 1);
    delete [] temp;
    // for(int i=0; i<N; i++){
    //     cout<< rho_i[i] <<endl;
    //     cout << "Density calculated" << endl;
    // }

}

// calculate pressure force
void SPH_serial::calPre(){
    // calculate p
    fill(p, p + N, -k * rho_0);
    F77NAME(daxpy)(N, k, rho_i, 1, p, 1);

    // calculate phi_p
    for (int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            F77NAME(dcopy)(2, r[i*N +j], 1, phi_p[i*N+j], 1); // phi_p[i*N + j] = r[i*N +j]
            double scal_fac = 0.0;
            if (q[i*N + j] < 1 && i != j){
                scal_fac = -30 / M_PI / pow(h,3.0); 
                scal_fac *= pow((1-q[i*N + j]),2)/q[i*N + j]; // scal_fac = -30/(pi*h^3) * (1-q)^2/q
            }
            F77NAME(dscal)(2, scal_fac, phi_p[i*N + j], 1);
            // cout << phi_p[i*N + j][0] << " ";
            // cout << phi_p[i*N + j][1] << endl;
        }
    }
    // cout << "phi_p calculated" << endl;


    // calculate F_p
    for (int i=0; i<N; i++){
        F_p[i][0] = 0;
        F_p[i][1] = 0;

        for (int j = 0; j<N; j++){
            double scale_fac = -(m/rho_i[j]) * (p[i] + p[j])/2; 
            F77NAME(dscal)(2, scale_fac, phi_p[j*N + i], 1);
            F77NAME(daxpy)(2, 1.0, phi_p[j*N + i], 1, F_p[i], 1);
        }
//        cout << F_p[i][0] << endl;
//        cout << F_p[i][1] << endl;
//        cout << "pressure force calculated" << endl;
    }
}

// calculate Viscous Force
void SPH_serial::calVis(){

    // calculate Vij
    double ** vij = new double * [N * N];
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            vij[i*N + j] = new double[2];
            F77NAME(dcopy)(2, v[j], 1, vij[i*N + j], 1);
            F77NAME(daxpy) (2, -1.0, v[i], 1, vij[i*N + j], 1);
        }
    }

    // calculate phi_v
    for(int i = 0; i < N; i ++){
        for(int j = 0; j < N; j++ ){
            if (q[i*N + j] < 1 && i !=j ){
                phi_v[i*N + j] = (1 - q[i*N + j]) * 40/M_PI/pow(h, 4);
            }else{
                phi_v[i*N + j] = 0;
            }
        }
    }


    //calculate F_v
    for (int i=0; i<N; i++){
        F_v[i][0] = 0;
        F_v[i][1] = 0;

        for (int j = 0; j<N; j++){
            double scale_fac = -miu * (m/rho_i[j]) * phi_v[i*N + j]; 
            F77NAME(dscal)(2, scale_fac, vij[j*N + i], 1);
            F77NAME(daxpy)(2, 1.0, vij[j*N + i], 1, F_v[i], 1);
        }
        // cout << F_v[i][0] << endl;
        // cout << F_v[i][1] << endl;
        // cout << "viscous force calaulated" << endl;
    }

    for (int i=0; i< N*N; i++){
        delete [] vij[i];
        
    }
    delete [] vij;

}


// calculated gravity force
void SPH_serial::calGra(){
    for (int i=0; i<N; i++){
        F_g[i][0] = 0;
        F_g[i][1] = -rho_i[i] * g;

    }

}

// time integration
void SPH_serial::timeInte(){
    int t = 1;
    ofstream Fout;
    Fout.open("output_4_noise.txt");
    for (int i=0; i<N; i++){ 
            Fout << setw(12) << "a" << i+1 << "_x";
            Fout << setw(12) << "a" << i+1 << "_y";
            Fout << setw(12) << "v" << i+1 << "_x";
            Fout << setw(12) << "v" << i+1 << "_y";
            Fout << setw(12) << "x" << i+1 << "_x";
            Fout << setw(12) << "x" << i+1 << "_y";
    }
    Fout << endl;
    
    for (int i=0; i<N; i++){ 
            Fout << setw(15) << a[i][0];
            Fout << setw(15) << a[i][1];
            Fout << setw(15) << v[i][0];
            Fout << setw(15) << v[i][1];
            Fout << setw(15) << x[i][0];
            Fout << setw(15) << x[i][1];
    }
    Fout << endl;

    while(t <= 200000){
        calRho();
        if(t == 1){
            scaleRecal();
        
        }
        calPre();
        calVis();
        calGra();


        // calculate a
        for(int i=0; i< N; i++){
            F77NAME(dcopy)(2, F_p[i], 1, a[i], 1);
            F77NAME(daxpy)(2, 1.0, F_v[i], 1, a[i], 1);
            F77NAME(daxpy)(2, 1.0, F_g[i], 1, a[i], 1);
            F77NAME(dscal)(2, 1/rho_i[i], a[i], 1);
        }

        // time integration step
        if (t == 1){
            for(int i=0; i<N; i++){
                F77NAME(daxpy)(2, dt/2, a[i], 1, v[i], 1);
                F77NAME(daxpy)(2, dt, v[i], 1, x[i], 1);
            }
        }else{
            for(int i=0; i<N; i++){
                F77NAME(daxpy)(2, dt, a[i], 1, v[i], 1);
                F77NAME(daxpy)(2, dt, v[i], 1, x[i], 1);
            }
        }

        // check for BC
        for (int i=0; i<N; i++){
            if (x[i][0] < domain[0] + h){
                v[i][0] = -e*v[i][0];
                x[i][0] = domain[0] + h;
            }else if(x[i][0] > domain[1] - h){
                v[i][0] = -e*v[i][0];
                x[i][0] = domain[1] - h;
            }
            if (x[i][1] < domain[0] + h){
                v[i][1] = -e*v[i][1];
                x[i][1] = domain[0] + h;
            }else if(x[i][1] > domain[1] - h){
                v[i][1] = -e*v[i][1];
                x[i][1] = domain[1] - h;
            }
        }

//        // print result
//        for (int i=0; i<N; i++){
//            cout << "ax ay of particle: " << i+1 << endl;
//            cout << a[i][0] << " " << a[i][1] << endl;
//            cout << "vx vy xx xy of particle: " << i+1 << endl;
//            cout << v[i][0] << " " << v[i][1] << " ";
//            cout << x[i][0] << " " << x[i][1] << endl;
//        }
//        cout << "end of time step: " << t << endl;
//        cout << endl;
        if (t % 100 == 0){
            for (int i=0; i<N; i++){ 
                Fout << setw(15) << a[i][0];
                Fout << setw(15) << a[i][1];
                Fout << setw(15) << v[i][0];
                Fout << setw(15) << v[i][1];
                Fout << setw(15) << x[i][0];
                Fout << setw(15) << x[i][1];
            }
            Fout << endl;
        
            
        }

        t ++;

    }

}


int main(int argc, char const *argv[])
{
    SPH_serial test = SPH_serial(4);
    double loc[8] = {0.505, 0.5, 0.515, 0.5, 0.51, 0.45, 0.5, 0.45};

    // noise generation
    srand(time(0));
    for (int i=0; i<4 ; i++){
        double noise = (double) rand()/(RAND_MAX/2) - 1; // noise is in -1 to 1
        noise *= 0.01 / 10; // noise is scaled to -h/10 to h/10
        loc[i*2] += noise;
    }

    // for (int i=0; i<4; i++){
    //     // cout.precision(8); 
    //     cout << loc[i] << endl;
    // }
    test.inputLocation(loc);
    // test.calRho();
    // test.scaleRecal();
    // test.calPre();
    // test.calVis();
    // test.calGra();
    test.timeInte();
    return 0;
}
