#include "SPH_parallel.h"
#include <cmath>
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

    void F77NAME(dgemv)(const char & TRANS,
                        const int & M, const int & N,
                        const double & alpha, 
                        const double * A,
                        const int & LDA,
                        const double * x,
                        const int & incx,
                        const double & beta,
                        double * y,
                        const int & incy);

}

double * FindPair_brute(double * x, int N){
    double * r = new double [2 * N * N];
    for(int i = 0; i < N; i++){ // col
        for (int j = 0; j < N; j++){ // row
            F77NAME(dcopy)(2, x+ 2*j, 1, r + 2*(i*N +j) , 1);
            F77NAME(daxpy)(2 , -1.0, x+ 2*i, 1, r + 2*(i*N +j), 1);
        }
    }
    return r;
}


void decompose_particles(int N, int rank, int size, double * loc, 
                         int & N_proc, double ** substart){
    // if there are more process than particles, do not care about this for now
    if (size > N){
        throw runtime_error("too little particle");
    }
    // determine the no. of particles for each rank
    N_proc = N / size;
    // determine the starting point of input_loc for each rank
    *substart = loc + 2 * N_proc * rank;
    // if at last rank, put the remainder of particles in the last process
    if (rank == size - 1){
        N_proc += N % size;
    }
}

// constructor for SPH_parallel
SPH_parallel::SPH_parallel(int N, int rank, int N_proc){
    // input N and rank
    this -> N = N;
    this -> rank = rank;
    this -> N_proc = N_proc;
    // initialise class parameters, no double pointer is used
    x = new double[N_proc * 2]();          // x has N vectors of 2
    v = new double[N_proc * 2]();          // v has N vectors of 2
    r = new double[N_proc * N * 2]();      // r has N*N vectors of 2
    q = new double[N_proc * N]();          // q has N*N values
    p = new double[N_proc]();              // p has N values
    Fp = new double[N_proc * 2]();         // Fp has N vectors of 2
    Fv = new double[N_proc * 2]();         // Fv has N vectors of 2
    Fg = new double[N_proc * 2]();         // Fg has N vectors of 2
    rho = new double[N_proc]();            // rho has N values
    phi_d = new double[N_proc * N]();      // phi_d has N*N values
    phi_p = new double[N_proc * N * 2]();  // phi_p has N*N vectors of 2
    phi_v = new double[N_proc * N]();      // phi_v has N*N values
    a = new double[N_proc * 2]();          // a has N vectors of 2
}

// destructor for parallel
SPH_parallel::~SPH_parallel(){
    // Free all memory
    delete [] x;
    delete [] v;
    delete [] r;
    delete [] q;
    delete [] p;
    delete [] Fp;
    delete [] Fv;
    delete [] Fg;
    delete [] rho;
    delete [] phi_d;
    delete [] phi_p;
    delete [] phi_v;
    delete [] a;
    
}


// input location
void SPH_parallel::inputLocation(double * loc){
    // read locations from array loc(col major)
    for(int i = 0; i < N_proc; i++){      
        x[2*i] = loc[2*i];
        x[2*i + 1] = loc[2*i + 1]; 
        cout << x[2*i] << " " << rank << endl;
        cout << x[2*i + 1] << " " << rank << endl;      
    }
}

// calculate rho
void SPH_parallel::calRho(double * loc){
    double * ptr = phi_d + rank * N_proc;
    for(int i=0; i< N_proc; i++){
        *ptr = 4 / M_PI / h / h;
        ptr +=  N_proc + 1;
    }
    // printMatrix(phi_d, N_proc, N);
    
    double * temp = new double[N]; 
    fill(temp, temp + N, 1.0); // temp = [1,1,...1]
    F77NAME(dgemv)('N', N_proc, N, m, phi_d, N_proc, temp, 1, 0.0, rho, 1); // m * phi_d * temp
    delete [] temp;

    // for(int i=0; i< N_proc; i++){
    //    cout << rho[i] << " " << rank <<endl;
    // }

}
