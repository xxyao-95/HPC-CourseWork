// #pragma once



#include "SPH_parallel.h"
#include <cmath>
#include "mpi.h"
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
    x = new double[N_proc * 2]();          // x has N_proc vectors of 2
    v = new double[N_proc * 2]();          // v has N_proc vectors of 2
    r = new double[N_proc * N * 2]();      // r has N_proc*N vectors of 2
    q = new double[N_proc * N]();          // q has N_proc*N values
    p = new double[N_proc]();              // p has N_proc values
    Fp = new double[N_proc * 2]();         // Fp has N_proc vectors of 2
    Fv = new double[N_proc * 2]();         // Fv has N_proc vectors of 2
    Fg = new double[N_proc * 2]();         // Fg has N_proc vectors of 2
    rho = new double[N_proc]();            // rho has N_proc values
    phi_d = new double[N_proc * N]();      // phi_d has N_proc*N values
    phi_p = new double[N_proc * N * 2]();  // phi_p has N_proc*N vectors of 2
    phi_v = new double[N_proc * N]();      // phi_v has N_proc*N values
    a = new double[N_proc * 2]();          // a has N vectors of 2
    rho_global = new double[N];            // rho of all N particles
    p_global = new double[N];              // P of all N particles
    x_global = new double[N * 2];              // x of all N particles
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
void SPH_parallel::calRho(){
    // assgin phi_d values at postitions where rij = [0, 0];
    double * ptr = phi_d + rank * N_proc * 2;
    for(int i=0; i< N_proc; i++){
        *ptr = 4 / M_PI / h / h;
        ptr +=  N_proc + 1;
    }
    // calculate phi_d for places where q < 1
    if (!coordQ.empty()){
        for (int i: coordQ){
            phi_d[i] = pow((1-q[i]*q[i]), 3) * 4 / M_PI /h /h;
        }
    }
    // printMatrix(phi_d, N_proc, N);
    
    double * temp = new double[N]; 
    fill(temp, temp + N, 1.0); // temp = [1,1,...1]
    F77NAME(dgemv)('N', N_proc, N, m, phi_d, N_proc, temp, 1, 0.0, rho, 1); // m * phi_d * temp
    delete [] temp;

    // gather all rhos into rho global
    int recvcounts[3] = {2,2,3};
    int displs[3] = {0,2,4};
    MPI_Allgatherv(rho, N_proc, MPI_DOUBLE, rho_global, 
                    recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

}

// scale m and recalculate rho
void SPH_parallel::scaleRecal(){

    double sum_rho = F77NAME(dasum)(N, rho_global, 1);
    m = N*rho_0 / sum_rho;

    F77NAME(dscal)(N, m, rho_global, 1);
    F77NAME(dscal)(N_proc, m, rho, 1);

    // for(int i = 0; i < N; i++){
    //     cout << rho_global[i] << " " <<rank << endl;
    // }
}


// Pressure Force Calculation
void SPH_parallel::calPre(){
    // calculate p
    fill(p_global, p_global + N, -k * rho_0);
    F77NAME(daxpy)(N, k, rho_global, 1, p_global, 1);

    // calculate phi_p
    // double * ptr = phi_p;
    // F77NAME(dcopy)(2, r[i*N +j], 1, ptr, 1); // phi_p[i*N + j] = r[i*N +j]

    // calculate F_p
    // F77NAME(dscal)(2, scale_fac, phi_p[j*N + i], 1);
    // F77NAME(daxpy)(2, 1.0, phi_p[j*N + i], 1, F_p[i], 1);
    
}

// calculate Viscous Force
void SPH_parallel::calVis(){
    // calculate vij
    int recvcounts[3] = {4, 4, 6};
    int displs[3] = {0, 4, 8};
    double * v_global = new double[N * 2];
    MPI_Allgatherv(v, N_proc * 2, MPI_DOUBLE, v_global, 
                    recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);


    // // calculate phi_v
    // for(int i = 0; i < N; i ++){
    //     for(int j = 0; j < N; j++ ){
    //         if (q[i*N + j] < 1 && i !=j ){
    //             phi_v[i*N + j] = (1 - q[i*N + j]) * 40/M_PI/pow(h, 4);
    //         }else{
    //             phi_v[i*N + j] = 0;
    //         }
    //     }
    // }


    // // calculate F_v
    // double scale_fac = -miu * (m/rho_i[j]) * phi_v[i*N + j];
    // F77NAME(dscal)(2, scale_fac, vij[j*N + i], 1);
    // F77NAME(daxpy)(2, 1.0, vij[j*N + i], 1, F_v[i], 1);
}

// calculated gravity force
void SPH_parallel::calGra(){
    double * ptr = Fg + 1;
    F77NAME(dcopy)(N_proc, rho, 1, ptr, 2);
    F77NAME(dscal)(N_proc * 2, g, Fg, 1);

}

void SPH_parallel::timeInte(){
    int t = 1;
    while(t <= 20){
        if(t > 1){
            calRijQ();
        }
        calRho();
        if(t == 1){
            scaleRecal();
        }
        if (!coordQ.empty()){
            calPre();
            calVis();
        }
        calGra();
        F77NAME(daxpy)(2 * N_proc, 1.0, Fp, 1, a, 1);
        F77NAME(daxpy)(2 * N_proc, 1.0, Fv, 1, a, 1);
        F77NAME(daxpy)(2 * N_proc, 1.0, Fg, 1, a, 1);
        double * ptr = a;
        for(int i = 0; i < N_proc; i++){
             F77NAME(dscal)(2, 1/rho[i], ptr, 1);
             ptr += 2;
        }

        // for(int i = 0; i < N_proc*2; i++){
        //     cout << a[i] << " " << rank << endl;
        // }

        // time integration step
        if (t == 1){
            F77NAME(daxpy)(N_proc * 2, dt/2, a, 1, v, 1);
            F77NAME(daxpy)(N_proc * 2, dt, v, 1, x, 1);

        }else{
            F77NAME(daxpy)(N_proc * 2, dt, a, 1, v, 1);
            F77NAME(daxpy)(N_proc * 2, dt, v, 1, x, 1);
        }

        // check for BC
        for (int i=0; i<N_proc; i++){
            if (x[i*2] < domain[0] + h){
                v[i*2] = -e*v[i*2];
                x[i*2] = domain[0] + h;
            }else if(x[i*2] > domain[1] - h){
                v[i*2] = -e*v[i*2];
                x[i*2] = domain[1] - h;
            }
            if (x[i*2 + 1] < domain[0] + h){
                v[i*2 + 1] = -e*v[i*2 + 1];
                x[i*2 + 1] = domain[0] + h;
            }else if(x[i*2 + 1] > domain[1] - h){
                v[i*2 + 1] = -e*v[i*2 + 1];
                x[i*2 + 1] = domain[1] - h;
            }
        }
        cout << "finished " << t <<" rank " << rank <<endl; 
        t++;
        // send and gather all data of locations
        MPI_Barrier(MPI_COMM_WORLD);
        sendRecvLoc();
    }

}

void SPH_parallel::sendRecvLoc(){

    int recvcounts[3] = {4, 4, 6};
    int displs[3] = {0, 4, 8};
    MPI_Allgatherv(x, N_proc * 2, MPI_DOUBLE, x_global, 
                    recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    
}


void SPH_parallel::calRijQ(){
    double * rij = new double[2 * N_proc * N];
    double * ptr = rij;
    double * xgptr = x_global;
    for(int i =0; i < N; i++){
        F77NAME(dcopy)(N_proc * 2, x, 1, ptr, 1);
        for (int j = 0; j < N_proc; j++){
            F77NAME(daxpy)(2, -1.0, xgptr,1,ptr,1);
            q[i*N_proc + j] = F77NAME(dnrm2)(2, ptr, 1) / h;
            // use a vector to save of coordinate of q
            if (q[i*N_proc + j] < 1){
                coordQ.push_back(i*N_proc + j);
            }
            ptr += 2;
        }
        xgptr += 2;      
    }


}