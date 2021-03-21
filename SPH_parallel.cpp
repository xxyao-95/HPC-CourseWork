// #pragma once



#include "SPH_parallel.h"
#include <cmath>
#include "mpi.h"
#include <fstream>
#include <iomanip>

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


// constructor for SPH_parallel
SPH_parallel::SPH_parallel(int N, int size, int rank){
    // input N and rank
    this -> N = N;
    this -> rank = rank;
    info = new DecomposeInfo(N, size);
    this -> N_proc = info->N_proc[rank];
    cout << "Process "<< rank << " is working on "<< N_proc << " points" <<endl;

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
    rho_global = new double[N]();            // rho of all N particles
    p_global = new double[N]();              // P of all N particles
    x_global = new double[N * 2]();          // x of all N particles
    v_global = new double[N * 2]();          // v of all N particles
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
// every process need to know the location of all particles at start
// every process only take care of time integration of the particle it need to check
void SPH_parallel::inputLocation(double * loc){
    // read locations from array loc(col major)
    double * ptr = loc + 2*info->startloc[rank];
    for(int i = 0; i < N_proc; i++){      
        x[2*i] = ptr[2*i];
        x[2*i + 1] = ptr[2*i + 1]; 
        // cout << x[2*i] << " " << rank << endl;
        // cout << x[2*i + 1] << " " << rank << endl;      
    }
    for (int i = 0; i< 2*N; i++){
        x_global[i] = loc[i];
        // cout << x_global[i] << " " << rank << endl;
    }
}

// calculate rho
void SPH_parallel::calRho(){
    fill(phi_d, phi_d + (N_proc * N), 0.0);
    // assgin phi_d values at postitions where rij = [0, 0];
    double * ptr = phi_d + info->startloc[rank];
    for(int i=0; i< N_proc; i++){
        *ptr = 4 / M_PI / h / h;
        ptr +=  N_proc + 1;
    }
    // calculate phi_d for places where q < 1
    if (!coordQ.empty()){
        for (int i: coordQ){
            phi_d[i] = pow((1.0-pow(q[i], 2.0)), 3.0) * 4 / (M_PI * pow(h, 2.0));
        }
    }

    double * temp = new double[N]; 
    fill(temp, temp + N, 1.0); // temp = [1,1,...1]
    F77NAME(dgemv)('N', N_proc, N, m, phi_d, N_proc, temp, 1, 0.0, rho, 1); // m * phi_d * temp
    delete [] temp;
    
    // gather all rhos into rho global
    MPI_Allgatherv(rho, N_proc, MPI_DOUBLE, rho_global, 
                    info->recvcounts_scalar, info->displs_scalar, MPI_DOUBLE, MPI_COMM_WORLD);

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
    fill(Fp, Fp + (N_proc * 2), 0.0);
    fill(phi_p, phi_p + (N*N_proc*2), 0.0);
    // calculate p
    fill(p_global, p_global + N, -k * rho_0);
    F77NAME(daxpy)(N, k, rho_global, 1, p_global, 1);
    fill(p, p + N_proc, -k * rho_0);
    F77NAME(daxpy)(N_proc, k, rho, 1, p, 1);


    // calculate phi_p at the required location
    for (int i: coordQ){
        double * ptr_r = r + (2 * i); // pointer to r at the location of non-zero q
        double * ptr_phi_p = phi_p + (2 * i); // pointer to phi_p at the location of non-zero q
        double scal_fac = (-30/M_PI/pow(h,3)) * pow((1-q[i]),2)/q[i];
        F77NAME(dcopy)(2, ptr_r, 1, ptr_phi_p, 1);
        F77NAME(dscal)(2, scal_fac, ptr_phi_p, 1);
    }
    // calculate F_p
    for (int i:coordQ){
        int col = i / N_proc; // get col no. of the non-zero q
        int row = i % N_proc; // get row no. of the non-zero q
        double * ptr_phi_p = phi_p + (2 * i); // pointer that point to the phi_p calculated
        double * ptr_Fp = Fp + (2*row); // pointer that point to the particle
        double scal_fac =  -(m/rho_global[col]) * (p[row] + p_global[col])/2;
        F77NAME(dscal)(2, scal_fac, ptr_phi_p, 1);
        F77NAME(daxpy)(2, 1.0, ptr_phi_p, 1, ptr_Fp, 1);
    }
}

// calculate Viscous Force
void SPH_parallel::calVis(){
    fill(Fv, Fv + (N_proc * 2), 0.0);
    fill(phi_v, phi_v + (N_proc * N), 0.0);
    // calculate phi_v at the required non-zero q location
    for (int i: coordQ){
        phi_v[i] = (1-q[i]) * 40 / M_PI / pow(h, 4);
    }

    // calculate F_v
    for (int i: coordQ){
        int col = i / N_proc; // get col no. of the non-zero q
        int row = i % N_proc; // get row no. of the non-zero q
        double scale_fac = -miu * (m/rho_global[col]) * phi_v[i]; // calculate the scale factor -miu *m/rho_i
        double *vij = new double[2];
        double * ptr_v = v +(2 * row);
        F77NAME(dcopy)(2, ptr_v, 1, vij, 1);
        double * ptr_v_global = v_global + (2 * col);
        F77NAME(daxpy) (2, -1.0, ptr_v_global, 1, vij, 1);
        F77NAME(dscal)(2, scale_fac, vij, 1);
        double * ptr_fv = Fv + (2 * row);
        F77NAME(daxpy)(2, 1.0, vij, 1,  ptr_fv, 1);
    }
    
}

// calculated gravity force
void SPH_parallel::calGra(){
    double * ptr = Fg + 1;
    F77NAME(dcopy)(N_proc, rho, 1, ptr, 2);
    F77NAME(dscal)(N_proc * 2, -g, Fg, 1);

}

void SPH_parallel::timeInte(){
    int t = 1;
    ofstream Fout;
    if (rank == 0){
        Fout.open("output_4_noise_para.txt");
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
            Fout << setw(15) << 0;
            Fout << setw(15) << 0;
            Fout << setw(15) << v_global[2*i];
            Fout << setw(15) << v_global[2*i + 1];
            Fout << setw(15) << x_global[2*i];
            Fout << setw(15) << x_global[2*i + 1];
        }
        Fout << endl;
    }
    
    while(t <= 200000){

        calRijQ();     
        calRho();
        if(t == 1){
            scaleRecal();
        }

        // if (rank == 0 && t == 3100){
        //     for(int i=0; i<N*N_proc; i++){
        //         cout << phi_d[i] << endl;
        //     }
        // }

        if (!coordQ.empty()){
            calPre();
            calVis();
        }
        calGra();

        // if (rank == 1){
        //     for (int i=0; i<N_proc; i++){
        //         cout << rho[i] << endl;
        //     }     
        // }

        F77NAME(daxpy)(2 * N_proc, 1.0, Fp, 1, a, 1);
        F77NAME(daxpy)(2 * N_proc, 1.0, Fv, 1, a, 1);
        F77NAME(daxpy)(2 * N_proc, 1.0, Fg, 1, a, 1);
        double * ptr = a;
        for(int i = 0; i < N_proc; i++){
             F77NAME(dscal)(2, 1/rho[i], ptr, 1);
             ptr += 2;
        }

        // if (rank == 2){
        //     for (int i=0; i<N_proc; i++){
        //         cout << a[2*i] << endl;
        //         cout << a[2*i + 1] << endl;
        //     }     
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

        // if (rank == 2){
        //     for (int i=0; i<N_proc; i++){
        //         cout << x[2*i] << endl;
        //         cout << x[2*i + 1] << endl;
        //     }
            
        // }
        // send and gather all data of locations
        MPI_Barrier(MPI_COMM_WORLD);
        sendRecvLoc();

        // write to file at root process
        if (rank == 0 and t % 100 == 0){
            for (int i=0; i<N; i++){ 
                Fout << setw(15) << 0;
                Fout << setw(15) << 0;
                Fout << setw(15) << v_global[2*i];
                Fout << setw(15) << v_global[2*i + 1];
                Fout << setw(15) << x_global[2*i];
                Fout << setw(15) << x_global[2*i + 1];
            }
            Fout << endl;
        }
        t++;
        MPI_Barrier(MPI_COMM_WORLD);
    }

}

void SPH_parallel::sendRecvLoc(){

    MPI_Allgatherv(x, N_proc * 2, MPI_DOUBLE, x_global, 
                   info->recvcounts_vec, info->displs_vec, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(v, N_proc * 2, MPI_DOUBLE, v_global, 
                   info->recvcounts_vec, info->displs_vec, MPI_DOUBLE, MPI_COMM_WORLD);
    
}


void SPH_parallel::calRijQ(){
    coordQ.clear();
    double * ptr = r;
    double * xgptr = x_global;
    for(int i =0; i < N; i++){
        F77NAME(dcopy)(N_proc * 2, x, 1, ptr, 1);
        for (int j = 0; j < N_proc; j++){
            F77NAME(daxpy)(2, -1.0, xgptr,1,ptr,1);
            q[i*N_proc + j] = F77NAME(dnrm2)(2, ptr, 1) / h;
            // use a vector to save of coordinate of q
            if (q[i*N_proc + j] < 1 && q[i*N_proc + j]>0){
                // cout << "Hi " << rank << endl;
                coordQ.push_back(i*N_proc + j);
            }
            ptr += 2;
        }
        xgptr += 2;      
    }

}