
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

    // calculate rho for all particles in the process
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

    // scale m and rho
    F77NAME(dscal)(N, m, rho_global, 1);
    F77NAME(dscal)(N_proc, m, rho, 1);

    // for(int i = 0; i < N; i++){
    //     cout << rho_global[i] << " " <<rank << endl;
    // }
}


// Pressure Force Calculation
void SPH_parallel::calPre(){
    // make Fp and phi_p zero at begining
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
        F77NAME(dcopy)(2, ptr_r, 1, ptr_phi_p, 1); // r -> phi_p
        F77NAME(dscal)(2, scal_fac, ptr_phi_p, 1); // phi_p * scal_fac
    }
    // calculate F_p
    for (int i:coordQ){
        int col = i / N_proc; // get col no. of the non-zero q, the location of p and rho that need to be included in the formula
        int row = i % N_proc; // get row no. of the non-zero q, the particle that we want to calculate force
        double * ptr_phi_p = phi_p + (2 * i); // pointer that points to the phi_p calculated
        double * ptr_Fp = Fp + (2 * row); // pointer that points to the particle
        double scal_fac =  -(m/rho_global[col]) * (p[row] + p_global[col])/2; // scale factor
        F77NAME(dscal)(2, scal_fac, ptr_phi_p, 1); // phi_p * scal_fac
        F77NAME(daxpy)(2, 1.0, ptr_phi_p, 1, ptr_Fp, 1); // Fp += phi_p*calfac
    }
}

// calculate Viscous Force
void SPH_parallel::calVis(){
    // make Fv and phi_v zero at the begining
    fill(Fv, Fv + (N_proc * 2), 0.0);
    fill(phi_v, phi_v + (N_proc * N), 0.0);
    // calculate phi_v at the required non-zero q location
    for (int i: coordQ){
        phi_v[i] = (1-q[i]) * 40 / M_PI / pow(h, 4); // calculate phi_v
    }

    // calculate F_v
    for (int i: coordQ){
        int col = i / N_proc; // get col no. of the non-zero q, the location of p and rho that need to be included in the formula
        int row = i % N_proc; // get row no. of the non-zero q, the particle that we want to calculate force
        double scale_fac = -miu * (m/rho_global[col]) * phi_v[i]; // calculate the scale factor -miu *m/rho_i
        double *vij = new double[2](); // create a vector that hold vij for that particular coodinate
        double * ptr_v = v +(2 * row); // pointer that points to vi
        F77NAME(dcopy)(2, ptr_v, 1, vij, 1); // vi -> vij
        double * ptr_v_global = v_global + (2 * col); // pointer that points to vj
        F77NAME(daxpy) (2, -1.0, ptr_v_global, 1, vij, 1); // vij - vj
        F77NAME(dscal)(2, scale_fac, vij, 1); // vij * scalfactor
        double * ptr_Fv = Fv + (2 * row); // point to the location of F_v that we want to calculate
        F77NAME(daxpy)(2, 1.0, vij, 1,  ptr_Fv, 1); // Fv += vij * scal_pac
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
    ofstream Fout_loc;
    ofstream Fout_energy;
    // if (rank == 0){
    //     Fout_loc.open("output_droplet.txt");
    //     for (int i=0; i<N; i++){ 
    //         Fout_loc << setw(15) << "x_x";
    //         Fout_loc << setw(15) << "x_y";
    //     }
    //     Fout_loc << endl;

    //     for (int i=0; i<N; i++){ 
    //         Fout_loc << setw(15) << x_global[2*i];
    //         Fout_loc << setw(15) << x_global[2*i + 1];
    //     }
    //     Fout_loc << endl;
    // }

    if (rank = info->size -1){
        Fout_energy.open("enery.txt");
        Fout_energy << setw(15) << "Time Step";
        Fout_energy << setw(15) << "Kinetic Energy";
        Fout_energy << setw(15) << "Potential Energy";
        Fout_energy << setw(15) << "Total Energy";
        Fout_energy << endl;
    }

    while(t <= 200000){

        calRijQ();     
        calRho();
        if(t == 1){
            scaleRecal();
            if (rank = info->size -1){
                calEnergy();
                Fout_energy << setw(15) << 0;
                Fout_energy << setw(15) << Ek;
                Fout_energy << setw(15) << Ep;
                Fout_energy << setw(15) << E_total;
                Fout_energy << endl;
            }
        }

        if (!coordQ.empty()){
            calPre();
            calVis();
        }
        calGra();

        F77NAME(dcopy)(2 * N_proc, Fg, 1, a, 1);
        if (!coordQ.empty()){
            F77NAME(daxpy)(2 * N_proc, 1.0, Fv, 1, a, 1);
            F77NAME(daxpy)(2 * N_proc, 1.0, Fp, 1, a, 1);
        }
        double * ptr = a;
        for(int i = 0; i < N_proc; i++){
             F77NAME(dscal)(2, 1/rho[i], ptr, 1);
             ptr += 2;
        }

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

        // send and gather all data of locations
        MPI_Barrier(MPI_COMM_WORLD);
        sendRecvLoc();

        // write to file at root process
        // if (rank == 0 and t % 100 == 0){
        //     for (int i=0; i<N; i++){ 
        //         Fout_loc << setw(15) << x_global[2*i];
        //         Fout_loc << setw(15) << x_global[2*i + 1];
        //     }
        //     Fout_loc << endl;
        // }
        // collect energy and write to file
        if (rank == info->size - 1){
            calEnergy();
            Fout_energy << setw(15) << t;
            Fout_energy << setw(15) << Ek;
            Fout_energy << setw(15) << Ep;
            Fout_energy << setw(15) << E_total;
            Fout_energy << endl;  
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
    // need to clear the vector for positon of non-zero q at every time step
    coordQ.clear();
    // use pointers to traverse r and x_global
    double * ptr = r;
    double * xgptr = x_global;
    for(int i =0; i < N; i++){ // col of r
        F77NAME(dcopy)(N_proc * 2, x, 1, ptr, 1); // x -> r
        for (int j = 0; j < N_proc; j++){ // row of r
            F77NAME(daxpy)(2, -1.0, xgptr,1,ptr,1); // r - x_global
            q[i*N_proc + j] = F77NAME(dnrm2)(2, ptr, 1) / h; // find norm
            // use a vector to save of coordinate of q
            if (q[i*N_proc + j] < 1 && q[i*N_proc + j]>0){ // record position on non-zero q
                // cout << "Hi " << rank << endl;
                coordQ.push_back(i*N_proc + j);
            }
            ptr += 2; // increment pointer by 2 to move to next point
        }
        xgptr += 2; // increment pointer by 2 to move to next point     
    }
}

// energy calculation
void SPH_parallel::calEnergy(){
    Ek =0;
    Ep =0;
    E_total = 0;
    double * ptr_v = v_global;
    double * ptr_x = x_global+1;
    for(int i=0; i<N; i++){
        double vnorm = F77NAME(dnrm2)(2,ptr_v,1);
        Ek += 0.5 * m * pow(vnorm,2);
        Ep += m * g * ptr_x[0];
        ptr_v += 2;
        ptr_x += 2;
    }
    E_total = Ek + Ep;
}