/*
Parallel SPH source file made by xiyao
this file contains all the functions for the 
SPH_parallel class
*/

// includes
#include "SPH_parallel.h"

#define F77NAME(x) x##_

// BLAS routines
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
    // calculate number of grid
    // need to make sure the number of grid >= the largest y coordinate of grids 
    // the grid discretization include the upmost and righmost grid when
    // 1/h does not give an int 
    n_grid = (int) ceil((domain[1] -domain[0])/h);
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
    }
    // save the starting location into x_global
    for (int i = 0; i< 2*N; i++){
        x_global[i] = loc[i];
    }
}

// set dt h T
void SPH_parallel::setPara(double dt, double h, double T){
    this->dt = dt;
    this->h = h;
    this->T = T;
    // recalculate n_grid value changes
    n_grid = (int) ceil((domain[1] -domain[0])/h);
}

// calculate rho
void SPH_parallel::calRho(){
    fill(phi_d, phi_d + (N_proc * N), 0.0);
    // assgin phi_d values at postitions where rij = [0, 0];
    double * ptr = phi_d + info->startloc[rank] * N_proc; // find the starting location in global phi_d matrix
    for(int i=0; i< N_proc; i++){
        *ptr = 4 / M_PI / h / h; // phi_d = 4*pi/h^2 at i=j
        ptr +=  N_proc + 1; // increment to the next diagnal element
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
    // sum all the rhos
    double sum_rho = F77NAME(dasum)(N, rho_global, 1);
    // scale m
    m = N*rho_0 / sum_rho;
    // recalculate both global and local rho
    F77NAME(dscal)(N, m, rho_global, 1);
    F77NAME(dscal)(N_proc, m, rho, 1);
}


// Pressure Force Calculation
void SPH_parallel::calPre(){
    // make Fp and phi_p zero at begining
    fill(Fp, Fp + (N_proc * 2), 0.0);
    fill(phi_p, phi_p + (N*N_proc*2), 0.0);
    // calculate p
    fill(p_global, p_global + N, -k * rho_0); //-k * rho_0 -> p_global
    F77NAME(daxpy)(N, k, rho_global, 1, p_global, 1); // rho_global * k + p_global
    fill(p, p + N_proc, -k * rho_0); // -k * rho_0 -> p
    F77NAME(daxpy)(N_proc, k, rho, 1, p, 1); // rho * k + p


    // calculate phi_p at the required location
    for (int i: coordQ){
        // move pointers to the location
        double * ptr_r = r + (2 * i); // pointer to r at the location of non-zero q
        double * ptr_phi_p = phi_p + (2 * i); // pointer to phi_p at the location of non-zero q
        // calculate a scale factor
        double scal_fac = (-30/M_PI/pow(h,3)) * pow((1-q[i]),2)/q[i];
        // r -> phi_p
        ptr_phi_p[0] = ptr_r[0];
        ptr_phi_p[1] = ptr_r[1];
        // phi_p * scal_fac
        ptr_phi_p[0] *= scal_fac;
        ptr_phi_p[1] *= scal_fac;
    }
    // calculate F_p
    for (int i:coordQ){
        int col = i / N_proc; // get col no. of the non-zero q, the location of p and rho that need to be included in the formula
        int row = i % N_proc; // get row no. of the non-zero q, the particle that we want to calculate force
        double * ptr_phi_p = phi_p + (2 * i); // pointer that points to the phi_p calculated
        double * ptr_Fp = Fp + (2 * row); // pointer that points to the particle
        double scal_fac =  -(m/rho_global[col]) * (p[row] + p_global[col])/2; // scale factor
        // phi_p * scal_fac
        ptr_phi_p[0] *= scal_fac;
        ptr_phi_p[1] *= scal_fac;
        // Fp += phi_p*calfac
        ptr_Fp[0] += ptr_phi_p[0];
        ptr_Fp[1] += ptr_phi_p[1];
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
    double vij[2] = {0,0}; // a vector for vij
    for (int i: coordQ){
        int col = i / N_proc; // get col no. of the non-zero q, the location of p and rho that need to be included in the formula
        int row = i % N_proc; // get row no. of the non-zero q, the particle that we want to calculate force
        double scale_fac = -miu * (m/rho_global[col]) * phi_v[i]; // calculate the scale factor -miu *m/rho_i
        double * ptr_v = v +(2 * row); // pointer that points to vi
        double * ptr_v_global = v_global + (2 * col); // pointer that points to vj
        // vij = vi - vj
        vij[0] = ptr_v[0] - ptr_v_global[0];
        vij[1] = ptr_v[1] - ptr_v_global[1];
        // vij * scalfactor
        vij[0] *= scale_fac;
        vij[1] *= scale_fac;
        double * ptr_Fv = Fv + (2 * row); // point to the location of F_v that we want to calculate
        // Fv += vij * scal_fac
        ptr_Fv[0] += vij[0];
        ptr_Fv[1] += vij[1];
    }
    
}

// calculated gravity force
void SPH_parallel::calGra(){
    double * ptr = Fg + 1;
    F77NAME(dcopy)(N_proc, rho, 1, ptr, 2);
    F77NAME(dscal)(N_proc * 2, -g, Fg, 1);

}

// Perform time integration
void SPH_parallel::timeInte(){
    int t = 1;
    ofstream Fout_energy;
    // let last process write the output file
    if (rank == info->size -1){
        Fout_energy.open("enery.txt");
        Fout_energy << setw(15) << "Time Step";
        Fout_energy << setw(15) << "Kinetic Energy";
        Fout_energy << setw(15) << "Potential Energy";
        Fout_energy << setw(15) << "Total Energy";
        Fout_energy << endl;
    }


    int timesteps = (int) T/dt;
    while(t <= timesteps){

        calRijQ();     
        calRho();
        if(t == 1){
            scaleRecal();
            // write energy at initial time step
            // let last process do the work
            if (rank == info->size -1){
                calEnergy();
                Fout_energy << setw(15) << 0;
                Fout_energy << setw(15) << Ek;
                Fout_energy << setw(15) << Ep;
                Fout_energy << setw(15) << E_total;
                Fout_energy << endl;
            }
        }

        // calculate gravity first
        calGra();
        // copy gravity to a
        F77NAME(dcopy)(2 * N_proc, Fg, 1, a, 1);

        if (!coordQ.empty()){
            calPre(); // calculate pressure force only when needed
            calVis(); // calculate viscous force only when needed
            F77NAME(daxpy)(2 * N_proc, 1.0, Fv, 1, a, 1);
            F77NAME(daxpy)(2 * N_proc, 1.0, Fp, 1, a, 1);
        }

        double * ptr_a = a;
        for(int i = 0; i < N_proc; i++){
            // a = F/rho
            ptr_a[0] /= rho[i];
            ptr_a[1] /= rho[i];
            ptr_a += 2;
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
        sendRecvLoc();

        // collect energy and write to file at each time step
        if (rank == info->size - 1){
            calEnergy();
            Fout_energy << setw(15) << t;
            Fout_energy << setw(15) << Ek;
            Fout_energy << setw(15) << Ep;
            Fout_energy << setw(15) << E_total;
            Fout_energy << endl;  
        }
        // increment time step
        t++;
    }

    Fout_energy.close();
    // write the final location
    ofstream Fout_loc;
    if (rank == info->size - 1){
        Fout_loc.open("output.txt");
        Fout_loc << setw(15) << "Particle";
        Fout_loc << setw(15) << "x";
        Fout_loc << setw(15) << "y";
        Fout_loc << endl;
        for (int i=0; i<N; i++){ 
            Fout_loc << setw(15) << i+1;
            Fout_loc << setw(15) << x_global[2*i];
            Fout_loc << setw(15) << x_global[2*i + 1];
            Fout_loc << endl;
        }
    }

}

// Function to send and receive x and v
void SPH_parallel::sendRecvLoc(){

    MPI_Allgatherv(x, N_proc * 2, MPI_DOUBLE, x_global, 
                   info->recvcounts_vec, info->displs_vec, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(v, N_proc * 2, MPI_DOUBLE, v_global, 
                   info->recvcounts_vec, info->displs_vec, MPI_DOUBLE, MPI_COMM_WORLD);
    
}

// Function to calculate rij and q this function implements a proximity seach algorithm by 
// first discretizing the domian into grids and put all particles inthe grids
// each particle only need to check for 9 grids and this appraoch is roughly O(n)

void SPH_parallel::calRijQ(){
    fill(r, r + (N*N_proc*2), 0.0);
    fill(q, q + (N*N_proc), 0.0);
    // clear the vector for positon of non-zero q at every time step
    // clear gridMap at every time step
    coordQ.clear();
    gridMap.clear();
    // int n_grid = (int) ceil(1.0/h);
    putIntoGrid(gridMap, x_global, N, h, n_grid);
    // pointer to x
    double * ptr_x = x;
    for (int j=0; j < N_proc; j++){
        // for each point, get the grid
        int grid_x = (int) floor(ptr_x[0] / h);
        int grid_y = (int) floor(ptr_x[1] / h);
        int grid_coor = grid_x * n_grid + grid_y;
        // calculate the global index of j
        int j_global = info->startloc[rank] + j;
        // the array for checking the neighbours, there are total 9 grid boxes to check 
        int check[9] = {grid_coor, grid_coor-1, grid_coor+1, 
                        grid_coor-n_grid, grid_coor+n_grid, grid_coor-1-n_grid,
                        grid_coor-1+n_grid, grid_coor+1-n_grid, grid_coor+1+n_grid};
        // iterate through the grid boxes
        for(int c: check){
            // if there is a match
            if (gridMap.find(c) != gridMap.end()){
                // iterate through the points in that box
                for(point p : gridMap[c]){
                    int i = p.loc; // get the index of particle in global
                    if (i != j_global){
                        // calculate rij
                        double rijx = ptr_x[0] - p.x_loc[0];
                        double rijy = ptr_x[1] - p.x_loc[1];
                        r[2*(i*N_proc + j)] = rijx;
                        r[2*(i*N_proc + j)+1] = rijy;
                        // calculate q
                        q[i*N_proc + j] = sqrt(rijx*rijx + rijy*rijy) / h;
                        if (q[i*N_proc + j]<1 && q[i*N_proc + j]>0){
                            coordQ.push_back(i*N_proc + j);
                        }
                    }
                    
                }
            }
        }
        // increment pointer to x
        ptr_x += 2;
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
        // v_sqr = v*v
        double v_sqr = ptr_v[0]*ptr_v[0] + ptr_v[1]*ptr_v[1];
        Ek += 0.5 * m * v_sqr;
        Ep += m * g * ptr_x[0];
        // increment the pointers
        ptr_v += 2;
        ptr_x += 2;
    }
    E_total = Ek + Ep;
}