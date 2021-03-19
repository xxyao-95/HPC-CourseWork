#include "SPH_parallel.h"
#include "mpi.h"
using namespace std;

int main(int argc, char *argv[])
{
    // initialise MPI
    MPI_Init(&argc, &argv);
    // determine rank, size
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // create N and loc to receive input parameters from root process
    int  N;         // no. of particles
    double * loc;   // coordinate of particles
    // if rank = 0 (root rank), parse argument
    if (rank == 0){
        cout << "root will parse argument" << endl;
        cout << "root will read in loc input" << endl;
        N = 7;
        loc = new double [N * 2];
        double testnum = 0.7;
        for (int i = 0; i < N; i++){
            loc[2*i] = testnum;
            loc[2*i + 1] = testnum;
            testnum -= 0.1;
        }
    }
    // broadcast no. of particles
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0){
        loc = new double[N*2];
    }
    // broadcast coordinate of particles
    MPI_Bcast(loc, N*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // decompose the particles into different SPH process
    int N_proc;         // particle per porcess
    double * substart;  // where to start read from loc
    decompose_particles(N, rank, size, loc, N_proc, &substart);
    // instanciate SPH_parallel class
    SPH_parallel test1 = SPH_parallel(N, rank, N_proc);
    cout << "process: " << rank << " initialized " << N_proc << " points" << endl;
    // input location for all processes
    test1.inputLocation(substart);

    // now start the algorithm at the 1st time step
    // need to calculate rij to determine q which decide whether the particles are colliding
    // let root calculalte rij
    if (rank == 0){
        double * r;
        r = FindPair_brute(loc, N);
        cout << r[2] << endl;
        cout << r[3] << endl;
    }

    test1.calRho(loc);

    
        
  

    
    // Finailze MPI
    MPI_Finalize();
    return 0;
}
