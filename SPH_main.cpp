#include "SPH_parallel.h"
#include "cases/generate_cases.h"
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[]){
    // initialise MPI
    MPI_Init(&argc, &argv);
    // determine rank, size
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // create N and loc to receive input parameters from root process
    int  N;         // no. of particles
    double * loc;   // coordinate of particles
    vector<double> locvec; // vector to store input
    // if rank = 0 (root rank), parse argument
    if (rank == 0){
        cout << "root will parse argument" << endl;
        cout << "root will read in loc input" << endl;
        // generate the case base on input argument
        // generate_validation(N, locvec, 0.01);
        generate_dambreak(N,locvec,0.02);
        add_noise(N, locvec);
        if (size > N){
            throw runtime_error("More process than particle");
        }

    }
    MPI_Barrier(MPI_COMM_WORLD);
    // broadcast no. of particles
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // make array for holding input locations
    loc = new double[2*N];
    if (rank == 0){
        // input the locations
        for(int i=0; i<2*N; i++){
            loc[i] = locvec[i];
        }
    }
    // broadcast coordinate of particles
    MPI_Bcast(loc, N*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // instanciate SPH_parallel class for each process
    SPH_parallel SPH = SPH_parallel(N, size, rank);
    // input location for all processes
    SPH.inputLocation(loc);
    // time integration
    SPH.timeInte();
    // Finailze MPI
    MPI_Finalize();
    cout << "finished" << endl;
    return 0;
}
