#include "SPH_parallel.h"
#include "cases/generate_cases.h"
#include <cstdlib>
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
        vector<double> locvec;
        generate_droplet(N,locvec,0.01);
        if (size > N){
            throw runtime_error("More process than particle");
        }

        loc = new double[2*N];
        for(int i=0; i<2*N; i++){
            loc[i] = locvec[i];
            // cout << locvec[i] << endl;
        }
        // noise generation
        // srand(time(0));
        // for (int i=0; i<4 ; i++){
        //     double noise = (double) rand()/(RAND_MAX/2) - 1; // noise is in -1 to 1
        //     noise *= 0.01 / 10; // noise is scaled to -h/10 to h/10
        //     loc[i*2] += noise;
        // }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // broadcast no. of particles
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0){
        loc = new double[N*2];
    }
    // broadcast coordinate of particles
    MPI_Bcast(loc, N*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // instanciate SPH_parallel class
    SPH_parallel SPH = SPH_parallel(N, size, rank);
    
    // input location for all processes
    SPH.inputLocation(loc);

    // SPH.calRho();
  
    // SPH.scaleRecal();
    // test1.calPre();
    // test1.calVis();
    // test1.calGra();
    SPH.timeInte();
    // Finailze MPI
    MPI_Finalize();
    cout << "finished" << endl;
    return 0;
}
