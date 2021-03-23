#include "SPH_parallel.h"
#include "cases/generate_cases.h"
#include <cstdlib>
#include "parse_argument/parse_argument.h"
using namespace std;

int main(int argc, char *argv[]){
    // initialise MPI
    MPI_Init(&argc, &argv);
    // determine rank, size
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // parameters for input
    po::options_description opts( "HPC course work, SPH algorithm made by Xiyao Liu");
    po::variables_map vm;

    // parse argument
    try{
        parse_argument(&argc, &argv, opts, vm);
    }catch (boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::unknown_option> > & e){
        if(rank == 0){
            cout << e.what() << endl;
        }
        MPI_Finalize();
        return 0;
    }
    // parse_argument(&argc, &argv, opts, vm);
    // check for help
    if (vm.count("help")) {
        if (rank == 0){
            cout << "Performs SPH algorithm to model behaviour of  particles" << endl;
            cout << opts << endl;
        }
        MPI_Finalize();
        return 0;
    }
    // This parameter specify which case to run
    int case_name;
    // create N and loc to receive input parameters from root process
    int  N;         // no. of particles
    double * loc;   // coordinate of particles
    vector<double> locvec; // vector to store input
    // readin some parameters
    const double dt = vm["dt"].as<double>();
    const double T = vm["T"].as<double>();
    const double h = vm["h"].as<double>();

    try{
        check_argument(vm, case_name);
        generate_case(case_name, size, N, locvec, h, rank);
        if (rank == 0){
            add_noise(N, locvec);
        }
    }catch(logic_error & e){
        if(rank == 0){
            cout << e.what() << endl;
        }
        MPI_Finalize();
        return 0;
    }catch(const char* msg){
        if (rank == 0){
            cout << msg << endl;
            cout << "Please ensure no. of process < no. of particles" << endl;
        }
        MPI_Finalize();
        return 0;
    }
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
    // set dt T h
    SPH.setPara(dt, h, T);
    // time integration
    SPH.timeInte();
    // Finailze MPI
    cout << "Process " << rank <<" finished" << endl;
    MPI_Finalize();
    return 0;
}
