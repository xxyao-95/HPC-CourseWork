/*
Data structure to hold world information
and information for receiving communications
*/
# pragma once
class DecomposeInfo{
public:
    int N;                      // total number of particles
    int size;                   // 
    int * N_proc;               // number of particles for each process
    int * startloc;             // starting location of each process to read particle locations    
    int * recvcounts_scalar;    // the receive counts for scalar quantities
    int * displs_scalar;        // the starting postion to hold communication for scalar quantities
    int * recvcounts_vec;       // the receive counts for vector quantities
    int * displs_vec;           // the starting position to hold communicaiton for vector quantities

    // constructor based on size and no. of particles
    DecomposeInfo(int N, int size);
    // destructor
    ~DecomposeInfo();

};

