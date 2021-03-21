#include "DataStructure.h"

// construct based on no. of particles and world size
DecomposeInfo::DecomposeInfo(int N, int size){

    this->N = N;
    this->size = size;
    N_proc = new int[size]();
    startloc = new int[size]();
    recvcounts_scalar = new int[size]();
    displs_scalar = new int[size]();
    recvcounts_vec = new int[size]();   
    displs_vec = new int[size]();      

    startloc[0] = 0;
    for(int i=0; i < size; i++){
        N_proc[i] = N / size; // each process hold at least N/size particles
        if (i < N % size){
            N_proc[i] += 1; // if there is remainder, evenly distribute remainder
        }
        if (i < size - 1){
            startloc[i+1] = startloc[i]+N_proc[i]; // calculate starting position
        }
        // set receive counts and the starting postion to hold communications
        recvcounts_scalar[i] = N_proc[i];
        displs_scalar[i] = startloc[i];
        recvcounts_vec[i] = N_proc[i] * 2;
        displs_vec[i] = startloc[i] * 2;
    }
}

DecomposeInfo::~DecomposeInfo(){
    delete[] N_proc;               
    delete[] startloc;                 
    delete[] recvcounts_scalar;    
    delete[] displs_scalar;        
    delete[] recvcounts_vec;       
    delete[] displs_vec;           
}
