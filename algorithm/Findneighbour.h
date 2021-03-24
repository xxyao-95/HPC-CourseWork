/*
Function to find neighbouring particles made by xiyao

This algorithm will split the domain into grids of dimention h x h
each particle is put into the grid then to check for neighbours, 
only look for particles in adjacent grids, it is possible 
that more than 1 particle occupy a grid, so need to take care of this.

This file contains the funtion declairations
*/

// header guard
#pragma once
// includes
#include <vector>
#include <unordered_map>
#include <cmath>
#include "../datastructure/DataStructure.h"
// namespace
using namespace std;

// some simple structure to hold information of a particle
typedef struct particle_info{
    // coordinate of the point
    double x_loc[2];
    // index of the point in global frame
    int loc;
}point;

// function to generate grid map based on total no. of particles and their coordinates
void putIntoGrid(unordered_map <int, vector<point>> & res, double * x_global, const int & N, const double & h, const int & n_grid);

// the implementation of this algorithm in the calRijQ() function in SPH_parallel class
