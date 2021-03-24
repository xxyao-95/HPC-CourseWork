/*
Efficient way of Finding neighbouring particles made by xiyao

This algorithm will split the domain into grids of dimention h x h
each particle is put into the grid then to check for neighbours, 
only look for particles in adjacent grids, it is possible 
that more than 1 particle occupy a grid, so need to take care of this with a vector
since we do not know the number of particles in one grid.

This file contains the function implementations
*/

// includes
#include "Findneighbour.h"

// function to generate grid map based on total no. of particles and their coordinates
void putIntoGrid(unordered_map <int, vector<point>> & res, double * x_global, 
                const int & N, const double & h, const int & n_grid){
    // pointer to x_global
    double * ptr = x_global;
    for(int i=0; i < N; i++){
        // calcualte the coordinate of the point in grid map
        int grid_x = (int) floor(ptr[0] / h);
        int grid_y = (int) floor(ptr[1] / h);
        int grid_coor = grid_x * n_grid + grid_y;
        // put the information of particle into p
        point p;
        p.x_loc[0] = ptr[0];
        p.x_loc[1] = ptr[1];
        p.loc = i;
        // put particle information into the map at the corresponding grid coordinate
        if (res.find(grid_coor)!=res.end()){
            // if the key exist, append to the vector
            res[grid_coor].push_back(p);
        }else{
            // if key not exist, 
            vector<point> data;
            data.push_back(p);
            res[grid_coor] = data;
        }
        // increment the pointer to check for next point
        ptr += 2;
    }
}
