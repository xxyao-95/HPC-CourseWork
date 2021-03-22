/*
Case generation source file made by Xiyao
This file contains the functions for generating 
cases and noise
*/

// includes
#include "generate_cases.h"
#include <cmath>
#include <fstream>
#include <iomanip>

// generate the validation cases
void generate_validation(const int & N, vector<double> & loc, const double & h){
    // location of particles of the validations cases as listed in the assignment
    double case1[2] = {0.5, 0.5};
    double case2[4] = {0.5, 0.5, 0.5, h};
    double case3[6] = {0.5, 0.5, 0.495, h, 0.505, h};
    double case4[8] = {0.505, 0.5, 0.515, 0.5, 0.51, 0.45, 0.5, 0.45};
    // switch on the no. of particles
    switch (N){
        case 1:
            for(double i: case1){loc.push_back(i);}
            break;
        case 2:
            for(double i: case2){loc.push_back(i);}
            break;
        case 3:
            for(double i: case3){loc.push_back(i);}
            break;
        case 4:
            for(double i: case4){loc.push_back(i);}
            break;        
    }
}

// generate dam break case
void generate_dambreak(int & N, vector<double> & loc, const double &  distance){
    // initial conditions
    N = 0;
    double bottom = distance;
    double top = 0.2;
    double left = distance;
    double right = 0.2;
    // use 2 while loop to locate and input the point
    double y = bottom;
    while (y <= top){
        double x = left;
        while(x <= right){
            loc.push_back(x);
            loc.push_back(y);
            N += 1;
            x += distance;
        }
        y += distance;
    }     
}

// generate block drop case
void generate_blockdrop(int & N, vector<double> & loc, const double &  distance){
    N = 0;
    double bottom = 0.3;
    double top = 0.6;
    double left = 0.1;
    double right = 0.3;
    // use 2 while loop to locate and input the point
    double y = bottom;
    while (y <= top){
        double x = left;
        while(x <= right){
            loc.push_back(x);
            loc.push_back(y);
            N += 1;
            x += distance;
        }
        y += distance;
    }     
}

// generate droplet case
void generate_droplet(int & N, vector<double> & loc, const double & distance){
    double center[2] = {0.5, 0.7};
    double top = center[1] + 0.1;
    double bottom = center[1] - 0.1;
    N = 0;
    // use 2 while loop to locate and input the point
    double y = bottom;
    while(y <= top){
        // calculate distance from center for x
        double dy = abs(center[1] - y);
        double dx = sqrt(0.1*0.1 - dy*dy);
        double left_bound = center[0] - dx;
        double right_bound = center[0] + dx;
        // fill right part
        double x = center[0];
        while (x <= right_bound){
            loc.push_back(x);
            loc.push_back(y);
            N += 1;
            x += distance;
        }
        // fill left part
        x = center[0]-distance; // avoid duplicates
        while (x >= left_bound){
            loc.push_back(x);
            loc.push_back(y);
            N += 1;
            x -= distance;
        }
        y += distance;
    }

}

// add some noise at the starting location
void add_noise(const int & N, vector<double> & loc){
    srand(time(0));
    // add a bit of noise for both x and y locations
    for (int i=0; i<N*2 ; i++){
        double noise = (double) rand()/(RAND_MAX/2) - 1; // noise is in -1 to 1
        noise *= 0.01 / 100; // noise is scaled to -h/100 to h/100
        loc[i] += noise;
    }
}


// // test generate case functions
// int main(int argc, char const *argv[])
// {
//     int N = 0;
//     vector<double> loc;
//     generate_droplet(N, loc, 0.02);
//     cout << N << endl;
//     ofstream Fout;
//     Fout.open("test_droplet.txt");
//     int count = 1;
//     for(double i : loc){
//         Fout << setw(15) << i;
//         if(count % 2 == 0){
//             Fout << endl;
//         }
//         count += 1;
//     }

//     return 0;
// }
