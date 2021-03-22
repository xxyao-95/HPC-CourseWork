/*
Case generation header file made by Xiyao
This file contains the declaration of functions 
for generating cases and noise
*/

// header guard
# pragma once
// includes
#include <iostream>
#include <vector>
using namespace std;

// generate the validation case
void generate_validation(const int & N, vector<double> & loc, const double & h);
// generate the dambreak case
void generate_dambreak(int & N, vector<double> & loc, const double &  distance);
// generate the block drop case
void generate_blockdrop(int & N, vector<double> & loc, const double &  distance);
// generate the droplet case
void generate_droplet(int & N, vector<double> & loc, const double &  distance);
// function to add noise
void add_noise(const int & N, vector<double> & loc);