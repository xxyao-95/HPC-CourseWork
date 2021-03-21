#include <iostream>
#include <vector>
# pragma once
using namespace std;

void generate_droplet(int & N, vector<double> & loc, const double &  distance);

void generate_blockdrop(int & N, vector<double> & loc, const double &  distance);

void generate_dambreak(int & N, vector<double> & loc, const double &  distance);

void add_noise(const int & N, vector<double> & loc);