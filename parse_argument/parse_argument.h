/*
This document contains function declairations for parsing
command line arguments
made by Xiyao Liu

*/

// header guard
#pragma once

// includes
#include <iostream>
#include <vector>
#include <boost/program_options.hpp>

// name space
using namespace std;
namespace po = boost::program_options;

// function to parse in argument
void parse_argument(int * argc, char** argv[], po::options_description & opts,
                    po::variables_map & vm);

// function to check argument
void check_argument(po::variables_map & vm, int & case_name);

// funtion to generate case based on input argument
void generate_case(int case_name, int size, int & N, vector<double> & locvec, 
                   const double & h, const int & rank);