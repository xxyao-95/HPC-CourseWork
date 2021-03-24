/*
This document contains function for parsing
command line arguments
made by Xiyao Liu

*/

// includes
#include "parse_argument.h"
#include "../cases/generate_cases.h"

// function to parse in argument
void parse_argument(int * argc, char** argv[], po::options_description & opts,
                    po::variables_map & vm){
    // specify options
    opts.add_options()
    ("ic-dam-break", "Use dam-break initial condition.")
    ("ic-block-drop", "Use block-drop inital condition.")
    ("ic-droplet", "Use droplet initial condition.")
    ("ic-one-particle", "Use one particle validation case initial condition.")
    ("ic-two-particles", "Use two particles validation case initial condition.")
    ("ic-three-particles", "Use three particles validation case initial condition.")
    ("ic-four-particles", "Use four particles validation case initial condition.")
    ("dt", po::value<double>()->default_value(1E-4), "Time-step to use.")
    ("T", po::value<double>()->default_value(2.0), "Total integration time.")
    ("h", po::value<double>()->default_value(0.01), "Radius of influence of each particle.")
    ("help", "Print help message");
    // parse argument to vm
    po::store(po::parse_command_line(*argc, *argv, opts), vm);
    po::notify(vm);

}

// function to check argument
void check_argument(po::variables_map & vm, int & case_name){
    if (vm.count("ic-dam-break")) {
        case_name = 5;
    }else if(vm.count("ic-block-drop")){
        case_name = 6;
    }else if(vm.count("ic-droplet")){
        case_name = 7;
    }else if(vm.count("ic-one-particle")){
        case_name = 1;
    }else if(vm.count("ic-two-particles")){
        case_name = 2;
    }else if(vm.count("ic-three-particles")){
        case_name = 3;
    }else if(vm.count("ic-four-particles")){
        case_name = 4;
    }else{
        throw logic_error("Please provide valid inputs !!!!");
    }
}

// funtion to generate case based on input argument
void generate_case(int case_name, int size, int & N, vector<double> & locvec, const double & h, const int & rank){
    switch (case_name)
    {
    case 1:
        N = 1;
        generate_validation(N, locvec, h);
        if (rank == 0){
            cout << "Running one particle case..." << endl;
        }
        
        break;
    case 2:
        N = 2;
        generate_validation(N, locvec, h);
        if (rank == 0){
            cout << "Running two particle case..." << endl;
        }
        
        break;
    case 3:
        N = 3;
        generate_validation(N, locvec, h);
        if (rank == 0){
            cout << "Running three particle case..." << endl;
        }
        
        break;
    case 4:
        N = 4;
        generate_validation(N, locvec, h);
        if (rank == 0){
            cout << "Running four particle case..." << endl;
        }
        
        break;
    case 5:
        generate_dambreak(N, locvec, h);
        if (rank == 0){
            cout << "Running dam break case..." << endl;
        }
        
        break;
    case 6:
        generate_blockdrop(N, locvec, h);
        if (rank == 0){
            cout << "Running block drop case..." << endl;
        }
        break;
    case 7:
        generate_droplet(N, locvec, h);
        if (rank == 0){
            cout << "Running droplet case..." << endl;
        }
        break;
    }

    // throw an error when size > N
    if (size > N){
        throw "There are more process than input particles !!!!";
    }
}