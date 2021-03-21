#include "generate_cases.h"
#include <cmath>
// #include <fstream>
// #include <iomanip>

void generate_droplet(int & N, vector<double> & loc, const double & distance){
    double center[2] = {0.5, 0.7};
    double top = center[1] + 0.1;
    double bottom = center[1] - 0.1;
    N = 0;
    double y = bottom;
    while(y <= top){
        double dy = abs(center[1] - y);
        double dx = sqrt(0.1*0.1 - dy*dy);
        double left = center[0] - dx;
        double right = center[0] + dx;
        double x = left;
        while (x <= right){
            loc.push_back(x);
            loc.push_back(y);
            N += 1;
            x += distance;
        }
        y += distance;
    }

}

// int main(int argc, char const *argv[])
// {
//     int N = 0;
//     vector<double> loc;
//     generate_droplet(N, loc, 0.01);
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
