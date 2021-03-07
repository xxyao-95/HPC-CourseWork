#include <iostream>

using namespace std;


class SPM




int main(int argc, char const *argv[])
{
    double dt = 1E-3;
    double h = 0.01;
    double g = 9.81;
    double x = 0.5;
    long double result;
    for(int i=0; i < 100000; i++){
        x = (long double)x - (g*dt*dt/2);
        if (x >= h){
            break;
        }
        cout << x << endl;
    }
    
    return 0;
}
