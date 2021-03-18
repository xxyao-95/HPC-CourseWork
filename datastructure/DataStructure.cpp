/*
Some potential data structurs

*/



#include<iostream>
#include<cmath>


class VectorOfVector{
    private:
    double * x = nullptr;
    
    public:
    VectorOfVector(int no_points){
        // initialise x
        x = new double[2 * no_points];
    }
    void inputLocal(double * loc){
        // input initial location of x
        x = loc;
    }
    double * getVecAt(int point){
        return x + (point * 2 - 2);
    }

};




class MatrixOfVector{
    private:
    double * r = nullptr;
    int N;

    public:
    MatrixOfVector(int no_points){
        N = no_points;
        r = new double[no_points * (no_points - 1) + 2];
    }

    void update(VectorOfVector x);

    double * getVecAt(int row, int col){
        double * result;

        return result;
    }

};