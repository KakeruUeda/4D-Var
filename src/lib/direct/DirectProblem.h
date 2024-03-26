#ifndef DIRECT_H
#define DIRECT_H

#include <iostream>
#include <omp.h>
#include <mpi.h>
#include "Grid.h"

struct EstimatedVariable
{
    public:
        Array1D<double> u, v ,w;
};

class DirectProblem
{
    public:
        DirectProblem(Config &conf);
        ~DirectProblem(){}
     
        Array1D<double> u, v, w;
        double Re, pho, mu, nu;

    private:
        Grid grid;


};



#endif

