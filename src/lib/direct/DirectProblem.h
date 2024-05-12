#ifndef DIRECT_H
#define DIRECT_H

#include <iostream>
#include <omp.h>
#include <mpi.h>
#include "Grid.h"
#include "Boundary.h"
#include "PetscSolver.h"
#include "Config.h"

extern MyMPI mpi;

struct EstimatedVariable
{
    public:
        Array1D<double> u, v ,w;
};

class DirectProblem
{
    public:
        DirectProblem(Config conf);
        ~DirectProblem(){}

        Grid grid;
        PetscSolver petsc;
     
        Array1D<double> u, v, w;
        double Re, pho, mu, nu;

        void runSimulation();

    private:
        void prepareSerialMatrix();
        void prepareParallelMatrix();


};



#endif

