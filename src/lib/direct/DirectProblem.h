#ifndef DIRECT_H
#define DIRECT_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <set>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <mpi.h>
#include <omp.h>
#include <algorithm>
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
        DirectProblem(Config &conf);
        ~DirectProblem(){}

        int dim, nOMP;
        std::string outputDir;

        Grid grid;
        PetscSolver petsc;
     
        Array1D<double> u, v, w;
        double Re, pho, mu, nu;

        void runSimulation();
        void preprocess();

    private:
        void prepareSerialMatrix();
        void visualizeDomain();

};



#endif

