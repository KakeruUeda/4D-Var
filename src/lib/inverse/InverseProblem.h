#ifndef INVERSEPROBLEM_H
#define INVERSEPROBLEM_H

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
#include "Gauss.h"
#include "ShapeFunction.h"
#include "MathFEM.h"

extern MyMPI mpi;

struct EstimatedVariable
{
    public:
        std::vector<double> u, v ,w;
};

class InverseProblem
{
    public:
        InverseProblem(Config &conf);
        ~InverseProblem(){}
        void runSimulation();

    private:
        void prepareMatrix();
};

#endif