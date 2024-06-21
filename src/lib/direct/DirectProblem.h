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
#include "Gauss.h"
#include "ShapeFunction.h"
#include "MathFEM.h"

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

        // Pysical parameter
        double Re, rho, mu, nu;

        // Time parameter
        double dt;
        int timeMax;
        int pulsatileFlow;
        int pulseBeginItr;
        double T;

        // Darcy parameter
        double alpha, resistance;

        void runSimulation();
        void preprocess();
        void solveUSNS();
        void matrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, 
                                const int ic, const int tItr);
        void DarcyMatrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, 
                                const int ic, const int tItr);

    private:
        void prepareSerialMatrix();
        void visualizeDomain();
        void setVelocityValue(double (&vel)[3], double (&advel)[3], double (&dvdx)[3][3],
                              std::vector<double> &N, std::vector<std::vector<double>> &dNdx, 
                              const int ic, const int tItr);
        void updataValiables(const int tItr);

};



#endif

