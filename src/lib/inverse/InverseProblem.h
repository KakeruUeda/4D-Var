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
#include "DirectProblem.h"

extern MyMPI mpi;

struct EstimatedVariable
{
    public:
        std::vector<double> u, v ,w;
};

struct CostFunction
{
    double term1,term2,term3,total;
    std::vector<double> history;
    void sum(){
        total = term1 + term2 + term3;
    }
};

class Adjoint
{
    public:
        Adjoint(Config &conf):
        grid(conf){}

        Grid grid;
        PetscSolver petsc;
};

class InverseProblem
{
    public:
        InverseProblem(Config &conf);
        ~InverseProblem(){}
        
        int dim, nOMP;
        std::string outputDir;

        Application app;
        DataGrid data;

        DirectProblem main;
        Adjoint adjoint;
        CostFunction costFunction;

        int timeMax;
        double dt;
        double rho, mu;
        double alpha, resistance;
        double aCF, bCF1, bCF2, gCF;
        int loopMax;
        int nControlNodesInCell;
        
        void initialize(Config &conf);
        void runSimulation();

        void calcCostFunction();
        void GaussIntegralRegTerm1(std::vector<double> &N, std::vector<std::vector<double>> &dNdr,
                                   std::vector<std::vector<double>> &xCurrent, double &value, 
                                   const double weight, const int ic, const int t);
        void GaussIntegralRegTerm2(std::vector<double> &N, std::vector<std::vector<double>> &dNdr,
                                   std::vector<std::vector<double>> &xCurrent, double &value, 
                                   const double weight, const int ic, const int t);
    private:

};

#endif