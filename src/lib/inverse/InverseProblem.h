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
        grid(conf), planeDir(conf.planeDir), timeMax(conf.timeMax), 
        rho(conf.rho), mu(conf.mu), dt(conf.dt),
        alpha(conf.alpha), resistance(conf.resistance)
        {}
        
        Grid grid;
        PetscSolver petsc;
        
        int timeMax;
        double dt;
        double rho, mu, nu, Re;
        double alpha, resistance;
        
        std::vector<int> planeDir;

        void solveAdjointEquation(DirectProblem &main, std::string outputDir, 
                                  std::vector<std::vector<std::vector<double>>> &feedbackForce);
        void matrixAssemblyAdjointUSNS(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal, 
                                       std::vector<std::vector<std::vector<double>>> &feedbackForce,
                                       int &st, const int ic, const int t);
        void boundaryIntegral(DirectProblem &main, MatrixXd &Klocal,
                              const int ic, const int ib);
        void updateVariables(std::string output, const int dim, const int t);

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

        double aCF, bCF1, bCF2, gCF;
        int loopMax;
        int nControlNodesInCell;

        std::vector<int> planeDir;
        std::vector<std::vector<std::vector<double>>> feedbackForce;
        std::vector<std::vector<std::vector<double>>> gradWholeNode;
        std::vector<std::vector<std::vector<double>>> grad;
        
        
        void initialize(Config &conf);
        void runSimulation();

        void calcCostFunction();
        void GaussIntegralRegTerm1(std::vector<double> &N, std::vector<std::vector<double>> &dNdr,
                                   std::vector<std::vector<double>> &xCurrent, double &value, 
                                   const double weight, const int ic, const int t);
        void GaussIntegralRegTerm2(std::vector<double> &N, std::vector<std::vector<double>> &dNdr,
                                   std::vector<std::vector<double>> &xCurrent, double &value, 
                                   const double weight, const int ic, const int t);
        void calcFeedbackForce();
        void calcEdgeValue(std::vector<std::vector<std::vector<std::vector<double>>>> &vEX, const int t);
        void calcInterpolatedFeeback(std::vector<std::vector<double>> &xCurrent, double (&feedback)[3], 
                                     std::vector<std::vector<std::vector<std::vector<double>>>> &vEX, 
                                     double (&point)[3]);
        void feedbackGaussIntegral(std::vector<double> &N, double (&feedback)[3], 
                                   const double detJ, const double weight, const int ic, const int t);
        void calcFeedbackForce2();
        void feedbackGaussIntegral2(std::vector<double> &N, std::vector<std::vector<double>> &xCurrent,
                                    std::vector<std::vector<double>> &velCurrent, const double detJ,
                                    const double weight, const int voxelId, const int cellId, const int t);

        void calcOptimalCondition();
        void GaussIntegralOpttimalConditionTerm1(std::vector<double> &N, std::vector<std::vector<double>> &dNdr, 
                                                 std::vector<std::vector<double>> &xCurrent, 
                                                 double (&value)[4][3], const double weight, 
                                                 const int ic, const int t);
        void GaussIntegralOpttimalConditionTerm2(std::vector<double> &N, std::vector<std::vector<double>> &dNdr, 
                                                 std::vector<std::vector<double>> &xCurrent, 
                                                 double (&value)[4][3], const double weight, 
                                                 const int ic, const int t);
        void GaussIntegralOpttimalConditionTerm3(std::vector<double> &N, std::vector<std::vector<double>> &dNdr, 
                                                 std::vector<std::vector<double>> &xCurrent, 
                                                 double (&value)[4][3], const double weight, 
                                                 const int ic, const int t);
        double armijoCriteria(const double fk);
    
    private:

};

#endif