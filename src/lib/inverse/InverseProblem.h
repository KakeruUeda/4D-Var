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
#include "Tool.h"
#include "Spline.h"
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


class Adjoint : public MathFEM
{
    public:
        Adjoint(Config &conf):
        grid(conf), dim(conf.dim), planeDir(conf.planeDir), 
        timeMax(conf.timeMax), rho(conf.rho), mu(conf.mu), dt(conf.dt),
        alpha(conf.alpha), resistance(conf.resistance)
        {}
        
        Grid grid;
        PetscSolver petsc;
        
        int dim;
        int timeMax;
        double dt;
        double rho, mu, nu, Re;
        double alpha, resistance;

        std::vector<double> vk, vk1, vk2;
        std::vector<double> advk1, advk2;
        std::vector<double> wk1, wk2;
        std::vector<std::vector<double>> dvkdx, dvk1dx, dvk2dx;
        std::vector<std::vector<double>> dwk1dx, dwk2dx;
        
        std::vector<int> planeDir;

        void initializeFEM();
        void solveAdjoint(DirectProblem &main, std::string outputDir, 
                          std::vector<std::vector<std::vector<double>>> &feedbackForceT,  
                          const int nData, const int loop);
        void solveAdjointDO(DirectProblem &main, std::string outputDir,
                             std::vector<std::vector<std::vector<double>>> &feedbackForceT,
                             const int nData, const int loop);
        void setValue(DirectProblem &main, std::vector<double> &N, 
                      std::vector<std::vector<double>> &dNdx, const int ic, const int t);
        void matrixAssemblyAdjointUSNS(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal,
                                       const int ic, const int t);
        void matrixAssemblyAdjointDO(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal,
                                     const int ic, const int t);
        void boundaryIntegral(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal,
                              const int ic, const int ib);
        void boundaryInGaussIntegral(MatrixXd &Klocal, double (&dxdr2D)[2][2], const double weight,
                                     const int ii, const int jj);
        void updateVariables(std::string output, const int dim, const int t, const int loop);
        void adjointGaussIntegralLHS(DirectProblem &main, MatrixXd &Klocal, const double &f, const int &ii, const int &jj);
        void adjointGaussIntegralRHS(DirectProblem &main, VectorXd &Flocal, const double &f, const int &ii);

    private:
        void setVariablesZero(const int dim);
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

        VoxelVelocity vvox;

        double aCF, bCF1, bCF2, gCF;
        int loopMax;
        int nControlNodesInCell;

        std::vector<int> planeDir;
        std::vector<std::vector<std::vector<double>>> feedbackForce;
        std::vector<std::vector<std::vector<double>>> feedbackForceT;
        std::vector<std::vector<std::vector<double>>> gradWholeNode;
        std::vector<std::vector<std::vector<double>>> grad;
        std::vector<std::vector<std::vector<double>>> X;
        
        void initialize(Config &conf);
        void runSimulation();

        void output(const int loop);
        void guessInitialCondition();
        void compCostFunction();
        void GaussIntegralRegTerm1(std::vector<double> &N, std::vector<std::vector<double>> &dNdr,
                                   std::vector<std::vector<double>> &xCurrent, double &value, 
                                   const double weight, const int ic, const int t);
        void GaussIntegralRegTerm2(std::vector<double> &N, std::vector<std::vector<double>> &dNdr,
                                   std::vector<std::vector<double>> &xCurrent, double &value, 
                                   const double weight, const int ic, const int t);
        void compFeedbackForce();
        void compEdgeValue(std::vector<std::vector<std::vector<std::vector<double>>>> &vEX, const int t);
        void compInterpolatedFeeback(std::vector<std::vector<double>> &xCurrent, double (&feedback)[3], 
                                     std::vector<std::vector<std::vector<std::vector<double>>>> &vEX, 
                                     double (&point)[3]);
        void compTimeInterpolatedFeedbackForce();
        void feedbackGaussIntegral(std::vector<double> &N, double (&feedback)[3], 
                                   const double detJ, const double weight, const int ic, const int t);
        void compFeedbackForce2();
        void feedbackGaussIntegral2(std::vector<double> &N, std::vector<std::vector<double>> &xCurrent,
                                    std::vector<std::vector<double>> &velCurrent, const double detJ,
                                    const double weight, const int voxelId, const int cellId, const int t);

        void compOptimalCondition();
        void GaussIntegralOptimalConditionTerm1(std::vector<double> &N, std::vector<std::vector<double>> &dNdr, 
                                                std::vector<std::vector<double>> &xCurrent, 
                                                double (&value)[4][3], const double weight, 
                                                const int ic, const int t);
        void GaussIntegralOptimalConditionTerm2(std::vector<double> &N, std::vector<std::vector<double>> &dNdr, 
                                                std::vector<std::vector<double>> &xCurrent, 
                                                double (&value)[4][3], const double weight, 
                                                const int ic, const int t);
        void GaussIntegralOptimalConditionTerm3(std::vector<double> &N, std::vector<std::vector<double>> &dNdr, 
                                                std::vector<std::vector<double>> &xCurrent, 
                                                double (&value)[4][3], const double weight, 
                                                const int ic, const int t);
        double armijoCriteria(const double fk);
        void updataControlVariables(DirectProblem &main, const double alpha);
    
    private:

};

#endif