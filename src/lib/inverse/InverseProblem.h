/**
 * @file InverseProblem.h
 * @author K.U.
 * @date July, 2024
 */

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
#include "Function.h"

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

        int IU, IV, IW, IP;
        int ILU, ILV, ILW;
        int JU, JV, JW, JP;
        int JLU, JLV, JLW;

        double tau;
        std::vector<double> vgp;
        std::vector<double> advgp;
        std::vector<std::vector<double>> dvgpdx;

        std::vector<double> vk, vk1, vk2;
        std::vector<double> advk1, advk2;
        std::vector<double> wk1, wk2;
        std::vector<std::vector<double>> dvkdx, dvk1dx, dvk2dx;
        std::vector<std::vector<double>> dwk1dx, dwk2dx;
        
        std::vector<int> planeDir;

        void solveAdjoint(DirectProblem &main, std::string outputDir,
                             std::vector<std::vector<std::vector<double>>> &feedbackForceT,
                             const int nData, const int loop);
        void setValue(DirectProblem &main, Function &func, const int ic, const int t);
        void matrixAssemblyAdjoint(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal,
                                   Function &func, const int ic, const int t);
        void boundaryIntegral(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal,
                              Function &func, const int ic, const int ib);
        void boundaryInGaussIntegral(MatrixXd &Klocal, Function &func, const int ii, const int jj);
        void updateVariables(std::string output, const int dim, const int t, const int loop);
        void adjointGaussIntegralLHS(DirectProblem &main, MatrixXd &Klocal, Function &func, 
                                     const double f, const int ii, const int jj);
        void adjointGaussIntegralRHS(DirectProblem &main, VectorXd &Flocal, Function &func, 
                                     const double f, const int ii);

    private:
        void setVariablesZero(const int dim);
        void updateRowIndex(const int ii, const int ic);
        void updateColumnIndex(const int ii, const int ic);
        void updateRowIndexPlane(const int jj, const int ic);
        void updateColumnIndexPlane(const int jj, const int ic);
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
        void GaussIntegralRegTerm2(Function &func, double &value, const int ic, const int t);
        void GaussIntegralRegTerm3(Function &func, double &value, const int ic, const int t);
        void compFeedbackForce();
        void compInterpolatedFeeback(double (&feedback)[3], double (&point)[3]);
        void compTimeInterpolatedFeedbackForce();
        void feedbackGaussIntegral(Function &func, double (&feedback)[3], const int ic, const int t);
        void compOptimalCondition();
        void GaussIntegralOptimalConditionTerm1(Function &func, double (&value)[4][3], const int ic, const int t);
        void GaussIntegralOptimalConditionTerm2(Function &func, double (&value)[4][3], const int ic, const int t);
        void GaussIntegralOptimalConditionTerm3(Function &func, double (&value)[4][3], const int ic, const int t);
        double armijoCriteria(const double fk);
        void updataControlVariables(DirectProblem &main, const double alpha);
    
    private:
        void assembleFeedbackForce(Function &func, const int ic, const int t);

};

#endif