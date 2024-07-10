/**
 * @file DirectProblem.h
 * @author K.Ueda
 * @date Jun, 2024
 */

#ifndef DIRECTPROBLEM_H
#define DIRECTPROBLEM_H

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
#include <algorithm>
#include "Grid.h"
#include "Boundary.h"
#include "PetscSolver.h"
#include "Config.h"
#include "Gauss.h"
#include "ShapeFunction.h"
#include "MathFEM.h"
#include "Tool.h"
#include "Function.h"

extern MyMPI mpi;

class DirectProblem
{
    public:
        DirectProblem(Config &conf);
        ~DirectProblem(){}

        int dim, nOMP;
        std::string outputDir;

        Application app;
        Grid grid;
        PetscSolver petsc;
        SnapShot snap;

        int IU, IV, IW, IP;
        int ILU, ILV, ILW;
        int JU, JV, JW, JP;
        int JLU, JLV, JLW;

        // Pysical parameter
        double Re, rho, mu, nu;

        // Time parameter
        double dt;
        int timeMax;
        int pulsatileFlow;
        int pulseBeginItr;
        double T;

        double NRtolerance;

        // Darcy parameter
        double alpha, resistance;

        double tau;
        std::vector<double> vgp;
        std::vector<double> advgp;
        std::vector<std::vector<double>> dvgpdx;

        void initialize(Config &conf);
        void runSimulation();
        void outputDomain();
        void solveUSNS(Application &app);
        void solveUSNS(std::vector<std::map<int, std::vector<double>>> &vDirichletTmp,
                       std::vector<std::map<int, double>> &pDirichletTmp);
        void compInitialCondition(std::vector<std::map<int, std::vector<double>>> &vDirichletTmp,
                                  std::vector<std::map<int, double>> &pDirichletTmp);
        void matrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, Function &func, const int ic, const int t);

    private:
        void mainGaussIntegralLHS(MatrixXd &Klocal, Function &func, const double f, const int ii, const int jj);
        void mainGaussIntegralRHS(VectorXd &Flocal, Function &func, const double f, const int ii);
        void setVelocityValue(Function &func, const int ic, const int t);
        void updateVariables(const int t);
        void assignTimeVariables(const int t);
        void outputSolution(const int t);
        void setVariablesZero();
        void updateRowIndex(const int ii, const int ic);
        void updateColumnIndex(const int ii, const int ic);
};

#endif

