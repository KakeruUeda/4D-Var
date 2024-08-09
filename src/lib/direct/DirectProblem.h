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
#include "Tool.h"
#include "Function.h"
#include "VTK.h"
#include "MathTool.h"

extern MyMPI mpi;

class DirectProblem
{
public:
  DirectProblem(Config &conf);
  ~DirectProblem() {}

  Application app;
  Grid grid;
  PetscSolver petsc;
  SnapShot snap;

  int dim, nOMP;
  std::string outputDir;

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
                 std::vector<std::map<int, double>> &pDirichletTmp, std::vector<std::vector<double>> &v0Tmp);
  void compInitialCondition(std::vector<std::map<int, std::vector<double>>> &vDirichletTmp,
                            std::vector<std::map<int, double>> &pDirichletTmp);
  void matrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, MathTools3D &math, const int ic, const int t);

  void updateSolutionsVTI();
  void updateSolutionsVTI(const int t);
  void outputSolutionsVTI(const std::string &dir, const int t);
  void outputSolutionsVTU(const std::string &dir, const int t);
  void outputSolutionsVTI(const std::string &dir, const int t, const int loop);
  void outputSolutionsVTU(const std::string &dir, const int t, const int loop);

private:
  void setValuesInGaussIntegral(MathTools3D &math, Gauss &g2, const double he,
                                const int i1, const int i2, const int i3, const int ic, const int t);
  void usnsGaussIntegralLHS(MatrixXd &Klocal, MathTools3D &math, const double f, const int ii, const int jj);
  void usnsGaussIntegralRHS(VectorXd &Flocal, MathTools3D &math, const double f, const int ii);
  void setVelocityValue(MathTools3D &math, const int ic, const int t);
  void updateSolutions();
  void updateTimeSolutions(const int t);
  void setVariablesZero();
  void updateRowIndex(const int ii, const int ic);
  void updateColumnIndex(const int ii, const int ic);
  void compVorticity(const int t);
};

#endif
