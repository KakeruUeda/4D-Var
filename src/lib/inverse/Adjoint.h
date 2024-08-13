/**
 * @file InverseProblem.h
 * @author K.Ueda
 * @date July, 2024
 */

#ifndef ADJOINT_H
#define ADJOINT_H

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
#include "Tool.h"
#include "Spline.h"
#include "ShapeFunction.h"
#include "DirectProblem.h"
#include "Function.h"

class Adjoint
{
public:
  Adjoint(Config &conf) : grid(conf), dim(conf.dim), planeDir(conf.planeDir), timeMax(conf.timeMax),
                          rho(conf.rho), mu(conf.mu), dt(conf.dt), alpha(conf.alpha), resistance(conf.resistance)
  {
  }

  Grid grid;
  PetscSolver petsc;

  int dim;
  int timeMax;
  double dt;
  double rho, mu, nu, Re;
  double alpha, resistance;

  std::string outputDir;

  int IU, IV, IW, IP;
  int ILU, ILV, ILW;
  int JU, JV, JW, JP;
  int JLU, JLV, JLW;

  double tau;
  std::vector<double> vgp;
  std::vector<double> advgp;
  std::vector<std::vector<double>> dvgpdx;

  std::vector<double> vk, vk1, vk2;
  std::vector<double> advk1, advk2, advk3;
  std::vector<double> dpkdx, dpk1dx, dpk2dx;
  std::vector<double> wk, wk1, wk2;
  std::vector<std::vector<double>> dvkdx, dvk1dx, dvk2dx;
  std::vector<std::vector<double>> dwkdx, dwk1dx, dwk2dx;
  std::vector<double> dqkdx, dqk1dx, dqk2dx;

  std::vector<int> planeDir;

  void solveAdjoint(DirectProblem &main, std::string outputDir,
                    std::vector<std::vector<std::vector<double>>> &feedbackForceT);
  void setValue(DirectProblem &main, Function &func, const int ic, const int t);
  void matrixAssemblyAdjoint(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal,
                             Function &func, const int ic, const int t);
  void boundaryIntegral(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal,
                        Function &func, const int ic, const int ib);
  void boundaryInGaussIntegral(MatrixXd &Klocal, Function &func, const int ii, const int jj);
  void setValuesInGaussIntegral(DirectProblem &main, Function &func, Gauss &g2, const double he,
                                const int i1, const int i2, const int i3, const int ic, const int t);
  void adjointGaussIntegralLHS(DirectProblem &main, MatrixXd &Klocal, Function &func,
                               const double f, const int ii, const int jj);
  void adjointGaussIntegralRHS(DirectProblem &main, VectorXd &Flocal, Function &func,
                               const double f, const int ii);
  void updateSolutionsVTI();
  void updateSolutionsVTI(const int t);
  void outputSolutionsVTU(const std::string &dir, const int t);
  void outputSolutionsVTI(const std::string &dir, const int t);
  void outputSolutionsVTU(const std::string &dir, const int t, const int loop);
  void outputSolutionsVTI(const std::string &dir, const int t, const int loop);

private:
  void setVariablesZero(const int dim);
  void updateRowIndex(const int ii, const int ic);
  void updateColumnIndex(const int ii, const int ic);
  void updateRowIndexPlane(const int jj, const int ic);
  void updateColumnIndexPlane(const int jj, const int ic);
  void updateSolutions();
  void updateTimeSolutions(const int t);
};

#endif
