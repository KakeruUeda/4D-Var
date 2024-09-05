/**
 * @file InverseProblem.h
 * @author K.Ueda
 * @date July, 2024
 */

#ifndef ADJOINT_H
#define ADJOINT_H

#include "Boundary.h"
#include "Config.h"
#include "DirectProblem.h"
#include "Function.h"
#include "Gauss.h"
#include "Grid.h"
#include "PetscSolver.h"
#include "ShapeFunction.h"
#include "Spline.h"
#include "Tool.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <set>
#include <string>
#include <sys/stat.h>

class Adjoint : public FEM
{
public:
  Adjoint(Config &conf)
      : FEM(conf), grid(conf), dim(conf.dim), planeDir(conf.planeDir)
  {
  }

  Grid grid;
  Dirichlet dirichlet;
  PetscSolver petsc;

  MathTools3D mt3d;
  MathTools2D mt2d;

  int dim;

  std::string outputDir;

  std::vector<double> vgp;
  std::vector<double> advgp;
  std::vector<std::vector<double>> dvgpdx;

  std::vector<int> planeDir;

  Array3D<double> feedbackForce;
  Array3D<double> feedbackForceT;

  void solveAdjoint(DirectProblem &main, ControlBoundary &cb);
  void setValue(DirectProblem &main, const int ic, const int t);
  void matrixAssemblyAdjoint(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t);
  void boundaryIntegral(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal, ControlBoundary &cb, const int ic, const int ib);
  void boundaryInGaussIntegral(MatrixXd &Klocal, const int ii, const int jj);
  void setValuesInGaussIntegral(DirectProblem &main, Gauss &g2, const double he, const int i1, const int i2, const int i3, const int ic, const int t);
  void adjointGaussIntegralLHS(DirectProblem &main, MatrixXd &Klocal, const double f, const int ii, const int jj);
  void adjointGaussIntegralRHS(DirectProblem &main, VectorXd &Flocal, const double f, const int ii);
  void updateSolutionsVTI();
  void updateSolutionsVTI(const int t);
  void outputSolutionsVTU(const std::string &dir, const int t);
  void outputSolutionsVTI(const std::string &dir, const int t);
  void outputSolutionsVTU(const std::string &dir, const int t, const int loop);
  void outputSolutionsVTI(const std::string &dir, const int t, const int loop);

private:
  void setVariablesZero(const int dim);
  void updateSolutions();
  void updateTimeSolutions(const int t);
};

#endif
