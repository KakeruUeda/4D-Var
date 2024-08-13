
#ifndef DIRECTPROBLEM_H
#define DIRECTPROBLEM_H

#include "Boundary.h"
#include "Config.h"
#include "Export.h"
#include "FEM.h"
#include "FileIO.h"
#include "Gauss.h"
#include "Grid.h"
#include "Import.h"
#include "MathCommon.h"
#include "MathTool.h"
#include "PetscSolver.h"
#include "ShapeFunction.h"
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

extern MyMPI mpi;

class DirectProblem : public virtual FEM
{
public:
  DirectProblem(Config &conf);
  ~DirectProblem()
  {
  }

  Application app;
  Grid grid;
  Dirichlet dirichlet;
  PetscSolver petsc;
  SnapShot snap;

  int dim, nOMP;
  std::string outputDir;

  /*
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
  */

  std::vector<double> vgp;
  std::vector<double> advgp;
  std::vector<std::vector<double>> dvgpdx;

  void initialize(Config &conf);
  void resize();
  void runSimulation();
  void outputDomain();
  void solveUSNS(Application &app);
  void solveUSNS(std::vector<std::map<int, std::vector<double>>> &vDirichletTmp,
                 std::vector<std::map<int, double>> &pDirichletTmp, std::vector<std::vector<double>> &v0Tmp);
  void compInitialCondition(std::vector<std::map<int, std::vector<double>>> &vDirichletTmp,
                            std::vector<std::map<int, double>> &pDirichletTmp);
  void matrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t);
  void assembleLocalMatrixAndVector(MatrixXd &Klocal, VectorXd &Flocal, MathTools3D &tools, double f,
                                    const int ic);
  void updateSolutionsVTI();
  void updateSolutionsVTI(const int t);
  void outputSolutions(const int t);
  void outputSolutionsVTI(const std::string &dir, const int t);
  void outputSolutionsVTU(const std::string &dir, const int t);
  void outputSolutionsVTI(const std::string &dir, const int t, const int loop);
  void outputSolutionsVTU(const std::string &dir, const int t, const int loop);
  void outputSolutionsBIN(const std::string &dir, const int t);

private:
  void setValuesInGaussIntegral(MathTools3D &tools, Gauss &g2, const double he, const int i1, const int i2,
                                const int i3, const int ic, const int t);
  void usnsGaussIntegralLHS(MatrixXd &Klocal, MathTools3D &tools, const double f, const int ii, const int jj);
  void usnsGaussIntegralRHS(VectorXd &Flocal, MathTools3D &tools, const double f, const int ii);
  void setVelocityValue(MathTools3D &tools, const int ic, const int t);
  void updateSolutions();
  void updateTimeSolutions(const int t);
  void setVariablesZero();
  void compVorticity(const int t);

  // add
  void solveNavierStokes();
  void solveFowardNavierStokes(Array2D<double> &X0, Array3D<double> &X);

  void updateInitialVelocity(Array2D<double> &X0);
};

#endif
