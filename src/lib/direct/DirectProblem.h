
#ifndef DIRECTPROBLEM_H
#define DIRECTPROBLEM_H

#include "Boundary.h"
#include "Config.h"
#include "Export.h"
#include "FEM.h"
#include "Gauss.h"
#include "Grid.h"
#include "Import.h"
#include "MathTool.h"
#include "PetscSolver.h"
#include "ShapeFunction.h"
#include "Tool.h"
#include "MySolver.h"
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

  void initialize(Config &conf);
  void runSimulation();
  void outputDomain();
  void matrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t);
  void setLocalValuesInGauss(MatrixXd &Klocal, VectorXd &Flocal, MathTools3D &tools, double f,
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
  void resizeVar();
  void initializeVarZero();
  void settingInGauss(MathTools3D &tools, Gauss &g2, const double he, const int i1, const int i2,
                                const int i3, const int ic, const int t);
  void usnsGaussIntegralLHS(MatrixXd &Klocal, MathTools3D &tools, const double f, const int ii, const int jj);
  void usnsGaussIntegralRHS(VectorXd &Flocal, MathTools3D &tools, const double f, const int ii);
  void setVelocityValue(MathTools3D &tools, const int ic, const int t);
  void updateSolutions();
  void updateTimeSolutions(const int t);
  void setVariablesZero();
  void compVorticity(const int t);

public:
  // add
  void solveNavierStokes();
  void solveNavierStokes(Array2D<double> &X0, Array3D<double> &X);

  void updateInitialVelocity(Array2D<double> &X0);
  void solveNaveirStokes(const int stepMax);
};

#endif
