
#ifndef DIRECTPROBLEM_H
#define DIRECTPROBLEM_H

#include "Boundary.h"
#include "Config.h"
#include "DataGrid.h"
#include "Export.h"
#include "FEM.h"
#include "Gauss.h"
#include "Grid.h"
#include "Import.h"
#include "MathTool.h"
#include "MySolver.h"
#include "PetscSolver.h"
#include "Polynomial.h"
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

class DirectProblem : public FEM
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
  void setLocalValuesInGauss(MatrixXd &Klocal, VectorXd &Flocal, MathTools3D &tools, double f, const int ic);

  // ----- output for VTI ------
  void updateVelocityVTI();
  void updatePressureVTI();
  void updateVelocityVTI(const int t);
  void updatePressureVTI(const int t);
  // -----------------------
  void outputSolutions(const int t);

  // ----- output VTI ------
  void outputVelocityVTI(const std::string &dir, const int t);
  void outputPressureVTI(const std::string &dir, const int t);
  void outputVelocityVTI(const std::string &dir, const int t, const int loop);
  void outputPressureVTI(const std::string &dir, const int t, const int loop);
  // -----------------------

  // ----- output VTU -----
  void outputVelocityVTU(const std::string &dir, const int t);
  void outputPressureVTU(const std::string &dir, const int t);
  void outputVelocityVTU(const std::string &dir, const int t, const int loop);
  void outputPressureVTU(const std::string &dir, const int t, const int loop);
  void outputTemporaryVelocityVTU(const std::string &dir, const int t);
  void outputTemporaryPressureVTU(const std::string &dir, const int t);
  // ----------------------

  // ----- output BIN -----
  void outputVelocityBIN(const std::string &dir, const int t);
  void outputPressureBIN(const std::string &dir, const int t);
  // ----------------------

private:
  void resizeVar();
  void initializeVarZero();
  void settingInGauss(MathTools3D &tools, Gauss &g2, const double he, const int i1, const int i2, const int i3,
                      const int ic, const int t);
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
  void solveNavierStokes(const int stepMax, std::vector<std::array<double, 2>> &velArr);
  void solveNavierStokes(std::vector<std::array<double, 2>> &velArr,
                         std::vector<std::map<int, std::vector<double>>> &vectorVelocitySet);

  void solve_NavierStokes_polynmial_BCs(DataGrid &data, ControlBoundary &cb, const int c_cycle,
                                        std::vector<std::map<int, std::vector<double>>> &vectorVelocitySet);
  void solve_NavierStokes_polynmial_BCs(DataGrid &data, ControlBoundary &cb, const int t_max);
  void solve_NavierStokes_polynmial_BCs(DataGrid &data, ControlBoundary &cb,
                                        std::vector<std::map<int, std::vector<double>>> &vectorVelocitySet);
  void solve_NavierStokes_optimized_BCs();
  double interpolated_velocity(Array2D<double> &vel_opt, double px, double pz, int dir);
  void interpolateGrid(DataGrid &data, const VectorXd &coeff, std::vector<double> &x_new, std::vector<double> &y_new,
                       std::vector<double> &z_new);


};

#endif
