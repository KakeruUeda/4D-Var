/**
 * @file InverseProblem.h
 * @author K.Ueda
 * @date July, 2024
 */

#ifndef INVERSEPROBLEM_H
#define INVERSEPROBLEM_H

#include "Adjoint.h"
#include "Boundary.h"
#include "Config.h"
#include "DataGrid.h"
#include "DirectProblem.h"
#include "Function.h"
#include "Gauss.h"
#include "Grid.h"
#include "Import.h"
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

extern MyMPI mpi;

struct CostFunction
{
  double term1, term2, term3, term4, term5;  // Reg term for X
  double term6, term7;                       // Reg term for X0
  double total;
  std::vector<double> history;
  void sum()
  {
    total = term1 + term2 + term3 + term4 + term5 + term6 + term7;
  }
};

class InverseProblem
{
public:
  InverseProblem(Config &conf);
  ~InverseProblem()
  {
  }

  int dim, nOMP;
  std::string outputDir;

  Application app;

  // Do not change the order of the class declarations
  DirectProblem main;
  Adjoint adjoint;
  DataGrid data;
  ControlBoundary inletCB;
  CostFunction costFunction;

  MathTools2D mt2d;
  MathTools3D mt3d;

  VoxelVelocity vvox;

  double aCF, bCF, gCF;
  int loopMax;
  int outputItr;

  double alphaX0, alphaX;
  bool isConverged_X, isConverged_X0;

  std::vector<int> planeDir;

  Array2D<double> gradX0;
  Array3D<double> gradX;
  Array3D<double> X;
  Array2D<double> X0;
  Array3D<double> Xvti;
  Array2D<double> X0vti;

  void initialize(Config &conf);
  void resize();
  void initializeVarZero();
  void runSimulation();

private:
  void compInitialOptimalVelocityField();

  // Cost function
  void compCostFunction();
  void RegTerm2_inGaussIntegral(double &value, const int nc, const int ic, const int t);
  void RegTerm3_inGaussIntegral(double &value, const int nc, const int ic, const int t);
  void RegTerm4_inGaussIntegral(double &value, const int nc, const int ic, const int t);
  void RegTerm5_inGaussIntegral(double &value, const int nc, const int ic, const int t);
  void RegTerm6_inGaussIntegral(double &value, const int ic);
  void RegTerm7_inGaussIntegral(double &value, const int ic);

  // Feedback force
  void compFeedbackForce();
  void compInterpolatedFeeback(double (&feedback)[3], double (&point)[3], const int t);
  void compTimeInterpolatedFeedbackForce();
  void assembleFeedbackForce(const int ic, const int t);
  void feedbackGaussIntegral(double (&feedback)[3], const int ic, const int t);
  
  // Optimal condition
  void compOptimalCondition();
  void OptCondX_Term1_inGaussIntegral(std::vector<std::vector<double>> &value, const int nc, const int ic, const int t);
  void OptCondX_Term2_inGaussIntegral(std::vector<std::vector<double>> &value, const int nc, const int ic, const int t);
  void OptCondX_Term3_inGaussIntegral(std::vector<std::vector<double>> &value, const int nc, const int ic, const int t);
  void OptCondX_Term4_inGaussIntegral(std::vector<std::vector<double>> &value, const int nc, const int ic, const int t);
  void OptCondX_Term5_inGaussIntegral(std::vector<std::vector<double>> &value, const int nc, const int ic, const int t);
  void OptCondX0_Term1_inGaussIntegral(std::vector<std::vector<double>> &value, const int ic);
  void OptCondX0_Term2_inGaussIntegral(std::vector<std::vector<double>> &value, const int ic);
  void OptCondX0_Term3_inGaussIntegral(std::vector<std::vector<double>> &value, const int ic);
  void setValue(const int ic);

  // Armijo criteria
  void armijoCriteriaX(const double fk);
  void armijoCriteriaX0(const double fk);

  // Output functions
  void outputFowardSolutions(const int loop);
  void outputAdjointSolutions(const int loop);
  void outputFeedbackForce(const int loop);
  void outputControlVariables(const int loop);
  void outputVelocityData(const int loop);
  void outputVelocityBIN(const int loop);
  void outputOptimizedVariables();

  void updateControlVariablesVTI();
  bool checkConvergence(std::ofstream &cf, const int loop);

};

#endif