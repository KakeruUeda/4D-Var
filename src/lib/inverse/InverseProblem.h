/**
 * @file InverseProblem.h
 * @author K.Ueda
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
#include "Adjoint.h"
#include "Import.h"

extern MyMPI mpi;

struct EstimatedVariable
{
public:
  std::vector<double> u, v, w;
};

struct CostFunction
{
  double term1, term2, term3, term4, term5; // Reg for X
  double term6, term7;                      // Reg for X0
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
  ~InverseProblem() {}

  int dim, nOMP;
  std::string outputDir;

  Application app;
  DataGrid data;

  DirectProblem main;
  Adjoint adjoint;
  CostFunction costFunction;

  VoxelVelocity vvox;

  double aCF, bCF, gCF;
  int loopMax;
  int outputItr;

  double alphaX0, alphaX;

  bool isConverged_X, isConverged_X0;

  std::vector<int> planeDir;
  std::vector<std::vector<std::vector<double>>> feedbackForce;
  std::vector<std::vector<std::vector<double>>> feedbackForceT;
  std::vector<std::vector<std::vector<double>>> gradWholeNode;
  std::vector<std::vector<double>> gradInitVel;
  std::vector<std::vector<std::vector<double>>> grad;
  std::vector<std::vector<std::vector<double>>> X;
  std::vector<std::vector<double>> X0;
  std::vector<std::vector<std::vector<double>>> Xvti;
  std::vector<std::vector<double>> X0vti;

  void initialize(Config &conf);
  void runSimulation();
  void guessInitialCondition();
  void compCostFunction();
  void GaussIntegralRegTerm2(Function &func, double &value, const int ic, const int t);
  void GaussIntegralRegTerm3(Function &func, double &value, const int ic, const int t);
  void GaussIntegralRegTerm4(Function &func, double &value, const int ic, const int t);
  void GaussIntegralRegTerm5(Function &func, double &value, const int ic, const int t);
  void GaussIntegralRegTerm6(Function &func, double &value, const int ic);
  void GaussIntegralRegTerm7(Function &func, double &value, const int ic);
  void compFeedbackForce();
  void compInterpolatedFeeback(double (&feedback)[3], double (&point)[3]);
  void compTimeInterpolatedFeedbackForce();
  void feedbackGaussIntegral(Function &func, double (&feedback)[3], const int ic, const int t);
  void compOptimalCondition();
  void GaussIntegralOptimalConditionXTerm1(Function &func, std::vector<std::vector<double>> &value, const int ic, const int t);
  void GaussIntegralOptimalConditionXTerm2(Function &func, std::vector<std::vector<double>> &value, const int ic, const int t);
  void GaussIntegralOptimalConditionXTerm3(Function &func, std::vector<std::vector<double>> &value, const int ic, const int t);
  void GaussIntegralOptimalConditionXTerm4(Function &func, std::vector<std::vector<double>> &value, const int ic, const int t);
  void GaussIntegralOptimalConditionXTerm5(Function &func, std::vector<std::vector<double>> &value, const int ic, const int t);
  void GaussIntegralOptimalConditionX0Term1(Function &func, std::vector<std::vector<double>> &value, const int ic);
  void GaussIntegralOptimalConditionX0Term2(Function &func, std::vector<std::vector<double>> &value, const int ic);
  void GaussIntegralOptimalConditionX0Term3(Function &func, std::vector<std::vector<double>> &value, const int ic);
  double armijoCriteria(const double fk);
  double armijoCriteriaX_tmp(const double fk);
  double armijoCriteriaX0_tmp(const double fk);
  void armijoCriteriaX(const double fk);
  void armijoCriteriaX0(const double fk);
  void updataControlVariables(DirectProblem &main);
  void updataControlVariables(DirectProblem &main, const double alphaX, const double alphaX0);
  void setValue(Function &func, const int ic);

private:
  void assembleFeedbackForce(Function &func, const int ic, const int t);
  bool checkConvergence(std::ofstream &cf, const int loop);
  void outputFowardSolutions(const int loop);
  void outputAdjointSolutions(const int loop);
  void outputFeedbackForce(const int loop);
  void outputControlVariables(const int loop);
  void outputVelocityData(const int loop);
  void outputVelocityBIN(const int loop);
  void outputOptimizedVariables();

  void updateControlVariablesVTI();
};

#endif