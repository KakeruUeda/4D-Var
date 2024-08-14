/**
 * @file Boundary.h
 * @author K.Ueda
 * @date May, 2024
 */

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Array.h"
#include "Cell.h"
#include "Config.h"
#include "Node.h"
#include "PetscSolver.h"
#include <iostream>
#include <memory>
#include <vector>

class Boundary
{
public:
  virtual void assignBCs(Node &node, const int t) = 0;
  virtual ~Boundary()
  {
  }
};

class Dirichlet : public Boundary
{
public:
  Dirichlet()
  {
  }

  Array1D<double> initialValues;
  Array1D<double> values;

  void initialize(Config &conf);
  void getNewArray(std::vector<int> mapNew);
  void setValuesZero(int n);
  void assignBCs(Node &node, const int t) override;
  void assignPulsatileBCs(const double pulse, const int nDofsGlobal);
  void applyBCs(Cell &cell, PetscSolver &petsc);

  void updateValues(Array3D<double> &X, const int t);

public:
  std::map<int, std::vector<double>> velocitySet;
  std::map<int, double> pressureSet;

  std::map<int, std::vector<double>> velocitySetNew;
  std::map<int, double> pressureSetNew;
};

class Neumann : public Boundary
{
public:
  Neumann(double gradient) : gradient(gradient)
  {
  }

private:
  double gradient;
};

class DirichletBoundary
{
public:
  DirichletBoundary()
  {
  }
  DirichletBoundary(Config &conf)
  {
  }
  virtual ~DirichletBoundary()
  {
  }

  int nNodesVelocity, nNodesPressure, nControlNodesInCell;
  int nControlCellsGlobal, nControlNodesGlobal;

  std::vector<double> dirichletBCsValue;
  std::vector<double> dirichletBCsValueNew;
  std::vector<double> dirichletBCsValueInit;
  std::vector<double> dirichletBCsValueNewInit;

  std::vector<std::map<int, std::vector<double>>> vDirichlet;
  std::vector<std::map<int, std::vector<double>>> vDirichletWall;
  std::vector<std::map<int, double>> pDirichlet;
  std::vector<std::map<int, std::vector<double>>> vDirichletNew;
  std::vector<std::map<int, std::vector<double>>> vDirichletWallNew;
  std::vector<std::map<int, double>> pDirichletNew;

  std::vector<int> controlBoundaryMap;
  std::vector<int> controlCellMap;
  std::vector<std::vector<int>> controlNodeInCell;

  std::vector<bool> isBoundaryEdge;

  void initialize(Config &conf);
  void initializeAdjoint(Config &conf);

  void assignDirichletBCs(std::vector<std::map<int, std::vector<double>>> &vDirichletNew,
                          std::vector<std::map<int, double>> &pDirichletNew, Node &node, int &dim, const int t);
  void assignConstantDirichletBCs(std::vector<std::map<int, std::vector<double>>> &vDirichletNew,
                                  std::vector<std::map<int, double>> &pDirichletNew, Node &node, int &dim, const int t);
  void assignPulsatileBCs(const int t, const double dt, const double T, const int pulseBeginItr, const int nDofsGlobal);
  void applyDirichletBCs(Cell &cell, PetscSolver &petsc);
  void applyDirichletBCsAdjoint(Cell &cell, PetscSolver &petsc);

private:
};

#endif
