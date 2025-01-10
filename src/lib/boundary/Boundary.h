/**
 * @file Boundary.h
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Array.h"
#include "Cell.h"
#include "Config.h"
#include "Node.h"
#include "PetscSolver.h"
#include "Spline.h"
#include <iostream>
#include <memory>
#include <vector>
#include <set>

class ControlBoundary
{
public:
  ControlBoundary() 
  {
  }
  ControlBoundaryFace inletFace;
  std::vector<bool> isBoundaryEdge;
  std::vector<int> CBNodeMap;
  std::vector<int> CBEdgeNodeMap;
  std::vector<int> CBCellMap;
  std::vector<std::vector<int>> CBNodeMapInCell;

  void initialize(Config &conf);
private:
};


class Boundary
{
public:
  virtual void assignBCs(Node &node) = 0;
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

  struct Velocity
  {
    double time;
    double value;
  };

  Array1D<double> initialValues;
  Array1D<double> values;

  void initialize(Config &conf);
  void getNewArray(std::vector<int> mapNew);
  void setValuesZero(int n);
  void assignBCs(Node &node) override;
  void assignPulsatileBCs(const double pulse);
  void applyBCs(Cell &cell, PetscSolver &petsc);
  void updateValues(Array3D<double> &X, const int t);
  void eraseControlNodes(Cell &cell, ControlBoundary &cb);
  double comp_pulse(const double step);
  double comp_pulse2(double timeNow, std::vector<std::array<double, 2>> &velArr);

public:
  std::map<int, std::vector<double>> velocitySet;
  std::map<int, double> pressureSet;

  std::map<int, std::vector<double>> velocitySetInit;

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

#endif
