/**
 * @file Boundary.h
 * @author K.Ueda
 * @date May, 2024
 */

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <iostream>
#include <vector>
#include <memory>
#include "Config.h"
#include "Array.h"
#include "Node.h"
#include "Cell.h"
#include "PetscSolver.h"

class Boundary
{
public:
  virtual void assignBCs(Node &node, const int t) = 0;
  virtual ~Boundary() {}
};

class Dirichlet : public Boundary
{
public:
  Dirichlet(){}

  void assignBCs(Node &node, const int t) override;

private:
  std::map<int, double[3]> vDirichletSet;
  std::map<int, double> pDirichletSet;
  Array1D<double> dirichletValue;
};

class Neumann : public Boundary
{
public:
  Neumann(double gradient) : gradient(gradient) {}

private:
  double gradient;
};


class DirichletBoundary
{
public:
  DirichletBoundary() {}
  DirichletBoundary(Config &conf) {}
  virtual ~DirichletBoundary() {}

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
                          std::vector<std::map<int, double>> &pDirichletNew, Node &node,
                          int &dim, const int t);
  void assignConstantDirichletBCs(std::vector<std::map<int, std::vector<double>>> &vDirichletNew,
                                  std::vector<std::map<int, double>> &pDirichletNew, Node &node,
                                  int &dim, const int t);
  void assignPulsatileBCs(const int t, const double dt, const double T,
                          const int pulseBeginItr, const int nDofsGlobal);
  void applyDirichletBCs(Cell &cell, PetscSolver &petsc);
  void applyDirichletBCsAdjoint(Cell &cell, PetscSolver &petsc);

private:
};

class StructuredBoundaryFace
{
public:
  int size;
  StructuredBoundaryFace(std::string face) : bdFaceStr(face) {};
  ~StructuredBoundaryFace() {};

  std::vector<int> node;
  std::vector<std::string> dirichletType;
  std::vector<std::vector<double>> dirichletValue;

  int getNodeSize()
  {
    return node.size();
  };

  void setSize(int n)
  {
    size = n;
  };

  std::string bdFaceStr;
  void setNodesOnBoundaryFace(int nxNodes, int nyNodes, int nzNodes);

  void setDirichletInfo(std::vector<std::string> bdType,
                        std::vector<std::vector<double>> bdValue,
                        int dim, int bdIndex);

private:
};

#endif
