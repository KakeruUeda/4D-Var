/**
 * @file Grid.h
 * @author K.Ueda
 * @date Jun, 2024
 */

#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <set>
#include "metis.h"
#include "Gauss.h"
#include "ShapeFunction.h"
#include "Array.h"
#include "Boundary.h"
#include "Config.h"

class Grid
{
public:
  Grid() {}
  Grid(Config &conf);
  virtual ~Grid() {}

  GridType gridType;
  Cell cell;
  Node node;
  DirichletBoundary dirichlet;

  int nx, ny, nz;
  double lx, ly, lz;
  double dx, dy, dz;

  int dim;
  int extractFluid;
  int rowStart, rowEnd;

  int nNodesGlobal, nCellsGlobal, nDofsGlobal;
  int nNodesLocal, nCellsLocal, nDofsLocal;

  std::vector<int> vecFluidUniqueNodes;

  void prepareMatrix(PetscSolver &petsc, std::string outputDir, const int timeMax);
  void setForSerial();
  void distributeToLocal(const int timeMax);

private:
};

#endif