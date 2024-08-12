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

  void setStructuredGrid(
      const int nxCells, const int nyCells, const int nzCells,
      const int nxNodes, const int nyNodes, const int nzNodes,
      const double dx, const double dy, const double dz,
      const int nNodesInCell, const int dim,
      Cell &cell, Node &node);

  void prepareMatrix(PetscSolver &petsc, std::string outputDir, const int timeMax);
  void setForSerial();
  void divideWholeGrid();
  void distributeToLocal(const int timeMax);

private:
  int structuredGridNodeSet(
      const int nxNodes, const int nyNodes, const int nzNodes,
      const int i, const int j, const int k, const int p);
  double structuredGridCoordinateSet(
      const double dx, const double dy, const double dz,
      const int i, const int j, const int k, const int d);
};

#endif