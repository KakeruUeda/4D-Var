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
#include "Export.h"

class Grid
{
public:
  Grid() {}
  Grid(Config &conf);
  virtual ~Grid() {}

  GridType gridType;
  Cell cell;
  Node node;

  int nx, ny, nz;
  double lx, ly, lz;
  double dx, dy, dz;

  int dim;
  int extractFluid;
  int nodeStart, nodeEnd;
  int rowStart, rowEnd;

  int nNodesGlobal, nCellsGlobal, nDofsGlobal;
  int nNodesLocal, nCellsLocal, nDofsLocal;

  std::vector<int> vecFluidUniqueNodes;

  void prepareMatrix(Dirichlet &dirichletBC, PetscSolver &petsc, std::string outputDir, const int timeMax);
  void initialSetting(Dirichlet &dirichletBC);
  void serialSetting();
  void parallelSetting(Dirichlet &dirichletBC);
  void collectNodeMap();
  void distributeToLocalDofs();
  void getNodeNewMap(Dirichlet &dirichletBC);
  void getCellNewMap();
  void petscMatrixSetting(PetscSolver &petsc);

private:
};

#endif // GRID_H
