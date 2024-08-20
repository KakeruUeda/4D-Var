/**
 * @file DataGrid.h
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef DATAGRID_H
#define DATAGRID_H

#include "Cell.h"
#include "Config.h"
#include "Function.h"
#include "Gauss.h"
#include "Grid.h"
#include "MathTool.h"
#include "PetscSolver.h"
#include "ShapeFunction.h"
#include "Tool.h"
#include <iostream>

struct VoxelInfo
{
  double center[3];
  double minX, minY, minZ;
  double maxX, maxY, maxZ;
  vector<int> cells;
  Array2D<double> v_cfd, v_mri, v_err;  // (time, dim)
};

class DataGrid
{
public:
  DataGrid(Config &conf, Grid &grid, SnapShot &snap);

  VoxelVelocity vvox;
  MathTools3D mt3d;
  Gauss gauss;
  Grid &grid;
  SnapShot &snap;

  int nxData, nyData, nzData;
  double dxData, dyData, dzData;
  double lxData, lyData, lzData;
  double xOrigin, yOrigin, zOrigin;

  int nDataCellsGlobal;
  int nDataNodesGlobal;
  int nDataNodesInCell;

  Array3D<VoxelInfo> voxel;

  void initialize(Config &conf);

  void setVoxelCenters();
  void setVoxelBoundaries();
  void collectCellsInVoxel();
  void collectCellsInCircle(const int radious);
  bool isCellCenterIncludedInVoxel(const int iv, const int ic);
  bool isCellBoundaryIncludedInVoxel(const int iv, const int ic);
  bool isCellIncludedInCircle(const double radious, const int iv, const int ic);

  void average(const int iv, const int t);
  void weightedAverage(const int iv, const int t);
  void interpolate();

  void importDAT(const std::string &filename, const int step);
  void exportDAT(const std::string &filename, const int step);
  void exportVTI(const std::string &filename, const int t);

private:
  double coeff;
  double weightIntegral;

  Array2D<double> velCurrent;
  Array1D<double> smoothing;

  void compSmoothing();
};

#endif
