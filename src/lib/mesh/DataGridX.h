/**
 * @file DataGrid.h
 * @author K.Ueda
 * @date July, 2024
 */

#ifndef DATAGRIDX_H
#define DATAGRIDX_H

#include "Cell.h"
#include "Config.h"
#include "Function.h"
#include "Gauss.h"
#include "Grid.h"
#include "MathCommon.h"
#include "MathTool.h"
#include "PetscSolver.h"
#include "ShapeFunction.h"
#include "Tool.h"
#include <iostream>

struct Vec3d
{
  double x, y, z;
};

struct VoxelInfoX
{
  double center[3];
  double minX, minY, minZ;
  double maxX, maxY, maxZ;
  vector<int> cells;
  Array2D<double> v_cfd, v_mri, v_err;
};

class DataGridX
{
public:
  DataGridX(Config &conf, Grid &grid, SnapShot &snap);

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

  Array3D<VoxelInfoX> voxel;

  void setVoxelCenters();
  void setVoxelBoundaries();
  void collectCellsInVoxel();
  void collectCellsInCircle(const int radious);
  bool isCellsIncludedInVoxel(const int iv, const int ic);
  bool isCellsIncludedInCircle(const double radious, const int iv, const int ic);

  void average(const int iv, const int t);
  void interpolate();

private:
  void compSmoothing();
  double coeff;
  double weightIntegral;

  Array2D<double> velCurrent;
  Array1D<double> smoothing;
};

#endif
