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
#include "Octree.h"
#include <iostream>

struct VoxelInfo
{
  double mask;
  double center[3];
  double minX, minY, minZ;
  double maxX, maxY, maxZ;
  vector<int> cells;
  vector<int> refinedVoxelId;
  Array2D<double> v_cfd, v_mri, v_err;  // (time, dim)
};

struct RefinedVoxelInfo
{
  std::vector<int> node;
  std::vector<std::vector<double>> x;
};

struct BoundingBox {
  double minX, maxX, minY, maxY, minZ, maxZ;
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

  /* data grid */
  int nxData, nyData, nzData;
  double dxData, dyData, dzData;
  double lxData, lyData, lzData;
  double xOrigin, yOrigin, zOrigin;

  int nDataCellsGlobal;
  int nDataNodesGlobal;
  int nDataNodesInCell;

  /* refined voxel for avarage vel computation*/
  int nxRefinedData, nyRefinedData, nzRefinedData;
  double dxRefinedData, dyRefinedData, dzRefinedData;

  int nRefinedCellGlobal;
  int nRefinedNodeGlobal;

  int refinementFactor;

  std::vector<std::vector<double>> refinedCoord;
  std::vector<int> CFDCellId;

  Array3D<VoxelInfo> voxel;
  Array3D<RefinedVoxelInfo> refinedVoxel;

  void initialize(Config &conf);

  void setVoxelCenters();
  void setVoxelBoundaries();
  void collectCellsInVoxel();
  void collectCellsInCircle(const int radious);
  bool isCellCenterIncludedInVoxel(const int iv, const int ic);
  bool isCellBoundaryIncludedInVoxel(const int iv, const int ic);
  bool isCellIncludedInCircle(const double radious, const int iv, const int ic);

  void average(const int iv, const int t);
  void averageTmp(const int iv, const int t);
  void weightedAverage(const int iv, const int t);
  void interpolate();

  void importDAT(const std::string &filename, const int step);
  void exportDAT(const std::string &filename, const int step);
  void exportVTI(const std::string &filename, const int t);
  void exportMaskVTI(const std::string &filename);

  void importMask(const std::string &filename);
  
  void setRefinedGrid();
  void collectRefinedVoxelId();
  void collectCFDCellId();
  bool isRefinedNodeIncludedInCFDCell(const int in, const int ic);

private:
  double coeff;
  double weightIntegral;

  Array2D<double> velCurrent;
  Array1D<double> smoothing;

  void compSmoothing();
  void setBoundingBox();
  bool isRefinedNodeOutside(const int in);

  BoundingBox box;
};

#endif
