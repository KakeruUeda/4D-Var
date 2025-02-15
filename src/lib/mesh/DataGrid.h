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

  int CFDCellId;
  int subId;

  vector<int> cells;
  vector<int> refinedVoxelId;
  Array2D<double> v_cfd, v_mri, v_err;  // (time, dim)
  Array2D<double> v_cfd_refined, v_mri_refined, v_err_refined;
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

  VoxelVelocity vel_space;
  VoxelVelocity vel_time;
  MathTools1D mt1d;
  MathTools3D mt3d;
  Gauss gauss, g1;
  Grid &grid;
  SnapShot &snap;

  /* data grid */
  int nxData, nyData, nzData;
  double dxData, dyData, dzData;
  double lxData, lyData, lzData;
  double xOrigin, yOrigin, zOrigin;

  double dt_mri, dt_cfd;
  double n_mri_step, n_cfd_step;

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

  std::vector<std::vector<int>> cfd_time_steps_id;

  std::vector<std::vector<double>> u_mri_slice;
  std::vector<std::vector<double>> v_mri_slice;
  std::vector<std::vector<double>> w_mri_slice;

  std::vector<std::vector<double>> u_interpolated_slice;
  std::vector<std::vector<double>> v_interpolated_slice;
  std::vector<std::vector<double>> w_interpolated_slice;

  std::vector<double> x_slice;
  std::vector<double> y_slice;
  std::vector<double> z_slice;

  Array3D<VoxelInfo> voxel;
  Array3D<RefinedVoxelInfo> refinedVoxel;

  void initialize(Config &conf);

  void comp_v_mri_refined(const double dt, const int timeMax);
  void comp_v_err_refined();

  void setVoxelCenters();
  void setVoxelBoundaries();
  void collectCellsInVoxel();
  void collectCellsInCircle(double radious);
  bool isCellCenterIncludedInVoxel(const int iv, const int ic);
  bool isCellBoundaryIncludedInVoxel(const int iv, const int ic);
  bool isCellIncludedInCircle(const double radious, const int iv, const int ic);

  void collectCFDTimeStep();

  void extract_slice();

  void average(const int iv, const int t);
  void interpolate(const int iv, const int t);

  void average_space(Array3D<double> &vt, const int iv, const int t);
  void weighted_average_space(Array3D<double> &vt, const int iv, const int t);
  void linear_interpolate_space(Array3D<double> &vt, const int iv, const int t);

  void average_time(const int iv, const int t);
  void weighted_average_time(const int iv, const int t);
  void linear_interpolate_time(const int iv, const int t);
  void spline_interpolate_time(const int iv, const int itm);
  
  void importDAT(const std::string &filename, const int step);
  void exportVelCFD_DAT(const std::string &filename, const int step);
  void exportVelMRI_DAT(const std::string &filename, const int step);
  void exportVelError_DAT(const std::string &filename, const int step);
  void exportVelCFD_t_DAT(const std::string &filename, const int step);
  void exportVelMRI_t_DAT(const std::string &filename, const int step);
  void exportVelError_t_DAT(const std::string &filename, const int step);
  void exportVTI(const std::string &filename, const int t);
  void exportVTI_t(const std::string &filename, const int t);
  void exportMaskVTI(const std::string &filename);

  void importMask(const std::string &filename);
  
  void setBoundingBox();
  void setRefinedGrid();
  void collectRefinedVoxelId();
  void collectCFDCellId();
  void collectVoxelCenterCFDCell();
  bool isRefinedNodeIncludedInCFDCell(const int in, const int ic);
  bool isVoxelCenterIncludedInCFDCell(const int iv, const int ic);

  void gather_voxel_info();

private:
  double coeff;
  double weightIntegral;
  double weight_integral_time;

  std::vector<double> local_voxel_info;
  std::vector<double> global_voxel_info;

  Array2D<double> velCurrent;
  Array1D<double> smoothing;
  Array2D<double> vel_current_time;
  Array1D<double> smoothing_time;

  //void compupte_smoothing(const int iv);
  void comp_smoothing_time(const int t);
  bool isRefinedNodeOutside(const int in);

  BoundingBox box;
};

#endif
