/**
 * @file VoxelDataCreation.cpp
 * @author
 * @date August, 2024
 */

#include "VoxelDataCreation.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

VoxelDataCreation::VoxelDataCreation(Config &conf)
    : grid(conf), snap(conf), data(conf, grid, snap), dim(conf.dim), timeMax(conf.timeMax), outputDir(conf.outputDir),
      inputDir(conf.inputDir), nx(conf.nx), ny(conf.ny), nz(conf.nz), nxOpt(conf.nxOpt), nyOpt(conf.nyOpt),
      nzOpt(conf.nzOpt), nxData(conf.nxData), nyData(conf.nyData), nzData(conf.nzData), dx(conf.dx), dy(conf.dy),
      dz(conf.dz), dxOpt(conf.dxOpt), dyOpt(conf.dyOpt), dzOpt(conf.dzOpt), dxData(conf.dxData), dyData(conf.dyData),
      dzData(conf.dzData), xOrigin(conf.xOrigin), yOrigin(conf.yOrigin), zOrigin(conf.zOrigin)
{
  createDirectories();
}

void VoxelDataCreation::createDirectories()
{
  std::string output = "output";
  mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  outputDir = "output/" + outputDir;
  mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

  std::string dir;
  dir = outputDir + "/bin";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/data";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/vtk";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}

void VoxelDataCreation::Initialize(Config &conf)
{
  ntInSnapshot = snap.snapInterval * (snap.nSnapShot - 1) + 1;
  initializeVectors(conf);
  initializeCells(conf);
  initializeNodes(conf);
  importVelocityData();
}

void VoxelDataCreation::initializeCells(Config &conf)
{
  grid.cell.nNodesInCell = conf.nNodesInCell;
  grid.cell.nCellsGlobal = conf.nCellsGlobal;
  grid.cell.resize(conf.nCellsGlobal);

  grid.cell.assignNodes(conf);
  grid.cell.assignCoordinates(conf);
  grid.cell.assignCellType(conf);

  if(conf.gridType == GridType::STRUCTURED) {
    grid.cell.nCellsStrGlobal = conf.nx * conf.ny * conf.nz;
  }
}

void VoxelDataCreation::initializeNodes(Config &conf)
{
  grid.node.nNodesGlobal = conf.nNodesGlobal;
  grid.node.x.resize(conf.nNodesGlobal, std::vector<double>(conf.dim));
  grid.node.assignCoordinates(conf);
}

void VoxelDataCreation::initializeVectors(Config &conf)
{
  vOrig.allocate(timeMax, grid.node.nNodesGlobal, 3);
  // uOrig.allocate(timeMax, grid.node.nNodesGlobal);
  // vOrig.allocate(timeMax, grid.node.nNodesGlobal);
  // wOrig.allocate(timeMax, grid.node.nNodesGlobal);
  vRef.allocate(ntInSnapshot, conf.nNodesOptGlobal, 3);
  vRefInit.allocate(conf.nNodesOptGlobal, 3);
  data.snap.vSnap.allocate(data.snap.nSnapShot, conf.nNodesGlobal, 3);
}

void VoxelDataCreation::importVelocityData()
{
  for(int step = 0; step < timeMax; step++) {
    std::string velFile = inputDir + "/velocity_" + std::to_string(step) + ".bin";
    vOrig.importBIN(velFile, step);
  }
}

void VoxelDataCreation::createRefAndData()
{
  takeSnapshots();
  createReferenceData();
  //createData();
}

void VoxelDataCreation::takeSnapshots()
{
  int snapCount = 0;
  for(int step = 0; step < timeMax; step++) {
    if((step >= snap.snapTimeBeginItr) && (snapCount < snap.nSnapShot)) {
      if((step - snap.snapTimeBeginItr) % snap.snapInterval == 0) {
        snap.takeSnapShot(vOrig, grid.node.nNodesGlobal, snapCount, step);
        snapCount++;
      }
    }
  }
}

void VoxelDataCreation::createReferenceData()
{
  for(int step = 0; step < timeMax; step++) {
    if((step >= snap.snapTimeBeginItr) && (step < (snap.snapTimeBeginItr + ntInSnapshot))) {
      compute(step);
    }
  }
}

void VoxelDataCreation::createData()
{
  std::string file_mask = "../../../example/inverse/4dvar/input/Ubend/bend1/data_half/mask.dat";
  data.importMask(file_mask);

  grid.cell.getBoundaries();
  grid.cell.getCenterCoordinates();
  data.setVoxelCenters();
  data.setBoundingBox();

  data.setRefinedGrid();
  data.collectRefinedVoxelId();
  data.collectCFDCellId();

  std::vector<int> cell_inside;

  for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
    if(grid.cell(ic).x[0][0] < data.xOrigin || grid.cell(ic).x[0][0] > data.xOrigin + data.lxData) continue;
    if(grid.cell(ic).x[0][1] < data.yOrigin || grid.cell(ic).x[0][1] > data.yOrigin + data.lyData) continue;
    if(grid.cell(ic).x[0][2] < data.zOrigin || grid.cell(ic).x[0][2] > data.zOrigin + data.lzData) continue;
    cell_inside.push_back(ic);
  }

  for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
    if(data.voxel(iv).mask < 0.5) continue;
    if(mpi.myId == 0) std::cout << "iv" << iv << std::endl;
    for(int ic = 0; ic < cell_inside.size(); ic++) {
      if(data.isCellIncludedInCircle(2 * data.dxData, iv, ic)) {
        data.voxel(iv).cells.push_back(ic);
      }
    }
  }

  //data.collectCellsInCircle(2 * dxData);
  data.collectCFDTimeStep();

  // ----- comp numerical velocity (space) ------
  switch(data.vel_space) {
  case VoxelVelocity::AVERAGE:
    for(int t = 0; t < data.n_cfd_step; t++) {
      if((t >= snap.snapTimeBeginItr) && (t < (snap.snapTimeBeginItr + ntInSnapshot))) {
        for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
          if(data.voxel(iv).mask < 1e-12) continue;
          if(data.voxel(iv).subId != mpi.myId) continue;
          data.average_space(vRef, iv, t);
        }
      }
    }
    break;

  case VoxelVelocity::WEIGHTED_AVERAGE:
    for(int t = 0; t < data.n_cfd_step; t++) {
      if((t >= snap.snapTimeBeginItr) && (t < (snap.snapTimeBeginItr + ntInSnapshot))) {
        int time = t - snap.snapTimeBeginItr;
        if(mpi.myId == 0) std::cout << "t = " << t << std::endl;
        for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
          if(data.voxel(iv).mask < 0.5) continue;
          data.weighted_average_space(vRef, iv, time);
        }
      }
    }
    break;

  case VoxelVelocity::INTERPOLATION:
    for(int t = 0; t < data.n_cfd_step; t++) {
      if((t >= snap.snapTimeBeginItr) && (t < (snap.snapTimeBeginItr + ntInSnapshot))) {
        for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
          if(data.voxel(iv).mask < 1e-12) continue;
          data.linear_interpolate_space(vRef, iv, t);
        }
      }
    }
    break;

  default:
    PetscPrintf(MPI_COMM_WORLD, "undefined VoxelVelocity method");
    exit(1);
  }
  // -----------------------

  if(mpi.myId == 0) std::cout << "v " << std::endl;
  // ----- comp numerical velocity (time) ------
  switch(data.vel_time) {
  case VoxelVelocity::AVERAGE:
    for(int t = 0; t < data.n_mri_step; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        //data.average_time(iv, t);
      }
    }
    break;

  case VoxelVelocity::WEIGHTED_AVERAGE:
    for(int t = 0; t < data.n_mri_step; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        if(data.voxel(iv).mask < 0.5) continue;
        data.weighted_average_time(iv, t);
      }
    }
    break;

  case VoxelVelocity::INTERPOLATION:
    for(int t = 0; t < data.n_mri_step; t++) {
      for(int iv = 0; iv < data.nDataCellsGlobal; iv++) {
        if(data.voxel(iv).mask < 1e-12) continue;
        data.spline_interpolate_time(iv, t);
      }
    }
    break;

  default:
    PetscPrintf(MPI_COMM_WORLD, "undefined VoxelVelocity method");
    exit(1);
  }
}

void VoxelDataCreation::compute(const int step)
{
  bool isInitial = (step == snap.snapTimeBeginItr);
  for(int k = 0; k < nzOpt + 1; k++) {
    for(int j = 0; j < nyOpt + 1; j++) {
      for(int i = 0; i < nxOpt + 1; i++) {
        auto [px, py, pz] = compReferenceGridPoints(i, j, k);
        auto [ix, iy, iz] = compOriginalGridIndeces(px, py, pz);
        auto [s, t, u] = calculateInterpolationParameters(px, py, pz, ix, iy, iz);

        validateInterpolationParameters(s, t, u);
        int n = calculateGlobalNodeIndex(i, j, k);

        int elm = calculateElementIndex(ix, iy, iz);
        Array1D<double> N(grid.cell.nNodesInCell);

        ShapeFunction3D::C3D8_N(N, s, t, u);
        interpolateReferenceData(N, elm, n, step);

        if(isInitial) {
          interpolateInitialReferenceData(N, elm, n);
        }
      }
    }
  }
}

std::tuple<double, double, double> VoxelDataCreation::compReferenceGridPoints(int i, int j, int k) const
{
  return {xOrigin + i * dxOpt, yOrigin + j * dyOpt, zOrigin + k * dzOpt};
}

std::tuple<int, int, int> VoxelDataCreation::compOriginalGridIndeces(double px, double py, double pz) const
{
  int ix = (px / grid.dx) + EPS;
  int iy = (py / grid.dy) + EPS;
  int iz = (pz / grid.dz) + EPS;

  if(ix >= nx) {
    ix--;
  }
  if(iy >= ny) {
    iy--;
  }
  if(iz >= nz) {
    iz--;
  }

  return {ix, iy, iz};
}

std::tuple<double, double, double> VoxelDataCreation::calculateInterpolationParameters(double px, double py, double pz,
                                                                                       int ix, int iy, int iz) const
{
  double s = (px - (ix * grid.dx + 5e-1 * grid.dx)) / (grid.dx / 2.0);
  double t = (py - (iy * grid.dy + 5e-1 * grid.dy)) / (grid.dy / 2.0);
  double u = (pz - (iz * grid.dz + 5e-1 * grid.dz)) / (grid.dz / 2.0);
  return {s, t, u};
}

void VoxelDataCreation::validateInterpolationParameters(double s, double t, double u) const
{
  if(s < -1 - EPS || s > 1 + EPS) {
    throw std::runtime_error("s interpolation error.");
  } else if(t < -1 - EPS || t > 1 + EPS) {
    throw std::runtime_error("t interpolation error.");
  } else if(u < -1 - EPS || u > 1 + EPS) {
    throw std::runtime_error("u interpolation error.");
  }
}

int VoxelDataCreation::calculateGlobalNodeIndex(int i, int j, int k) const
{
  return i + j * (nxOpt + 1) + k * (nxOpt + 1) * (nyOpt + 1);
}

int VoxelDataCreation::calculateElementIndex(int ix, int iy, int iz) const
{
  return ix + iy * grid.nx + iz * grid.nx * grid.ny;
}

void VoxelDataCreation::interpolateReferenceData(Array1D<double> &N, int elm, int n, const int step)
{
  int stepRef = step - snap.snapTimeBeginItr;
  for(int d = 0; d < dim; d++) {
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int node = grid.cell(elm).node[p];
      vRef(stepRef, n, d) += N(p) * vOrig(step, node, d);
    }
  }
}

void VoxelDataCreation::interpolateInitialReferenceData(Array1D<double> &N, int elm, int n)
{
  for(int d = 0; d < dim; d++) {
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int node = grid.cell(elm).node[p];
      vRefInit(n, d) += N(p) * vOrig(snap.snapTimeBeginItr - 1, node, d);
    }
  }
}

void VoxelDataCreation::outputVTK()
{
  for(int step = snap.snapTimeBeginItr; step < snap.snapTimeBeginItr + ntInSnapshot; step++) {
    std::string vtiFile = outputDir + "/vtk/" + "velocityOriginal_" + std::to_string(step) + ".vti";
    EXPORT::exportVectorPointDataVTI<double>(vtiFile, "velocityOriginal", vOrig, nx, ny, nz, dx, dy, dz, step);
  }
  for(int step = 0; step < ntInSnapshot; step++) {
    std::string vtiFile = outputDir + "/vtk/" + "velocityReference_" + std::to_string(step) + ".vti";
    EXPORT::exportVectorPointDataVTI<double>(vtiFile, "velocityReference", vRef, nxOpt, nyOpt, nzOpt, dxOpt, dyOpt,
                                             dzOpt, step);
  }
  // for(int step = 0; step < snap.nSnapShot; step++) {
  //   std::string vtiFile = outputDir + "/vtk/" + "data_" + std::to_string(step) + ".vti";
  //   data.exportVTI(vtiFile, step);
  // }
  std::string vtiFile = outputDir + "/vtk/" + "velocityReference_initial.vti";
  EXPORT::exportVectorPointDataVTI<double>(vtiFile, "velocityReference_initial", vRefInit, nxOpt, nyOpt, nzOpt, dxOpt,
                                           dyOpt, dzOpt);
}

void VoxelDataCreation::outputBIN()
{
  for(int step = 0; step < ntInSnapshot; step++) {
    std::string binFile = outputDir + "/bin/" + "velocityReference_" + std::to_string(step) + ".bin";
    vRef.exportBIN(binFile, step);
  }
  std::string binFile = outputDir + "/bin/" + "velocityReference_initial.bin";
  vRefInit.exportBIN(binFile);

  // for(int step = 0; step < snap.nSnapShot; step++) {
  //   std::string datFile = outputDir + "/data/" + "data_" + std::to_string(step) + ".dat";
  //   data.exportVelCFD_DAT(datFile, step);
  // }
}
