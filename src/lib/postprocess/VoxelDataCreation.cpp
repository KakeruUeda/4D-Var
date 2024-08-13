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
    : snap(conf), voxel(conf), dim(conf.dim), stepMax(conf.stepMax), outputDir(conf.outputDir),
      inputDir(conf.inputDir), nx(conf.nx), ny(conf.ny), nz(conf.nz), nxOpt(conf.nxOpt),
      nyOpt(conf.nyOpt), nzOpt(conf.nzOpt), nxData(conf.nxData), nyData(conf.nyData),
      nzData(conf.nzData), dx(conf.dx), dy(conf.dy), dz(conf.dz), dxOpt(conf.dxOpt),
      dyOpt(conf.dyOpt), dzOpt(conf.dzOpt), dxData(conf.dxData), dyData(conf.dyData),
      dzData(conf.dzData), xOrigin(conf.xOrigin), yOrigin(conf.yOrigin), zOrigin(conf.zOrigin)
{
  createDirectories();
}

void VoxelDataCreation::createDirectories()
{
  mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  outputDir = outputDir + "/" + outputDir;
  mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  mkdir((outputDir + "/bin").c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  mkdir((outputDir + "/vtk").c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}

void VoxelDataCreation::Initialize(Config &conf)
{
  ntInSnapshot = snap.snapInterval * (snap.nSnapShot - 1) + 1;

  initializeVectors(conf);
  initializeGrid(conf);
  importVelocityData();
}

void VoxelDataCreation::initializeVectors(const Config &conf)
{
  VecTool::resize(vRef, ntInSnapshot, conf.nNodesOptGlobal, conf.dim);
  VecTool::resize(vRefInit, conf.nNodesOptGlobal, conf.dim);
  VecTool::resize(snap.v, snap.nSnapShot, conf.nNodesGlobal, conf.dim);
  vOrig.resize(stepMax);
}

void VoxelDataCreation::initializeGrid(const Config &conf)
{
  grid.cell.resize(conf.nCellsGlobal);
  grid.node.x.resize(conf.nNodesGlobal, std::vector<double>(conf.dim));
  grid.cell.nCellsGlobal = conf.nCellsGlobal;
  grid.cell.nNodesInCell = conf.nNodesInCell;

  for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
    grid.cell(ic).node.resize(grid.cell.nNodesInCell);
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      grid.cell(ic).node[p] = conf.cell[ic][p];
    }
  }

  for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
    grid.cell(ic).x.resize(grid.cell.nNodesInCell, std::vector<double>(dim));
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      for(int d = 0; d < dim; d++) {
        grid.cell(ic).x[p][d] = conf.node[grid.cell(ic).node[p]][d];
      }
    }
  }

  for(int in = 0; in < conf.nNodesGlobal; in++) {
    for(int d = 0; d < dim; d++) {
      grid.node.x[in][d] = conf.node[in][d];
    }
  }
}

void VoxelDataCreation::importVelocityData()
{
  for(int step = 0; step < stepMax; step++) {
    std::string velFile = inputDir + "/velocity_" + std::to_string(step) + ".bin";
    IMPORT::importVectorDataBIN<double>(velFile, vOrig[step]);
  }
}

void VoxelDataCreation::createRefAndData()
{
  takeSnapshots();
  createReferenceData();
  createInitialReferenceData();
  createData();
}

void VoxelDataCreation::takeSnapshots()
{
  int snapCount = 0;
  for(int step = 0; step < stepMax; step++) {
    if((step >= snap.snapTimeBeginItr) && (snapCount < snap.nSnapShot)) {
      if((step - snap.snapTimeBeginItr) % snap.snapInterval == 0) {
        snap.takeSnapShot(vOrig[step], snapCount, grid.node.nNodesGlobal, dim);
        snapCount++;
      }
    }
  }
}

void VoxelDataCreation::createReferenceData()
{
  for(int step = 0; step < stepMax; step++) {
    if((step >= snap.snapTimeBeginItr) && (step < (snap.snapTimeBeginItr + ntInSnapshot))) {
      createReferenceDA(vOrig[step], vRef[step - snap.snapTimeBeginItr]);
    }
  }
}

void VoxelDataCreation::createInitialReferenceData()
{
  for(int step = 0; step < stepMax; step++) {
    if(step == (snap.snapTimeBeginItr - 1)) {
      createReferenceDA(vOrig[step], vRefInit);
    }
  }
}

void VoxelDataCreation::createData()
{
  voxel.range = 5e-1 * sqrt(voxel.dx * voxel.dx + voxel.dy * voxel.dy + voxel.dz * voxel.dz);
  initializeVoxelCenters();
  computeVoxelAverages();
}

void VoxelDataCreation::initializeVoxelCenters()
{
  for(int k = 0; k < voxel.nz; k++) {
    for(int j = 0; j < voxel.ny; j++) {
      for(int i = 0; i < voxel.nx; i++) {
        int ic = k * voxel.nx * voxel.ny + j * voxel.nx + i;
        voxel(ic).center[0] = xOrigin + (5e-1 + i) * voxel.dx;
        voxel(ic).center[1] = yOrigin + (5e-1 + j) * voxel.dy;
        voxel(ic).center[2] = zOrigin + (5e-1 + k) * voxel.dz;
        voxel(ic).setNearCell(grid.node, grid.cell, voxel.range, dim);
      }
    }
  }
}

void VoxelDataCreation::computeVoxelAverages()
{
  for(int t = 0; t < snap.nSnapShot; t++) {
    for(int ic = 0; ic < voxel.nCellsGlobal; ic++) {
      voxel(ic).average(grid.cell, snap.v[t], t, dim);
    }
  }
}

void VoxelDataCreation::createReferenceDA(std::vector<std::vector<double>> &vOrig,
                                          std::vector<std::vector<double>> &vRef)
{
  for(int k = 0; k < nzOpt + 1; k++) {
    for(int j = 0; j < nyOpt + 1; j++) {
      for(int i = 0; i < nxOpt + 1; i++) {
        auto [px, py, pz] = calculateVoxelCenter(i, j, k);
        auto [ix, iy, iz] = calculateGridIndices(px, py, pz);
        auto [s, t, u] = calculateInterpolationParameters(px, py, pz, ix, iy, iz);

        validateInterpolationParameters(s, t, u);

        int n = calculateGlobalNodeIndex(i, j, k);
        int elm = calculateElementIndex(ix, iy, iz);

        std::vector<double> N;
        VecTool::resize(N, grid.cell.nNodesInCell);

        ShapeFunction3D::C3D8_N(N, s, t, u);
        interpolateReferenceData(vOrig, vRef, N, elm, n);
      }
    }
  }
}

std::tuple<double, double, double> VoxelDataCreation::calculateVoxelCenter(int i, int j,
                                                                           int k) const
{
  return {xOrigin + i * dxOpt, yOrigin + j * dyOpt, zOrigin + k * dzOpt};
}

std::tuple<int, int, int> VoxelDataCreation::calculateGridIndices(double px, double py,
                                                                  double pz) const
{
  int ix = (px / dx) + EPS;
  int iy = (py / dy) + EPS;
  int iz = (pz / dz) + EPS;

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

std::tuple<double, double, double>
VoxelDataCreation::calculateInterpolationParameters(double px, double py, double pz, int ix, int iy,
                                                    int iz) const
{
  double s = (px - (ix * dx + 5e-1 * dx)) / (dx / 2.0);
  double t = (py - (iy * dy + 5e-1 * dy)) / (dy / 2.0);
  double u = (pz - (iz * dz + 5e-1 * dz)) / (dz / 2.0);
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
  return ix + iy * nx + iz * nx * ny;
}

void VoxelDataCreation::interpolateReferenceData(const std::vector<std::vector<double>> &vOrig,
                                                 std::vector<std::vector<double>> &vRef,
                                                 const std::vector<double> &N, int elm, int n)
{
  for(int d = 0; d < dim; d++) {
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int node = grid.cell(elm).node[p];
      vRef[node][d] += N[p] * vOrig[node][d];
    }
  }
}

void VoxelDataCreation::outputVTK()
{
  for(int step = 0; step < stepMax; step++) {
    exportVelocityVTI("velocityOriginal", vOrig[step], step, nx, ny, nz, dx, dy, dz);
  }
  for(int step = 0; step < ntInSnapshot; step++) {
    exportVelocityVTI("velocityReference", vRef[step], step, nxOpt, nyOpt, nzOpt, dxOpt, dyOpt,
                      dzOpt);
  }
  for(int step = 0; step < snap.nSnapShot; step++) {
    exportVoxelDataVTI("data", step);
  }
  exportVelocityVTI("velocityReference_initial", vRefInit, -1, nxOpt, nyOpt, nzOpt, dxOpt, dyOpt,
                    dzOpt);
}

void VoxelDataCreation::exportVelocityVTI(const std::string &prefix,
                                          std::vector<std::vector<double>> &data, int step, int nx,
                                          int ny, int nz, double dx, double dy, double dz)
{
  std::string vtiFile = outputDir + "/vtk/" + prefix + "_" + std::to_string(step) + ".vti";
  EXPORT::exportVectorPointDataVTI<double>(vtiFile, prefix.c_str(), data, nx, ny, nz, dx, dy, dz);
}

void VoxelDataCreation::exportVoxelDataVTI(const std::string &prefix, int step)
{
  std::string vtiFile = outputDir + "/vtk/" + prefix + "_" + std::to_string(step) + ".vti";
  EXPORT::exportVelocityDataVTI(vtiFile, voxel, step);
}

void VoxelDataCreation::outputBIN()
{
  for(int step = 0; step < ntInSnapshot; step++) {
    exportVelocityBIN("velocityReference", vRef[step], step);
  }
  exportVelocityBIN("velocityReference_initial", vRefInit, -1);

  std::vector<std::vector<std::vector<double>>> dataTmp;
  VecTool::resize(dataTmp, snap.nSnapShot, voxel.nCellsGlobal, dim);

  for(int step = 0; step < snap.nSnapShot; step++) {
    for(int ic = 0; ic < voxel.nCellsGlobal; ic++) {
      for(int d = 0; d < dim; d++) {
        dataTmp[step][ic][d] = voxel(ic).vCFD[step][d];
      }
    }
  }

  for(int step = 0; step < snap.nSnapShot; step++) {
    exportVoxelDataBIN("data", dataTmp[step], step);
  }
}

void VoxelDataCreation::exportVelocityBIN(const std::string &prefix,
                                          std::vector<std::vector<double>> &data, int step)
{
  std::string binFile = outputDir + "/bin/" + prefix + "_" + std::to_string(step) + ".bin";
  EXPORT::exportVectorDataBIN<double>(binFile, data);
}

void VoxelDataCreation::exportVoxelDataBIN(const std::string &prefix,
                                           std::vector<std::vector<double>> &data, int step)
{
  std::string binFile = outputDir + "/bin/" + prefix + "_" + std::to_string(step) + ".bin";
  EXPORT::exportVectorDataBIN<double>(binFile, data);
}