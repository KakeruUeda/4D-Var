/*===========================================
  Do not execute with multiple MPI processes.
=============================================*/
/**
 * @file main.cpp
 * @author K.Ueda
 * @date July, 2024
 */


#include <unistd.h>
#include "DirectProblem.h"
#include "MyMPI.h"
MyMPI mpi;

void createData(Config &conf, Cell &cell, Node &node, DataGrid &voxel, SnapShot &snap);
void createReferenceDA(Config &conf, std::vector<std::vector<double>> &vt,
                       std::vector<std::vector<double>> &velRef);

int main(int argc, char *argv[])
{
  std::string inputFile = argv[1];
  std::string appName = "VOXELDATA";
  Config conf(inputFile, appName);
  if (conf.isReadingError)
    return EXIT_FAILURE;

  std::string dir;
  std::string output = "output";
  mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  conf.outputDir = "output/" + conf.outputDir;
  mkdir(conf.outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = conf.outputDir + "/input_bin";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = conf.outputDir + "/vtk";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

  Cell cell;
  Node node;
  SnapShot snap(conf);
  DataGrid voxel(conf);

  std::vector<std::vector<std::vector<double>>> vOrig;
  std::vector<std::vector<std::vector<double>>> vRef;
  std::vector<std::vector<double>> vRefInit;

  int ntInSnapshot = snap.snapInterval * (snap.nSnapShot - 1) + 1;

  VecTool::resize(vRef, ntInSnapshot, conf.nNodesOptGlobal, conf.dim);
  VecTool::resize(vRefInit, conf.nNodesOptGlobal, conf.dim);
  VecTool::resize(snap.v, snap.nSnapShot, conf.nNodesGlobal, conf.dim);

  cell.resize(conf.nCellsGlobal);
  node.x.resize(conf.nNodesGlobal, std::vector<double>(conf.dim));

  cell.nCellsGlobal = conf.nCellsGlobal;
  cell.nNodesInCell = conf.nNodesInCell;

  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    cell(ic).node.resize(cell.nNodesInCell);
    for (int p = 0; p < cell.nNodesInCell; p++)
    {
      cell(ic).node[p] = conf.cell[ic][p];
    }
  }

  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    cell(ic).x.resize(cell.nNodesInCell, std::vector<double>(conf.dim));
  }

  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell.nNodesInCell; p++)
    {
      for (int d = 0; d < conf.dim; d++)
      {
        cell(ic).x[p][d] = conf.node[cell(ic).node[p]][d];
      }
    }
  }

  for (int in = 0; in < conf.nNodesGlobal; in++)
  {
    for (int d = 0; d < conf.dim; d++)
    {
      node.x[in][d] = conf.node[in][d];
    }
  }

  vOrig.resize(conf.stepMax);
  for (int step = 0; step < conf.stepMax; step++)
  {
    std::string velFile = conf.inputDir + "/velocity_" + to_string(step) + ".bin";
    try
    {
      BIN::importVectorDataBIN(velFile, vOrig[step]);
    }
    catch (const std::runtime_error &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
  }

  int snapCount = 0;
  for (int step = 0; step < conf.stepMax; step++)
  {
    if ((step >= conf.snapTimeBeginItr) && (snapCount < conf.nSnapShot))
    {
      if ((step - conf.snapTimeBeginItr) % conf.snapInterval == 0)
      {
        snap.takeSnapShot(vOrig[step], snapCount, conf.nNodesGlobal, conf.dim);
        snapCount++;
      }
    }
  }

  for (int step = 0; step < conf.stepMax; step++)
  {
    if ((step >= conf.snapTimeBeginItr) && (step < (conf.snapTimeBeginItr + ntInSnapshot)))
    {
      createReferenceDA(conf, vOrig[step], vRef[step - conf.snapTimeBeginItr]);
    }
  }

  for (int step = 0; step < conf.stepMax; step++)
  {
    if (step == (conf.snapTimeBeginItr - 1))
    {
      createReferenceDA(conf, vOrig[step], vRefInit);
    }
  }

  createData(conf, cell, node, voxel, snap);

  // output vtk
  for (int step = 0; step < conf.stepMax; step++)
  {
    std::string vtiFile = conf.outputDir + "/vtk/velocityOriginal_" + to_string(step) + ".vti";
    VTK::exportVectorPointDataVTI(vtiFile, "velOrig", vOrig[step], conf.nx, conf.ny, conf.nz, conf.dx, conf.dy, conf.dz);
  }
  for (int step = 0; step < ntInSnapshot; step++)
  {
    std::string vtiFile = conf.outputDir + "/vtk/velocityReference_" + to_string(step) + ".vti";
    VTK::exportVectorPointDataVTI(vtiFile, "velRef", vRef[step], conf.nxOpt, conf.nyOpt, conf.nzOpt, conf.dxOpt, conf.dyOpt, conf.dzOpt);
  }
  for (int step = 0; step < snap.nSnapShot; step++)
  {
    std::string vtiFile = conf.outputDir + "/vtk/data" + to_string(step) + ".vti";
    VTK::exportVelocityDataVTI(vtiFile, voxel, step);
  }
  std::string vtiFile = conf.outputDir + "/vtk/velocityReference_initial.vti";
  VTK::exportVectorPointDataVTI(vtiFile, "velRefInit", vRefInit, conf.nxOpt, conf.nyOpt, conf.nzOpt, conf.dxOpt, conf.dyOpt, conf.dzOpt);

  // output bin
  for (int step = 0; step < ntInSnapshot; step++)
  {
    std::string binFile = conf.outputDir + "/input_bin/velocityReference_" + to_string(step) + ".bin";
    BIN::exportVectorDataBIN(binFile, vRef[step]);
  }
  std::string binFile = conf.outputDir + "/input_bin/velocityReference_initial.bin";
  BIN::exportVectorDataBIN(binFile, vRefInit);

  std::vector<std::vector<std::vector<double>>> dataTmp;
  VecTool::resize(dataTmp, snap.nSnapShot, voxel.nCellsGlobal, conf.dim);
  for (int step = 0; step < snap.nSnapShot; step++)
  {
    for (int ic = 0; ic < voxel.nCellsGlobal; ic++)
    {
      for (int d = 0; d < conf.dim; d++)
      {
        dataTmp[step][ic][d] = voxel(ic).vCFD[step][d];
      }
    }
  }
  for (int step = 0; step < snap.nSnapShot; step++)
  {
    std::string binFile = conf.outputDir + "/input_bin/data_" + to_string(step) + ".bin";
    BIN::exportVectorDataBIN(binFile, dataTmp[step]);
  }

  std::cout << "Terminated." << std::endl;

  return EXIT_SUCCESS;
}

void createData(Config &conf, Cell &cell, Node &node, DataGrid &voxel, SnapShot &snap)
{
  voxel.range = 5e-1 * sqrt(voxel.dx * voxel.dx + voxel.dy * voxel.dy + voxel.dz * voxel.dz);

  for (int k = 0; k < voxel.nz; k++)
  {
    for (int j = 0; j < voxel.ny; j++)
    {
      for (int i = 0; i < voxel.nx; i++)
      {
        int ic = k * voxel.nx * voxel.ny + j * voxel.nx + i;
        voxel(ic).center[0] = conf.xOrigin + (5e-1 + i) * voxel.dx;
        voxel(ic).center[1] = conf.yOrigin + (5e-1 + j) * voxel.dy;
        voxel(ic).center[2] = conf.zOrigin + (5e-1 + k) * voxel.dz;
        voxel(ic).setNearCell(node, cell, voxel.range, conf.dim);
      }
    }
  }

  for (int t = 0; t < snap.nSnapShot; t++)
  {
    for (int ic = 0; ic < voxel.nCellsGlobal; ic++)
    {
      voxel(ic).average(cell, snap.v[t], t, conf.dim);
    }
  }
}

void createReferenceDA(Config &conf, std::vector<std::vector<double>> &vOrig,
                       std::vector<std::vector<double>> &vRef)
{
  double px, py, pz;
  int ix, iy, iz;
  double s, t, u;

  for (int k = 0; k < conf.nzOpt + 1; k++)
  {
    for (int j = 0; j < conf.nyOpt + 1; j++)
    {
      for (int i = 0; i < conf.nxOpt + 1; i++)
      {
        px = conf.xOrigin + i * conf.dxOpt;
        py = conf.yOrigin + j * conf.dyOpt;
        pz = conf.zOrigin + k * conf.dzOpt;

        ix = (px / conf.dx) + EPS;
        iy = (py / conf.dy) + EPS;
        iz = (pz / conf.dz) + EPS;

        if (ix >= conf.nx)
        {
          ix = ix - 1;
        }
        if (iy >= conf.ny)
        {
          iy = iy - 1;
        }
        if (iz >= conf.nz)
        {
          iz = iz - 1;
        }

        s = px - (ix * conf.dx + 5e-1 * conf.dx);
        t = py - (iy * conf.dy + 5e-1 * conf.dy);
        u = pz - (iz * conf.dz + 5e-1 * conf.dz);

        s = s / (conf.dx / 2e0);
        t = t / (conf.dy / 2e0);
        u = u / (conf.dz / 2e0);

        if (s < -1 - EPS || s > 1 + EPS)
        {
          std::cout << "s interpolation error." << std::endl;
        }
        else if (t < -1 - EPS || t > 1 + EPS)
        {
          std::cout << "t interpolation error." << std::endl;
        }
        else if (u < -1 - EPS || u > 1 + EPS)
        {
          std::cout << "u interpolation error." << std::endl;
        }

        int n = i + j * (conf.nxOpt + 1) + k * (conf.nxOpt + 1) * (conf.nyOpt + 1);
        int elm = ix + iy * conf.nx + iz * conf.nx * conf.ny;

        std::vector<double> N;
        VecTool::resize(N, conf.nNodesInCell);

        ShapeFunction3D::C3D8_N(N, s, t, u);

        for (int d = 0; d < conf.dim; d++)
        {
          for (int p = 0; p < conf.nNodesInCell; p++)
          {
            vRef[n][d] += N[p] * vOrig[conf.cell[elm][p]][d];
          }
        }
      }
    }
  }
}
