/**
 * @file VoxelDataCreation.h
 * @date August, 2024
 */

#ifndef VOXELDATACREATION_H
#define VOXELDATACREATION_H

#include "Config.h"
#include "DataGrid.h"
#include "Export.h"
#include "Grid.h"
#include "Import.h"
#include <sys/stat.h>

class VoxelDataCreation
{
public:
  VoxelDataCreation(Config &conf);
  ~VoxelDataCreation()
  {
  }

  Grid grid;
  SnapShot snap;
  DataGrid data;

  int dim;
  int timeMax;
  int ntInSnapshot;

  std::string outputDir;
  std::string inputDir;

  int nx, ny, nz;              // original grid size
  int nxOpt, nyOpt, nzOpt;     // reference grid size
  int nxData, nyData, nzData;  // data grid size

  double dx, dy, dz;
  double dxOpt, dyOpt, dzOpt;
  double dxData, dyData, dzData;

  double lx, ly, lz;
  double lxOpt, lyOpt, lzOpt;
  double lxData, lyData, lzData;

  double xOrigin, yOrigin, zOrigin;

  Array3D<double> vOrig;
  Array3D<double> vRef;
  Array2D<double> vRefInit;

  void Initialize(Config &conf);
  void createRefAndData();
  void compute(const int step);
  void createData();

  void outputVTK();
  void outputBIN();

private:
  void createDirectories();
  void initializeVectors(Config &conf);
  void initializeNodes(Config &conf);
  void initializeCells(Config &conf);
  void importVelocityData();
  void takeSnapshots();
  void createReferenceData();

  std::tuple<double, double, double> compReferenceGridPoints(int i, int j, int k) const;
  std::tuple<int, int, int> compOriginalGridIndeces(double px, double py, double pz) const;
  std::tuple<double, double, double> calculateInterpolationParameters(double px, double py, double pz, int ix, int iy,
                                                                      int iz) const;

  void validateInterpolationParameters(double s, double t, double u) const;
  int calculateGlobalNodeIndex(int i, int j, int k) const;
  int calculateElementIndex(int ix, int iy, int iz) const;
  void interpolateReferenceData(Array1D<double> &N, int elm, int n, const int step);
  void interpolateInitialReferenceData(Array1D<double> &N, int elm, int n);
};

#endif  // VOXELDATACREATION_H