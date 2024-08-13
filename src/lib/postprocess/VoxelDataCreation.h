/**
 * @file VoxelDataCreation.h
 * @date August, 2024
 */

#ifndef VOXELDATACREATION_H
#define VOXELDATACREATION_H

#include <sys/stat.h>
#include "Config.h"
#include "Grid.h"
#include "Export.h"
#include "Import.h"

class VoxelDataCreation
{
public:
  VoxelDataCreation(Config &conf);
  ~VoxelDataCreation() {}

  void Initialize(Config &conf);
  void createRefAndData();
  void createReferenceDA(std::vector<std::vector<double>> &vOrig,
                         std::vector<std::vector<double>> &vRef);
  void createData();

  void outputVTK();
  void outputBIN();

private:
  void createDirectories();
  void initializeVectors(const Config &conf);
  void initializeGrid(const Config &conf);
  void importVelocityData();
  void takeSnapshots();
  void createReferenceData();
  void createInitialReferenceData();
  void initializeVoxelCenters();
  void computeVoxelAverages();

  std::tuple<double, double, double> calculateVoxelCenter(int i, int j, int k) const;
  std::tuple<int, int, int> calculateGridIndices(double px, double py, double pz) const;
  std::tuple<double, double, double> calculateInterpolationParameters(double px, double py, double pz, int ix, int iy, int iz) const;

  void validateInterpolationParameters(double s, double t, double u) const;
  int calculateGlobalNodeIndex(int i, int j, int k) const;
  int calculateElementIndex(int ix, int iy, int iz) const;
  void interpolateReferenceData(const std::vector<std::vector<double>> &vOrig,
                                std::vector<std::vector<double>> &vRef,
                                const std::vector<double> &N, int elm, int n);

  void exportVelocityVTI(const std::string &prefix, std::vector<std::vector<double>> &data, int step,
                         int nx, int ny, int nz, double dx, double dy, double dz);
  void exportVoxelDataVTI(const std::string &prefix, int step);
  void exportVelocityBIN(const std::string &prefix, std::vector<std::vector<double>> &data, int step);
  void exportVoxelDataBIN(const std::string &prefix, std::vector<std::vector<double>> &data, int step);

  Grid grid;
  SnapShot snap;
  DataGrid voxel;

  int dim;
  int stepMax;
  int ntInSnapshot;

  std::string outputDir;
  std::string inputDir;

  int nx, ny, nz;             // original grid size
  int nxOpt, nyOpt, nzOpt;    // reference grid size
  int nxData, nyData, nzData; // data grid size

  double dx, dy, dz;
  double dxOpt, dyOpt, dzOpt;
  double dxData, dyData, dzData;

  double lx, ly, lz;
  double lxOpt, lyOpt, lzOpt;
  double lxData, lyData, lzData;

  double xOrigin, yOrigin, zOrigin;

  std::vector<std::vector<std::vector<double>>> vOrig;
  std::vector<std::vector<std::vector<double>>> vRef;
  std::vector<std::vector<double>> vRefInit;
};

#endif // VOXELDATACREATION_H