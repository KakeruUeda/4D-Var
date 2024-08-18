/**
 * @file Config.cpp
 * @author K.Ueda
 * @date May, 2024
 */
#ifndef CONFIG_H
#define CONFIG_H

#include "Array.h"
#include "Define.h"
#include "Import.h"
#include "MyMPI.h"
#include "TextParser.h"
#include "Tool.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

extern MyMPI mpi;

enum class Application
{
  SNS = 0,
  USNS = 1,
  VOXELDATACREATION = 2,
  FDVAR = 3,
  FLOWRATE = 4,
  MAE = 5,
  GRIDCREATION = 6
};

enum class GridType
{
  STRUCTURED = 0,
  UNSTRUCTURED = 1
};

enum class ControlBoundaryFace
{
  left = 0,
  right = 1,
  top = 2,
  bottom = 3,
  front = 4,
  back = 5
};


enum class VoxelVelocity
{
  POINTSPREAD = 0,
  AVERAGE = 1,
  INTERPOLATION = 2
};

enum class CrossSection
{
  XY = 0,
  YZ = 1,
  ZX = 2
};

class Config
{
public:
  Config(std::string inputFile, std::string appName);

  TextParser tp;
  Application app;
  GridType gridType;
  ControlBoundaryFace controlBoundary;
  ControlBoundaryFace inletCB;
  ControlBoundaryFace outletCB;
  VoxelVelocity vvox;
  CrossSection crossSection;

  // Basic parameter
  int dim, nOMP;
  int outputItr;
  std::string outputDir;

  // Physical parameter
  double rho, mu, Re;

  double NRtolerance;

  // CostFunction parameter
  double aCF, bCF, gCF;
  int loopMax;

  double alphaX0, alphaX;

  // Time parameter
  double dt;
  int timeMax;
  int pulsatileFlow;
  int pulseBeginItr;
  double T;

  // Darcy Parameter
  double alpha, resistance;

  // Grid parameter
  bool fluidExtraction;
  int extractFluid;
  int nx, ny, nz;
  double lx, ly, lz;
  double dx, dy, dz;
  int nxNodes, nyNodes, nzNodes;
  int nxCells, nyCells, nzCells;
  int nCellsGlobal, nNodesGlobal, nNodesInCell;
  int nStrCellsGlobal, nStrNodesGlobal;

  std::vector<int> cellId;
  std::vector<int> nodeId;

  // SnapShot parameter
  int nSnapShot;
  int snapInterval;
  int snapTimeBeginItr;

  // DataGrid parameter
  int nData;
  double xOrigin, yOrigin, zOrigin;
  int nxData, nyData, nzData;
  double lxData, lyData, lzData;
  double dxData, dyData, dzData;
  int nNodesInCellData;
  int nDataCellsGlobal;

  // Voxel creation
  int nxOpt, nyOpt, nzOpt;
  double lxOpt, lyOpt, lzOpt;
  double dxOpt, dyOpt, dzOpt;
  int nCellsOptGlobal;
  int nNodesOptGlobal;
  std::string inputDir;

  // For image data
  std::vector<double> phi;

  std::map<int, std::vector<double>> vDirichlet;
  std::map<int, std::vector<double>> vDirichletWall;
  std::map<int, double> pDirichlet;

  std::vector<int> controlBoundaryMap;
  std::vector<int> planeDir;

  std::vector<bool> isBoundaryEdge;

  double center[3];
  double R, Q, maxVelocity;

  std::vector<int> CBNodeMap;
  std::vector<int> CBCellMap;
  std::vector<std::vector<int>> CBNodeMapInCell;
  int CBExtraction;

  // For cell and node data
  std::vector<std::vector<double>> node;
  std::vector<std::vector<int>> cell;

  std::vector<int> sortCell;
  std::vector<int> sortNode;

  // Error parameter
  bool isReadingError = false;

  // data parameter
  int nControlNodesInCell;
  std::string dataDir;
  
  // post inverse parameter
  int nRef;
  int crossPoint;
  int flowRateVelDir;
  std::vector<std::vector<std::vector<double>>> velRef;
  std::vector<std::vector<std::vector<double>>> velOpt;

private:
  void setApplication(std::string appName);
  void tryOpenConfigFile(std::string inputFile);
  void tryReadConfigFile();
  void readConfigFile();

public:
  void setBoundaryVelocityValue(std::string face, double value[3]);
  void setBoundaryPressureValue(std::string face, const double value);
  void setBoundaryPoiseuilleValue(std::string face);
  void setControlBoundary();
  void setStrGrid();
  double setStrCoordinate(const int i, const int j, const int k, const int d);
  int setStrNode(const int i, const int j, const int k, const int p);

public:
  void setSolidDirichletValue();
  void filterFluidGrid();

  std::set<int> fluidUniqueCells;
  std::set<int> fluidUniqueNodes;
  std::set<int> fluidUniqueCBCells;
  std::set<int> fluidUniqueCBNodes;
  std::set<int> fluidUniqueCBCellsIdx;

  std::vector<int> vecFluidUniqueNodes;

  std::unordered_map<int, int> cellMapping;
  std::unordered_map<int, int> nodeMapping;

  void getUniqueCells();
  void getUniqueNodes();
  void getUniqueCBCells();
  void getUniqueCBNodes();
  void getUniqueCBCellsIdx();
  void getNewFilterdMap();
  void filterCell();
  void filterPhi();
  void filterNode();
  void filterMapCB();
  void filterMapCBCell();
  void filterMapCBInCell();
  void filterVelocityDirichlet();
  void filterPressureDirichlet();
  void applyMapping();
};


class TextReaderInterface
{
protected:
  std::shared_ptr<Config> ptr;

public:
  TextReaderInterface()
  {
  }
  virtual ~TextReaderInterface() = default;

  virtual void readBasicInfo(Config &conf);
  virtual void readGridInfo(Config &conf);
  virtual void readBoundaryInfo(Config &conf);
  virtual void readPhysicalInfo(Config &conf);
  virtual void readDarcyInfo(Config &conf);
  virtual void readTimeInfo(Config &conf);
};

class TextReaderGridCreation : public TextReaderInterface
{
public:
  TextReaderGridCreation()
  {
  }
  void readStructuredBoundaryInfo(Config &conf);
  void readBasicInfo(Config &conf) override;
  void readGridInfo(Config &conf) override;
};

class TextReaderUSNS : public TextReaderInterface
{
public:
  TextReaderUSNS()
  {
  }
  void readBasicInfo(Config &conf) override;
  void readGridInfo(Config &conf) override;
  void readBoundaryInfo(Config &conf) override;
  void readPhysicalInfo(Config &conf) override;
  void readDarcyInfo(Config &conf) override;
  void readTimeInfo(Config &conf) override;
};

class TextReaderVoxelDataCreation : public TextReaderInterface
{
public:
  TextReaderVoxelDataCreation()
  {
  }
  void readSnapInfo(Config &conf);
  void readOriginalInfo(Config &conf);
  void readBasicInfo(Config &conf) override;
  void readGridInfo(Config &conf) override;
};

class TextReader4DVar : public TextReaderInterface
{
public:
  TextReader4DVar()
  {
  }
  void readBasicInfo(Config &conf) override;
  void readGridInfo(Config &conf) override;
  void readBoundaryInfo(Config &conf) override;
  void readPhysicalInfo(Config &conf) override;
  void readDarcyInfo(Config &conf) override;
  void readTimeInfo(Config &conf) override;
  void readInverseInfo(Config &conf);
  void readDataInfo(Config &conf);
};

#endif