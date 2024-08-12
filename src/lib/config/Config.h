/**
 * @file Config.cpp
 * @author K.Ueda
 * @date May, 2024
 */
#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <string>
#include <fstream>
#include <sstream>
#include "TextParser.h"
#include "Array.h"
#include "MyMPI.h"
#include "Define.h"
#include "Tool.h"
#include "Import.h"

extern MyMPI mpi;

enum class Application
{
  SNS = 0,
  USNS = 1,
  VOXELDATA = 2,
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

enum class ControlBoundary
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
  ControlBoundary controlBoundary;
  ControlBoundary inletCB;
  ControlBoundary outletCB;
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
  int extractFluid;
  int nx, ny, nz;
  double lx, ly, lz;
  double dx, dy, dz;
  int nxNodes, nyNodes, nzNodes;
  int nxCells, nyCells, nzCells;
  int nCellsGlobal, nNodesGlobal, nNodesInCell;

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
  int nCellsDataGlobal;

  // Voxel creation
  int nxOpt, nyOpt, nzOpt;
  double lxOpt, lyOpt, lzOpt;
  double dxOpt, dyOpt, dzOpt;
  int nCellsOptGlobal;
  int nNodesOptGlobal;
  int stepMax;
  std::string inputDir;

  // Boundary parameter for stgrid
  std::vector<std::string> bdStr;
  std::vector<std::string> bdType;
  std::vector<std::vector<double>> bdValue;

  // For image data
  std::vector<double> phi;

  // For Dirichlet boundary data
  std::vector<int> vDirichletNode;
  std::vector<int> pDirichletNode;
  std::vector<std::vector<double>> vDirichletValue;
  std::vector<double> pDirichletValue;

  std::map<int, std::vector<double>> vDirichlet;
  std::map<int, std::vector<double>> vDirichletWall;
  std::map<int, double> pDirichlet;

  std::vector<int> controlBoundaryMap;
  std::vector<int> planeDir;

  std::vector<bool> isBoundaryEdge;

  double center[3];
  double R, Q, maxVelocity;

  std::vector<int> mapCB;
  std::vector<int> mapCBCell;
  std::vector<std::vector<int>> mapCBInCell;
  int extractCB;

  // For cell and node data
  std::vector<std::vector<double>> node;
  std::vector<std::vector<int>> cell;

  std::vector<int> sortCell;
  std::vector<int> sortNode;

  // Error parameter
  bool isReadingError = false;

  // data parameter
  int nControlNodesInCell;
  std::vector<std::vector<std::vector<double>>> velocityData;
  std::vector<std::vector<int>> controlNodeInCell;
  std::vector<int> controlCellMap;

  // post inverse parameter
  int nRef;
  int crossPoint;
  int flowRateVelDir;
  std::vector<std::vector<std::vector<double>>> velRef;
  std::vector<std::vector<std::vector<double>>> velOpt;

  void setSolidBoundary();
  void setFluidDomain();

private:
  void setApplication(std::string appName);
  void tryOpenConfigFile(std::string inputFile);
  void tryReadConfigFile();
  void readConfigFile();

  void readGridType();
  void readStrGridParameter();
  void readGridParameter();
  void readSubGridParameter();
  void readStrBoundaryParameter();
  void readBoundaryParameter();
  void readControlBoundaryParameter();
  void readBasicParameter();
  void readPysicalParameter();
  void readNRParameter();
  void readDarcyParameter();
  void readTimeParameter();
  void readInverseParameter();
  void readDataParameter();
  void readPostInverseBasicParameter();
  void readPostInverseVelocityParameter();
  void readPostInverseFlowRateParameter();
  void readVoxelCreationParameter();
  void readStrBoundaryValue(std::string face, std::string labelFace, std::string labelType, std::string labelValue);
  void setBoundaryVelocityValue(std::string face, double value[3]);
  void setBoundaryPressureValue(std::string face, const double value);
  void setBoundaryPoiseuilleValue(std::string face);
  void setControlBoundary();
  void setStrGrid();
  double setStrCoordinate(const int i, const int j, const int k, const int d);
  int setStrNode(const int i, const int j, const int k, const int p);
  

public:
  // Additional Setting for Structured Grid
  void setSolidDirichletValue();
  void filterFluidGrid();
  
  std::set<int> uniqueCells;
  std::set<int> uniqueNodes;
  std::set<int> uniqueCBCells;
  std::set<int> uniqueCBNodes;
  std::set<int> uniqueCBCellsIdx;

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

#endif