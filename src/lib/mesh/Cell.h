/**
 * @file Cell.h
 * @author K.Ueda
 * @date Jun, 2024
 */

#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <cassert>
#include <vector>
#include "Array.h"
#include "VTKCellType.h"
#include "Config.h"

struct CellInfo
{
public:
  VTKCellType cellType;
  int subId;
  double phi;

  std::vector<int> node, nodeNew;
  std::vector<int> dofsMap, dofsBCsMap;
  std::vector<int> dofStart;
  std::vector<int> dofStartPlane;
  std::vector<std::vector<double>> x;

  std::vector<int> dofsMapWall, dofsBCsMapWall;
};

class Cell
{
public:
  Cell();
  Cell(Config &conf);
  virtual ~Cell();

  inline CellInfo &operator()(int n)
  {
    return data[n];
  }

  inline int size() const
  {
    return data.size();
  }

  inline void resize(int n)
  {
    data.resize(n);
  }

  void initialize(Config &conf);
  void initializeAdjoint(Config &conf);

  void initializeDataStructures(Config &conf);
  void assignNodes(Config &conf);
  void assignCoordinates(Config &conf);
  void assignPhi(Config &conf);
  void assignSubId(Config &conf);
  void assignCellType(Config &conf);

  int nCellsGlobal;
  int nCellsLocal;
  int nNodesInCell;
  int nCellsStrGlobal;

private:
  std::vector<CellInfo> data;
};

#endif // CELL_H