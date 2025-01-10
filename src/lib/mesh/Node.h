/**
 * @file Node.h
 * @author K.Ueda
 * @date Jun, 2024
 */

#ifndef NODE_H
#define NODE_H

#include "Array.h"
#include "Config.h"
#include <cassert>
#include <iostream>
#include <vector>

class Node
{
public:
  Node() {};
  Node(Config &conf) : nNodesGlobal(conf.nNodesGlobal), nNodesLocal(0)
  {
    if(conf.gridType == GridType::STRUCTURED) {
      nStrNodesGlobal = conf.nStrNodesGlobal;
    }
  }
  virtual ~Node()
  {
  }

  int nNodesGlobal, nNodesLocal;
  int nStrNodesGlobal;

  std::vector<int> sortNode;
  std::vector<int> map, mapNew;
  std::vector<int> subId;
  std::vector<int> nDofsOnNode, nDofsOnNodeNew;
  std::vector<std::vector<int>> dofsMap, dofsBCsMap;
  std::vector<std::vector<int>> dofsMapNew, dofsBCsMapNew;
  std::vector<std::vector<bool>> isDirichlet, isDirichletNew;
  std::vector<std::vector<double>> type;
  std::vector<std::vector<double>> x;

  /// add ///
  std::vector<std::vector<int>> dofsMapWall, dofsBCsMapWall;
  std::vector<std::vector<int>> dofsMapWallNew, dofsBCsMapWallNew;
  std::vector<std::vector<bool>> isDirichletWall, isDirichletWallNew;
  ///////////

  // Main variable
  std::vector<std::vector<double>> v0, v, vPrev;
  std::vector<double> p;
  std::vector<std::vector<std::vector<double>>> vt;
  std::vector<std::vector<double>> pt;

  // Adjoint variable
  std::vector<std::vector<double>> w;
  std::vector<std::vector<double>> wPrev;
  std::vector<double> q;
  std::vector<double> qPrev;
  std::vector<std::vector<double>> l;
  std::vector<std::vector<std::vector<double>>> lt;
  std::vector<std::vector<std::vector<double>>> wt;
  std::vector<std::vector<double>> qt;

  // for vti visualization
  int nNodesStrGlobal;
  std::vector<std::vector<double>> vvti;
  std::vector<double> pvti;
  std::vector<std::vector<double>> wvti;
  std::vector<std::vector<double>> lvti;
  std::vector<double> qvti;

  int getDof(const int in, const int id)
  {
    return dofsMapNew[in][id];
  }
  void assignCoordinates(Config &conf);
  void initialize(Config &conf);
  void initializeAdjoint(Config &conf, std::vector<int> &controlBoundaryMap);
  void newMapping();
};

class SnapShot
{
public:
  SnapShot()
  {
  }
  SnapShot(Config &conf)
      : nSnapShot(conf.nSnapShot), snapInterval(conf.snapInterval), snapTimeBeginItr(conf.snapTimeBeginItr)
  {
    timeMax = (nSnapShot-1) * snapInterval + 1;
  }

  int nSnapShot;
  int snapInterval;
  int snapTimeBeginItr;
  int timeMax;

  std::vector<std::vector<std::vector<double>>> v;
  Array3D<double> vSnap;
  void takeSnapShot(Array2D<double> &v, const int nNodesGlobal, const int snapCount);
  void takeSnapShot(Array3D<double> &v, const int nNodesGlobal, const int snapCount, const int step);
};

#endif