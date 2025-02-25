/**
 * @file Node.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "Node.h"

void Node::initialize(Config &conf)
{
  // initialize variables
  x.resize(nNodesGlobal, std::vector<double>(conf.dim));
  nDofsOnNode.resize(nNodesGlobal);
  isDirichlet.resize(nNodesGlobal);
  subId.resize(nNodesGlobal);
  map.resize(nNodesGlobal);
  dofsMap.resize(nNodesGlobal);
  dofsBCsMap.resize(nNodesGlobal);
  sortNode.resize(nNodesGlobal);

  //// add ////
  isDirichletWall.resize(nNodesGlobal);
  dofsMapWall.resize(nNodesGlobal);
  dofsBCsMapWall.resize(nNodesGlobal);
  /////////////

  for(int in = 0; in < nNodesGlobal; in++) {
    subId[in] = conf.nodeId[in];
  }

  nNodesLocal = count(conf.nodeId.begin(), conf.nodeId.end(), mpi.myId);

  for(int in = 0; in < nNodesGlobal; in++) {
    nDofsOnNode[in] = conf.dim + 1;
    isDirichlet[in].resize(nDofsOnNode[in], false);
  }

  int tmp = 0;

  for(int in = 0; in < nNodesGlobal; in++) {
    dofsMap[in].resize(nDofsOnNode[in]);
    dofsBCsMap[in].resize(nDofsOnNode[in]);
    for(int id = 0; id < nDofsOnNode[in]; id++) {
      dofsMap[in][id] = tmp;
      dofsBCsMap[in][id] = tmp;
      tmp++;
    }
  }

  /// add ///
  for(int in = 0; in < nNodesGlobal; in++) {
    isDirichletWall[in].resize(nDofsOnNode[in], false);
  }
	
  tmp = 0;
  for(int in = 0; in < nNodesGlobal; in++) {
    dofsMapWall[in].resize(nDofsOnNode[in]);
    dofsBCsMapWall[in].resize(nDofsOnNode[in]);
    for(int id = 0; id < nDofsOnNode[in]; id++) {
      dofsMapWall[in][id] = tmp;
      dofsBCsMapWall[in][id] = tmp;
      tmp++;
    }
  }
  /////////////

  for(int in = 0; in < nNodesGlobal; in++) {
    for(int d = 0; d < conf.dim; d++) {
      x[in][d] = conf.node[in][d];
    }
  }

  if(conf.gridType == GridType::STRUCTURED) {
    nNodesStrGlobal = (conf.nx + 1) * (conf.ny + 1) * (conf.nz + 1);
    VecTool::resize(vvti, nNodesStrGlobal, conf.dim);
    VecTool::resize(pvti, nNodesStrGlobal);
  }
}

void Node::assignCoordinates(Config &conf)
{
  for(int in = 0; in < conf.nNodesGlobal; in++) {
    for(int d = 0; d < conf.dim; d++) {
      x[in][d] = conf.node[in][d];
    }
  }
}

void Node::newMapping()
{
  // initialize variables
  nDofsOnNodeNew.resize(nNodesGlobal);
  isDirichletNew.resize(nNodesGlobal);
  mapNew.resize(nNodesGlobal);
  dofsMapNew.resize(nNodesGlobal);
  dofsBCsMapNew.resize(nNodesGlobal);

  int n1;
  for(int in = 0; in < nNodesGlobal; in++) {
    n1 = map[in];
    mapNew[n1] = in;
  }

  for(int in = 0; in < nNodesGlobal; in++) {
    nDofsOnNodeNew[mapNew[in]] = nDofsOnNode[in];
  }

  for(int in = 0; in < nNodesGlobal; in++) {
    isDirichletNew[in].resize(nDofsOnNodeNew[in], false);
  }

  int tmp = 0;
  for(int in = 0; in < nNodesGlobal; in++) {
    dofsMapNew[in].resize(nDofsOnNodeNew[in]);
    dofsBCsMapNew[in].resize(nDofsOnNodeNew[in]);
    for(int id = 0; id < nDofsOnNodeNew[in]; id++) {
      dofsMapNew[in][id] = tmp;
      dofsBCsMapNew[in][id] = tmp;
      tmp++;
    }
  }
}

void Node::initializeAdjoint(Config &conf, std::vector<int> &CBNodeMap)
{
  // initialize variables
  x.resize(nNodesGlobal, std::vector<double>(conf.dim));
  nDofsOnNode.resize(nNodesGlobal);
  isDirichlet.resize(nNodesGlobal);
  subId.resize(nNodesGlobal);
  map.resize(nNodesGlobal);
  dofsMap.resize(nNodesGlobal);
  dofsBCsMap.resize(nNodesGlobal);
  sortNode.resize(nNodesGlobal);

  /// add ///
  isDirichletWall.resize(nNodesGlobal);
  dofsMapWall.resize(nNodesGlobal);
  dofsBCsMapWall.resize(nNodesGlobal);
  ///////////
  
	for(int in = 0; in < nNodesGlobal; in++) {
    subId[in] = conf.nodeId[in];
  }
	
	nNodesLocal = count(conf.nodeId.begin(), conf.nodeId.end(), mpi.myId);

  for(int in = 0; in < nNodesGlobal; in++) {
    nDofsOnNode[in] = conf.dim + 1;
  }

  for(int ib = 0; ib < CBNodeMap.size(); ib++) {
    nDofsOnNode[CBNodeMap[ib]] += conf.dim;
  }

  for(int in = 0; in < nNodesGlobal; in++) {
    isDirichlet[in].resize(nDofsOnNode[in], false);
  }

  int tmp = 0;

  for(int in = 0; in < nNodesGlobal; in++) {
    dofsMap[in].resize(nDofsOnNode[in]);
    dofsBCsMap[in].resize(nDofsOnNode[in]);
    for(int id = 0; id < nDofsOnNode[in]; id++) {
      dofsMap[in][id] = tmp;
      dofsBCsMap[in][id] = tmp;
      tmp++;
    }
  }

  /// add ///
  for(int in = 0; in < nNodesGlobal; in++) {
    isDirichletWall[in].resize(nDofsOnNode[in], false);
  }
  tmp = 0;
  for(int in = 0; in < nNodesGlobal; in++) {
    dofsMapWall[in].resize(nDofsOnNode[in]);
    dofsBCsMapWall[in].resize(nDofsOnNode[in]);
    for(int id = 0; id < nDofsOnNode[in]; id++) {
      dofsMapWall[in][id] = tmp;
      dofsBCsMapWall[in][id] = tmp;
      tmp++;
    }
  }
  /////////////

  for(int in = 0; in < nNodesGlobal; in++) {
    for(int d = 0; d < conf.dim; d++) {
      x[in][d] = conf.node[in][d];
    }
  }

  if(conf.gridType == GridType::STRUCTURED) {
    nNodesStrGlobal = (conf.nx + 1) * (conf.ny + 1) * (conf.nz + 1);
    VecTool::resize(wvti, nNodesStrGlobal, conf.dim);
    VecTool::resize(lvti, nNodesStrGlobal, conf.dim);
    VecTool::resize(qvti, nNodesStrGlobal);
  }
}
