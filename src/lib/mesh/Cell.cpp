/**
 * @file Cell.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "Cell.h"

Cell::Cell() : nCellsGlobal(0), nCellsLocal(0), nNodesInCell(0), nStrCellsGlobal(0)
{
}

Cell::Cell(Config &conf) : nNodesInCell(conf.nNodesInCell), nCellsGlobal(conf.nCellsGlobal), data(conf.nCellsGlobal)
{
  if(conf.gridType == GridType::STRUCTURED) {
    nStrCellsGlobal = conf.nStrCellsGlobal;
  }
}

Cell::~Cell()
{
}

void Cell::initialize(Config &conf)
{
  initializeDataStructures(conf);
  assignNodes(conf);
  assignCoordinates(conf);
  assignPhi(conf);
  assignSubId(conf);
  assignCellType(conf);

  nCellsLocal = std::count(conf.cellId.begin(), conf.cellId.end(), mpi.myId);

  if(conf.gridType == GridType::STRUCTURED) {
    nCellsStrGlobal = conf.nx * conf.ny * conf.nz;
  }
}

void Cell::initializeAdjoint(Config &conf)
{
  initializeDataStructures(conf);
  assignNodes(conf);
  assignCoordinates(conf);
  assignPhi(conf);
  assignCellType(conf);

  if(conf.gridType == GridType::STRUCTURED) {
    nCellsStrGlobal = conf.nx * conf.ny * conf.nz;
  }
}

void Cell::initializeDataStructures(Config &conf)
{
  for(int ic = 0; ic < nCellsGlobal; ++ic) {
    data[ic].node.resize(conf.nNodesInCell);
    data[ic].nodeNew.resize(conf.nNodesInCell);
  }
}

void Cell::assignNodes(Config &conf)
{
  for(int ic = 0; ic < nCellsGlobal; ++ic) {
    data[ic].node.resize(conf.nNodesInCell);
    for(int p = 0; p < conf.nNodesInCell; ++p) {
      data[ic].node[p] = conf.cell[ic][p];
    }
  }
}

void Cell::assignCoordinates(Config &conf)
{
  for(int ic = 0; ic < nCellsGlobal; ++ic) {
    data[ic].x.resize(conf.nNodesInCell);
    for(int p = 0; p < conf.nNodesInCell; ++p) {
      data[ic].x[p].resize(conf.dim);
      for(int d = 0; d < conf.dim; ++d) {
        data[ic].x[p][d] = conf.node[data[ic].node[p]][d];
      }
    }
  }
}

void Cell::assignPhi(Config &conf)
{
  for(int ic = 0; ic < nCellsGlobal; ++ic) {
    data[ic].phi = conf.phi[ic];
  }
}

void Cell::assignSubId(Config &conf)
{
  for(int ic = 0; ic < nCellsGlobal; ++ic) {
    data[ic].subId = conf.cellId[ic];
  }
}

void Cell::assignCellType(Config &conf)
{
  if(conf.nNodesInCell == 4) {
    for(int ic = 0; ic < nCellsGlobal; ++ic) {
      data[ic].cellType = VTK_QUAD;
    }
  } else if(conf.nNodesInCell == 8) {
    for(int ic = 0; ic < nCellsGlobal; ++ic) {
      data[ic].cellType = VTK_HEXAHEDRON;
    }
  }
}