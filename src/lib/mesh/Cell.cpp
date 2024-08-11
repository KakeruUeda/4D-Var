/**
 * @file Cell.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "Cell.h"

inline void CellInfo::setArrayZero(int n)
{
  assert(n != 0);
}

void Cell::initialize(Config &conf)
{
  for (int ic = 0; ic < nCellsGlobal; ic++)
    data[ic].node.resize(conf.nNodesInCell);

  for (int ic = 0; ic < nCellsGlobal; ic++)
    data[ic].nodeNew.resize(conf.nNodesInCell);

  for (int ic = 0; ic < nCellsGlobal; ic++)
    data[ic].x.resize(conf.nNodesInCell,
                      std::vector<double>(conf.dim));

  for (int ic = 0; ic < nCellsGlobal; ic++)
    for (int p = 0; p < conf.nNodesInCell; p++)
      data[ic].node[p] = conf.cell[ic][p];

  for (int ic = 0; ic < nCellsGlobal; ic++)
    for (int p = 0; p < conf.nNodesInCell; p++)
      for (int d = 0; d < conf.dim; d++)
        data[ic].x[p][d] = conf.node[data[ic].node[p]][d];

  for (int ic = 0; ic < nCellsGlobal; ic++)
    data[ic].phi = conf.phi[ic];

  for (int ic = 0; ic < nCellsGlobal; ic++)
    data[ic].subId = conf.cellId[ic];

  nCellsLocal = count(conf.cellId.begin(), conf.cellId.end(), mpi.myId);

  if (conf.nNodesInCell == 4)
    for (int ic = 0; ic < nCellsGlobal; ic++)
      data[ic].cellType = VTK_QUAD;
  if (conf.nNodesInCell == 8)
    for (int ic = 0; ic < nCellsGlobal; ic++)
      data[ic].cellType = VTK_HEXAHEDRON;

  if (conf.gridType == GridType::STRUCTURED)
  {
    nCellsStrGlobal = conf.nx * conf.ny * conf.nz;
  }
}

void Cell::initializeAdjoint(Config &conf)
{
  for (int ic = 0; ic < nCellsGlobal; ic++)
    data[ic].node.resize(conf.nNodesInCell);

  for (int ic = 0; ic < nCellsGlobal; ic++)
    data[ic].nodeNew.resize(conf.nNodesInCell);

  for (int ic = 0; ic < nCellsGlobal; ic++)
    data[ic].x.resize(conf.nNodesInCell,
                      std::vector<double>(conf.dim));

  for (int ic = 0; ic < nCellsGlobal; ic++)
    for (int p = 0; p < conf.nNodesInCell; p++)
      data[ic].node[p] = conf.cell[ic][p];

  for (int ic = 0; ic < nCellsGlobal; ic++)
    for (int p = 0; p < conf.nNodesInCell; p++)
      for (int d = 0; d < conf.dim; d++)
        data[ic].x[p][d] = conf.node[data[ic].node[p]][d];

  for (int ic = 0; ic < nCellsGlobal; ic++)
    data[ic].phi = conf.phi[ic];

  if (conf.nNodesInCell == 4)
    for (int ic = 0; ic < nCellsGlobal; ic++)
      data[ic].cellType = VTK_QUAD;
  if (conf.nNodesInCell == 8)
    for (int ic = 0; ic < nCellsGlobal; ic++)
      data[ic].cellType = VTK_HEXAHEDRON;

  if (conf.gridType == GridType::STRUCTURED)
  {
    nCellsStrGlobal = conf.nx * conf.ny * conf.nz;
  }
}