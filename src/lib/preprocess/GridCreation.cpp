/**
 * @file GridCreation.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "GridCreation.h"


/**********************
 * @brief tmp function.
 */
void GridCreation::initialize(Config &conf)
{
  // initialize grid
  grid.gridType = conf.gridType;

  // initialize cell
  grid.cell.nNodesInCell = conf.nNodesInCell;
  grid.cell.nCellsGlobal = conf.nCellsGlobal;

  grid.cell.resize(conf.nCellsGlobal);

  for (int ic = 0; ic < conf.nCellsGlobal; ic++)
  {
    grid.cell(ic).node.resize(conf.nNodesInCell);
  }

  for (int ic = 0; ic < conf.nCellsGlobal; ic++)
  {
    grid.cell(ic).x.resize(conf.nNodesInCell, std::vector<double>(conf.dim));
  }

  for (int ic = 0; ic < conf.nCellsGlobal; ic++)
  {
    for (int p = 0; p < conf.nNodesInCell; p++)
    {
      grid.cell(ic).node[p] = conf.cell[ic][p];
    }
  }

  for (int ic = 0; ic < conf.nCellsGlobal; ic++)
  {
    for (int p = 0; p < conf.nNodesInCell; p++)
    {
      for (int d = 0; d < conf.dim; d++)
      {
        grid.cell(ic).x[p][d] = conf.node[grid.cell(ic).node[p]][d];
      }
    }
  }

  for (int ic = 0; ic < conf.nCellsGlobal; ic++)
  {
    grid.cell(ic).phi = conf.phi[ic];
  }

  if (conf.nNodesInCell == 4)
  {
    for (int ic = 0; ic < conf.nCellsGlobal; ic++)
    {
      grid.cell(ic).cellType = VTK_QUAD;
    }
  }
  else if (conf.nNodesInCell == 8)
  {
    for (int ic = 0; ic < conf.nCellsGlobal; ic++)
    {
      grid.cell(ic).cellType = VTK_HEXAHEDRON;
    }
  }

  if (conf.gridType == GridType::STRUCTURED)
  {
    grid.cell.nCellsStrGlobal = conf.nx * conf.ny * conf.nz;
  }

  // initialize node
  grid.node.nNodesGlobal = conf.nNodesGlobal;
  grid.node.x.resize(conf.nNodesGlobal, std::vector<double>(conf.dim));
  grid.node.sortNode.resize(conf.nNodesGlobal);

  for (int in = 0; in < conf.nNodesGlobal; in++)
  {
    for (int d = 0; d < conf.dim; d++)
    {
      grid.node.x[in][d] = conf.node[in][d];
    }
    //grid.node.sortNode[in] = conf.sortNode[in];
  }

  cellId.resize(grid.cell.nCellsGlobal, 0);
  nodeId.resize(grid.node.nNodesGlobal, 0);

  // initialize dirichlet
  vDirichlet = conf.vDirichlet;
  pDirichlet = conf.pDirichlet;
}

/**********************
 * @brief tmp function.
 */
void GridCreation::divideWholeGrid()
{
  if (mpi.myId != 0) return;

  int kk = 0;
  int nparts = mpi.nId;

  std::unique_ptr<int[]> eptr;
  std::unique_ptr<int[]> eind;

  eptr = std::make_unique<int[]>(grid.cell.nCellsGlobal + 1);
  eind = std::make_unique<int[]>(grid.cell.nCellsGlobal * grid.cell.nNodesInCell);

  for (int in = 0; in < grid.cell.nCellsGlobal + 1; in++)
    eptr[in] = 0;

  for (int ic = 0; ic < grid.cell.nCellsGlobal; ic++)
    for (int p = 0; p < grid.cell.nNodesInCell; p++)
      eind[ic + p] = 0;

  for (int ic = 0; ic < grid.cell.nCellsGlobal; ic++)
  {
    eptr[ic + 1] = (ic + 1) * grid.cell.nNodesInCell;
    for (int p = 0; p < grid.cell.nNodesInCell; p++)
      eind[kk + p] = grid.cell(ic).node[p];
    kk += grid.cell.nNodesInCell;
  }
  // 8-noded hexa element
  int ncommon_nodes = 4;

  idx_t objval;
  idx_t options[METIS_NOPTIONS];

  METIS_SetDefaultOptions(options);

  // Specifies the partitioning method.
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;

  // Total communication volume minimization
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

  // C-style numbering is assumed that starts from 0
  options[METIS_OPTION_NUMBERING] = 0;

  // METIS partition routine
  int ret = METIS_PartMeshDual(&grid.cell.nCellsGlobal, &grid.node.nNodesGlobal, 
                               eptr.get(), eind.get(), NULL, NULL, &ncommon_nodes, 
                               &nparts, NULL, options, &objval, &cellId[0], &nodeId[0]);

  if (ret == METIS_OK)
  {
    std::cout << "METIS partition routine success " << std::endl;
  }
  else
  {
    std::cout << "METIS partition routine failed " << std::endl;
  }
}

/**********************
 * @brief tmp function.
 */
void GridCreation::collectLocalGrid()
{
  MPI_Bcast(cellId.data(), grid.cell.nCellsGlobal, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(nodeId.data(), grid.node.nNodesGlobal, MPI_INT, 0, MPI_COMM_WORLD);

  grid.nCellsLocal = count(cellId.begin(), cellId.end(), mpi.myId);
  grid.nNodesLocal = count(nodeId.begin(), nodeId.end(), mpi.myId);

  printf("nCellsLocal =  %5d \t numOfId = %5d \t myId = %5d \n", grid.nCellsLocal, mpi.nId, mpi.myId);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("nNodesLocal =  %5d \t numOfId = %5d \t myId = %5d \n", grid.nNodesLocal, mpi.nId, mpi.myId);
}

/**********************
 * @brief tmp function.
 */
void GridCreation::outputDat()
{
  if(mpi.myId != 0) return;

  std::string datFile;
  datFile = outputDir + "/cellId.dat";
  DAT::exportScalarDataDAT<int>(datFile, cellId);
  datFile = outputDir + "/nodeId.dat";
  DAT::exportScalarDataDAT<int>(datFile, nodeId);
}

/**********************
 * @brief tmp function.
 */
void GridCreation::outputVTU()
{
  if(mpi.myId != 0) return;

  std::string vtuFile;
  vtuFile = outputDir + "/cellId.vtu";
  VTKTMP::exportScalarCellDataVTU<int>(vtuFile, "cellId", grid.node, grid.cell, cellId);
  vtuFile = outputDir + "/nodeId.vtu";
  VTKTMP::exportScalarPointDataVTU<int>(vtuFile, "nodeId", grid.node, grid.cell, nodeId);
  vtuFile = outputDir + "/phi.vtu";
  VTK::exportPhiVTU(vtuFile, grid.node, grid.cell);

  std::vector<std::vector<double>> dirichlet;
  dirichlet.resize(grid.node.nNodesGlobal, std::vector<double>(3, 0e0));
  for (auto &entry : vDirichlet)
  {
    dirichlet[entry.first] = entry.second;
  }
  vtuFile = outputDir + "/vDirichlet.vtu";
  VTK::exportVectorPointDataVTU(vtuFile, "vDirichlet", grid.node, grid.cell, dirichlet);
}