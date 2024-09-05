/**
 * @file GridCreation.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "GridCreation.h"

/**
 * @brief Constructor.
 */
GridCreation::GridCreation(Config &conf)
{
  outputDir = conf.outputDir;
  std::string output = "output";
  mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  outputDir = "output/" + conf.outputDir;
  mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

  std::string dir;
  dir = outputDir + "/dat";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  dir = outputDir + "/vtu";
  mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}

/**
 * @brief Initialize grid.
 */
void GridCreation::initialize(Config &conf)
{
  grid.gridType = conf.gridType;
  initializeCells(conf);
  initializeNodes(conf);

  vDirichlet = conf.vDirichlet;
  pDirichlet = conf.pDirichlet;

	fluidUniqueNodes = conf.fluidUniqueNodes;
	phi = conf.phi;

  CBExtraction = conf.CBExtraction;
  if(CBExtraction == ON) {
    CBNodeMap = conf.CBNodeMap;
    CBCellMap = conf.CBCellMap;
    CBNodeMapInCell = conf.CBNodeMapInCell;
  }
}

/**
 * @brief Initialize CFD grid cells.
 */
void GridCreation::initializeCells(Config &conf)
{
  grid.cell.nNodesInCell = conf.nNodesInCell;
  grid.cell.nCellsGlobal = conf.nCellsGlobal;
  grid.cell.resize(conf.nCellsGlobal);

  grid.cell.assignNodes(conf);
  grid.cell.assignCoordinates(conf);
  grid.cell.assignCellType(conf);

  if(conf.gridType == GridType::STRUCTURED) {
    grid.cell.nCellsStrGlobal = conf.nx * conf.ny * conf.nz;
  }
}

/**
 * @brief Initialize CFD grid nodes.
 */
void GridCreation::initializeNodes(Config &conf)
{
  grid.node.nNodesGlobal = conf.nNodesGlobal;
  grid.node.x.resize(conf.nNodesGlobal, std::vector<double>(conf.dim));
  grid.node.sortNode.resize(conf.nNodesGlobal);
	grid.node.assignCoordinates(conf);

  cellId.resize(grid.cell.nCellsGlobal, 0);
  nodeId.resize(grid.node.nNodesGlobal, 0);
}

/**
 * @brief Initialize CFD cell type.
 */
void GridCreation::setCellType(int nNodesInCell)
{
  if(nNodesInCell == 4) {
    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      grid.cell(ic).cellType = VTK_QUAD;
    }
  } else if(nNodesInCell == 8) {
    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      grid.cell(ic).cellType = VTK_HEXAHEDRON;
    }
  }
}

/**
 * @brief Divide whole grid.
 */
void GridCreation::divideWholeGrid()
{
  if(mpi.myId != 0) {
    return;
  }

  int nparts = mpi.nId;

  std::vector<int> eptr(grid.cell.nCellsGlobal + 1, 0);
  std::vector<int> eind(grid.cell.nCellsGlobal * grid.cell.nNodesInCell, 0);

  int kk = 0;
  for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
    eptr[ic + 1] = (ic + 1) * grid.cell.nNodesInCell;
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      eind[kk + p] = grid.cell(ic).node[p];
    }
    kk += grid.cell.nNodesInCell;
  }

  partitionMesh(eptr.data(), eind.data(), nparts);
}

/**
 * @brief Partition mesh.
 */
void GridCreation::partitionMesh(int *eptr, int *eind, int nparts)
{
  int ncommon_nodes = 4;
  idx_t objval;
  idx_t options[METIS_NOPTIONS];

  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
  options[METIS_OPTION_NUMBERING] = 0;

  int ret = METIS_PartMeshDual(&grid.cell.nCellsGlobal, &grid.node.nNodesGlobal, eptr, eind, NULL, NULL, &ncommon_nodes,
                               &nparts, NULL, options, &objval, &cellId[0], &nodeId[0]);

  if(ret == METIS_OK) {
    std::cout << "METIS partition routine success " << std::endl;
  } else {
    std::cout << "METIS partition routine failed " << std::endl;
  }
}

/**
 * @brief Collect local grid.
 */
void GridCreation::collectLocalGrid()
{
  MPI_Bcast(cellId.data(), grid.cell.nCellsGlobal, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(nodeId.data(), grid.node.nNodesGlobal, MPI_INT, 0, MPI_COMM_WORLD);

  grid.nCellsLocal = std::count(cellId.begin(), cellId.end(), mpi.myId);
  grid.nNodesLocal = std::count(nodeId.begin(), nodeId.end(), mpi.myId);

  printf("nCellsLocal =  %5d \t numOfId = %5d \t myId = %5d \n", grid.nCellsLocal, mpi.nId, mpi.myId);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("nNodesLocal =  %5d \t numOfId = %5d \t myId = %5d \n", grid.nNodesLocal, mpi.nId, mpi.myId);
}

/**
 * @brief Output
 */
void GridCreation::outputDat()
{
  if(mpi.myId != 0) {
    return;
  }
  std::string datFile;

  datFile = outputDir + "/dat/cellId.dat";
  EXPORT::exportScalarDataDAT<int>(datFile, cellId);
  datFile = outputDir + "/dat/nodeId.dat";
  EXPORT::exportScalarDataDAT<int>(datFile, nodeId);
  datFile = outputDir + "/dat/image.dat";
  EXPORT::exportScalarDataDAT<double>(datFile, phi);
  datFile = outputDir + "/dat/velocityDirichlet.dat";
  EXPORT::exportMapVectorDataDAT<double>(datFile, vDirichlet);
  datFile = outputDir + "/dat/pressureDirichlet.dat";
  EXPORT::exportMapScalarDataDAT<double>(datFile, pDirichlet);
	datFile = outputDir + "/dat/fluidUniqueNodes.dat";
	EXPORT::exportScalarDataDAT<int>(datFile, fluidUniqueNodes);

  datFile = outputDir + "/dat/cell.dat";
  EXPORT::exportCellDataDAT(datFile, grid.cell);
  datFile = outputDir + "/dat/node.dat";
  EXPORT::exportNodeDataDAT(datFile, grid.node);

  if(CBExtraction == ON) {
    datFile = outputDir + "/dat/controlBoundaryNodeMap.dat";
    EXPORT::exportScalarDataDAT<int>(datFile, CBNodeMap);
    datFile = outputDir + "/dat/controlBoundaryCellMap.dat";
    EXPORT::exportScalarDataDAT<int>(datFile, CBCellMap);
    datFile = outputDir + "/dat/controlBoundaryNodeMapInCell.dat";
    EXPORT::exportVectorDataDAT<int>(datFile, CBNodeMapInCell);
  }
}

/**
 * @brief Output
 */
void GridCreation::outputVTU()
{
  if(mpi.myId != 0) {
    return;
  }
  std::string vtuFile;

  vtuFile = outputDir + "/vtu/cellId.vtu";
  EXPORT::exportScalarCellDataVTU<int>(vtuFile, "cellId", grid.node, grid.cell, cellId);
  vtuFile = outputDir + "/vtu/nodeId.vtu";
  EXPORT::exportScalarPointDataVTU<int>(vtuFile, "nodeId", grid.node, grid.cell, nodeId);
	vtuFile = outputDir + "/vtu/image.vtu";
	EXPORT::exportScalarCellDataVTU<double>(vtuFile, "image", grid.node, grid.cell, phi);

  std::vector<std::vector<double>> dirichlet;
  dirichlet.resize(grid.node.nNodesGlobal, std::vector<double>(3, 0e0));
  for (auto &entry : vDirichlet)
  {
     dirichlet[entry.first] = entry.second;
  }
  vtuFile = outputDir + "/vtu/vDirichlet.vtu";
  EXPORT::exportVectorPointDataVTU<double>(vtuFile, "vDirichlet", grid.node,grid.cell, dirichlet);
}

/*
 * std::vector<std::vector<double>> dirichlet;
 * dirichlet.resize(grid.node.nNodesGlobal, std::vector<double>(3, 0e0));
 * for (auto &entry : vDirichlet)
 * {
 *   dirichlet[entry.first] = entry.second;
 * }
 * vtuFile = outputDir + "/vtu/vDirichlet.vtu";
 * EXPORT::exportVectorPointDataVTU<double>(vtuFile, "vDirichlet", grid.node,
 * grid.cell, dirichlet);
 */