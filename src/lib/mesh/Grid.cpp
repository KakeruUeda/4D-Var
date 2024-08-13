/**
 * @file Grid.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "Grid.h"

Grid::Grid(Config &conf)
    : gridType(conf.gridType), cell(conf), node(conf), dirichlet(conf), nNodesGlobal(conf.nNodesGlobal),
      nCellsGlobal(conf.nCellsGlobal), nDofsGlobal(0), nDofsLocal(0), dim(conf.dim),
      vecFluidUniqueNodes(conf.vecFluidUniqueNodes)
{
  if(gridType == GridType::STRUCTURED) {
    nx = conf.nx;
    ny = conf.ny;
    nz = conf.nz;
    lx = conf.lx;
    ly = conf.ly;
    lz = conf.lz;
    dx = conf.dx;
    dy = conf.dy;
    dz = conf.dz;
  }
}

void Grid::prepareMatrix(PetscSolver &petsc, std::string outputDir, const int timeMax)
{
  for(auto &pair : dirichlet.vDirichlet[0]) {
    int count = 0;
    for(auto &value : pair.second) {
      node.isDirichlet[pair.first][count] = true;
      count++;
    }
  }

  for(auto &pair : dirichlet.pDirichlet[0])
    node.isDirichlet[pair.first][dim] = true;

  for(int in = 0; in < node.nNodesGlobal; in++)
    for(int id = 0; id < node.nDofsOnNode[in]; id++)
      if(node.isDirichlet[in][id])
        node.dofsBCsMap[in][id] = -1;

  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    VecTool::resize(cell(ic).dofStart, cell.nNodesInCell);
    for(int p = 1; p < cell.nNodesInCell; p++) {
      cell(ic).dofStart[p] = cell(ic).dofStart[p - 1] + node.nDofsOnNode[cell(ic).node[p - 1]];
    }
  }

  //// add //////////////////////////////////////
  for(auto &pair : dirichlet.vDirichletWall[0]) {
    int count = 0;
    for(auto &value : pair.second) {
      node.isDirichletWall[pair.first][count] = true;
      count++;
    }
  }
  for(int in = 0; in < node.nNodesGlobal; in++)
    for(int id = 0; id < node.nDofsOnNode[in]; id++)
      if(node.isDirichletWall[in][id])
        node.dofsBCsMapWall[in][id] = -1;
  ////////////////////////////////////////////////

  nDofsGlobal = 0;
  for(int in = 0; in < node.nNodesGlobal; in++) {
    for(int id = 0; id < node.nDofsOnNode[in]; id++) {
      nDofsGlobal++;
    }
  }

  if(mpi.nId == 1) {
    setForSerial();
  } else if(mpi.nId > 1) {
    distributeToLocal(timeMax);
  }

  for(int in = 0; in < node.nNodesGlobal; in++) {
    for(int id = 0; id < node.nDofsOnNodeNew[in]; id++) {
      node.dofsMapNew1D.push_back(node.dofsMapNew[in][id]);
      node.dofsBCsMapNew1D.push_back(node.dofsBCsMapNew[in][id]);
    }
  }

  int size;
  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    if(cell(ic).subId == mpi.myId) {
      size = 0;
      for(int p = 0; p < cell.nNodesInCell; p++) {
        size += node.nDofsOnNodeNew[cell(ic).nodeNew[p]];
      }
      cell(ic).dofsMap.resize(size);
      cell(ic).dofsBCsMap.resize(size);
      int i = 0;
      int j = 0;
      for(int p = 0; p < cell.nNodesInCell; p++) {
        j = cell(ic).nodeNew[p];
        for(int q = 0; q < node.nDofsOnNodeNew[cell(ic).nodeNew[p]]; q++) {
          cell(ic).dofsMap[i + q] = node.dofsMapNew[j][q];
          cell(ic).dofsBCsMap[i + q] = node.dofsBCsMapNew[j][q];
        }
        i += node.nDofsOnNodeNew[cell(ic).nodeNew[p]];
      }
    }
  }

  /// add ///
  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    if(cell(ic).subId == mpi.myId) {
      size = 0;
      for(int p = 0; p < cell.nNodesInCell; p++) {
        size += node.nDofsOnNodeNew[cell(ic).nodeNew[p]];
      }
      cell(ic).dofsMapWall.resize(size);
      cell(ic).dofsBCsMapWall.resize(size);
      int i = 0;
      int j = 0;
      for(int p = 0; p < cell.nNodesInCell; p++) {
        j = cell(ic).nodeNew[p];
        for(int q = 0; q < node.nDofsOnNodeNew[cell(ic).nodeNew[p]]; q++) {
          cell(ic).dofsMapWall[i + q] = node.dofsMapWallNew[j][q];
          cell(ic).dofsBCsMapWall[i + q] = node.dofsBCsMapWallNew[j][q];
        }
        i += node.nDofsOnNodeNew[cell(ic).nodeNew[p]];
      }
    }
  }
  //////////

  int *tt, tmpInt;
  int r, kk, nSize;
  int countDiag, countOffDiag;

  std::vector<std::set<int>> forAssyMatFluid;
  std::set<int>::iterator it;

  forAssyMatFluid.resize(nDofsGlobal);

  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    if(cell(ic).subId == mpi.myId) {
      tt = &(cell(ic).dofsMap[0]);
      nSize = cell(ic).dofsMap.size();

      for(int i = 0; i < nSize; i++) {
        r = tt[i];
        if(r != -1) {
          if(r >= rowStart && r <= rowEnd) {
            for(int j = 0; j < nSize; j++) {
              if(tt[j] != -1) {
                forAssyMatFluid[r].insert(tt[j]);
              }
            }
          }
        }
      }
    }
  }
  PetscMalloc1(nDofsLocal, &petsc.diag_nnz);
  PetscMalloc1(nDofsLocal, &petsc.offdiag_nnz);

  kk = 0;
  petsc.nnz_max_row = 0;
  for(int i = rowStart; i <= rowEnd; i++) {
    nSize = forAssyMatFluid[i].size();
    petsc.nnz_max_row = std::max(petsc.nnz_max_row, nSize);
    countDiag = 0, countOffDiag = 0;
    for(it = forAssyMatFluid[i].begin(); it != forAssyMatFluid[i].end(); it++) {
      tmpInt = *it;
      if(tmpInt >= rowStart && tmpInt <= rowEnd)
        countDiag++;
      else
        countOffDiag++;
    }
    petsc.diag_nnz[kk] = countDiag;
    petsc.offdiag_nnz[kk] = countOffDiag;
    kk++;
  }

  petsc.initialize(nDofsLocal, nDofsGlobal);
}

void Grid::setForSerial()
{
  nCellsLocal = cell.nCellsGlobal;
  nNodesLocal = node.nNodesGlobal;
  nDofsLocal = nDofsGlobal;
  rowStart = 0;
  rowEnd = nDofsGlobal - 1;

  node.isDirichletNew = node.isDirichlet;
  node.dofsBCsMapNew = node.dofsBCsMap;
  node.dofsMapNew = node.dofsMap;
  node.nDofsOnNodeNew = node.nDofsOnNode;
  dirichlet.vDirichletNew = dirichlet.vDirichletNew;

  for(int ic = 0; ic < cell.nCellsGlobal; ic++)
    cell(ic).nodeNew = cell(ic).node;
}

void Grid::distributeToLocal(const int timeMax)
{
  int kk = 0;
  std::vector<int> nodeListLocal(node.nNodesLocal);

  for(int in = 0; in < nNodesGlobal; in++) {
    if(node.subId[in] == mpi.myId) {
      nodeListLocal[kk++] = in;
    }
  }

  std::vector<int> nNodesLocalVector(mpi.nId);
  std::vector<int> nNodesLocalSum(mpi.nId);

  MPI_Allgather(&node.nNodesLocal, 1, MPI_INT, &nNodesLocalVector[0], 1, MPI_INT, MPI_COMM_WORLD);

  nNodesLocalSum = nNodesLocalVector;
  for(int i = 1; i < mpi.nId; i++)
    nNodesLocalSum[i] += nNodesLocalSum[i - 1];

  int nodeStart = 0;
  int nodeEnd = 0;

  if(mpi.myId > 0)
    nodeStart = nNodesLocalSum[mpi.myId - 1];
  nodeEnd = nNodesLocalSum[mpi.myId] - 1;

  printf("nodeStart = %5d \t nodeEnd = %5d \t myId = %5d \n", nodeStart, nodeEnd, mpi.myId);

  std::vector<int> displs(mpi.nId);

  displs[0] = 0;
  for(int i = 0; i < mpi.nId - 1; i++)
    displs[i + 1] = displs[i] + nNodesLocalVector[i];

  std::vector<int> tmp(nNodesGlobal);

  MPI_Allgatherv(&nodeListLocal[0], node.nNodesLocal, MPI_INT, &node.map[0], &nNodesLocalVector[0], &displs[0], MPI_INT,
                 MPI_COMM_WORLD);

  node.initializeNew();

  for(int ic = 0; ic < nCellsGlobal; ic++)
    for(int p = 0; p < cell.nNodesInCell; p++)
      cell(ic).nodeNew[p] = node.mapNew[cell(ic).node[p]];

  int n1;
  int count;

  dirichlet.vDirichletNew.resize(timeMax);
  dirichlet.pDirichletNew.resize(timeMax);

  /// add ///
  dirichlet.vDirichletWallNew.resize(timeMax);
  //////////

  for(int t = 0; t < timeMax; t++) {
    for(auto &pair : dirichlet.vDirichlet[t]) {
      std::vector<double> vecTmp;
      n1 = node.mapNew[pair.first];
      count = 0;
      for(auto &value : pair.second) {
        vecTmp.push_back(value);
        if(t == 0)
          node.isDirichletNew[n1][count] = true;
        count++;
      }
      dirichlet.vDirichletNew[t][n1] = vecTmp;
    }
    /// add ///
    for(auto &pair : dirichlet.vDirichletWall[t]) {
      std::vector<double> vecTmp;
      n1 = node.mapNew[pair.first];
      count = 0;
      for(auto &value : pair.second) {
        vecTmp.push_back(value);
        if(t == 0)
          node.isDirichletWallNew[n1][count] = true;
        count++;
      }
      dirichlet.vDirichletWallNew[t][n1] = vecTmp;
    }
    //////////

    for(auto &pair : dirichlet.pDirichlet[t]) {
      n1 = node.mapNew[pair.first];
      dirichlet.pDirichletNew[t][n1] = pair.second;
      if(t == 0)
        node.isDirichletNew[n1][dim] = true;
    }
  }

  for(int in = 0; in < node.nNodesGlobal; in++)
    for(int id = 0; id < node.nDofsOnNodeNew[in]; id++)
      if(node.isDirichletNew[in][id])
        node.dofsBCsMapNew[in][id] = -1;

  /// add ///
  for(int in = 0; in < node.nNodesGlobal; in++)
    for(int id = 0; id < node.nDofsOnNodeNew[in]; id++)
      if(node.isDirichletWallNew[in][id])
        node.dofsBCsMapWallNew[in][id] = -1;
  ///////////

  rowStart = 0;
  for(int in = 0; in < nodeStart; in++)
    rowStart += node.nDofsOnNodeNew[in];

  rowEnd = 0;
  for(int in = 0; in <= nodeEnd; in++)
    rowEnd += node.nDofsOnNodeNew[in];

  rowEnd = rowEnd - 1;

  for(int in = nodeStart; in <= nodeEnd; in++)
    for(int id = 0; id < node.nDofsOnNodeNew[in]; id++)
      nDofsLocal++;

  printf("nDofsLocal = %5d/%5d \t rowStart  = %5d \t rowEnd  = %5d \t myId  = %5d \n", nDofsLocal, nDofsGlobal,
         rowStart, rowEnd, mpi.myId);

  int nDofsGlobalCheck;
  MPI_Allreduce(&nDofsLocal, &nDofsGlobalCheck, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(nDofsGlobalCheck != nDofsGlobal)
    std::cout << "Sum of local problem sizes is not equal to global size" << std::endl;
}
