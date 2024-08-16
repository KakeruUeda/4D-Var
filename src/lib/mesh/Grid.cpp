/**
 * @file Grid.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "Grid.h"

Grid::Grid(Config &conf)
    : gridType(conf.gridType), cell(conf), node(conf), nNodesGlobal(conf.nNodesGlobal), nCellsGlobal(conf.nCellsGlobal),
      nDofsGlobal(0), nDofsLocal(0), dim(conf.dim), vecFluidUniqueNodes(conf.vecFluidUniqueNodes)
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

void Grid::prepareMatrix(Dirichlet &dirichletBC, PetscSolver &petsc, std::string outputDir, const int timeMax)
{
  nNodesLocal = node.nNodesLocal;
  nCellsLocal = cell.nCellsLocal;
  initialSetting(dirichletBC);

  if(mpi.nId == 1) {
    serialSetting();
  } else if(mpi.nId > 1) {
    parallelSetting(dirichletBC);
  }

  petscMatrixSetting(petsc);
  petsc.initialize(nDofsLocal, nDofsGlobal);
}

void Grid::initialSetting(Dirichlet &dirichletBC)
{
  for(auto &pair : dirichletBC.velocitySet) {
    int count = 0;
    for(auto &value : pair.second) {
      node.isDirichlet[pair.first][count] = true;
      count++;
    }
  }

  for(int in = 0; in < node.nNodesGlobal; in++) {
    for(int id = 0; id < node.nDofsOnNode[in]; id++) {
      if(node.isDirichlet[in][id]) {
        node.dofsBCsMap[in][id] = -1;
      }
    }
  }

  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    VecTool::resize(cell(ic).dofStart, cell.nNodesInCell);
    for(int p = 1; p < cell.nNodesInCell; p++) {
      cell(ic).dofStart[p] = cell(ic).dofStart[p - 1] + node.nDofsOnNode[cell(ic).node[p - 1]];
    }
  }

  nDofsGlobal = 0;
  for(int in = 0; in < node.nNodesGlobal; in++) {
    for(int id = 0; id < node.nDofsOnNode[in]; id++) {
      nDofsGlobal++;
    }
  }
}

void Grid::petscMatrixSetting(PetscSolver &petsc)
{
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
}

void Grid::serialSetting()
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

void Grid::parallelSetting(Dirichlet &dirichletBC)
{
  collectNodeMap();
  getNodeNewMap(dirichletBC);
  getCellNewMap();
  distributeToLocalDofs();
}

void Grid::collectNodeMap()
{
  int kk = 0;

  nodeStart = 0;
  nodeEnd = 0;

  std::vector<int> nodeListLocal(nNodesLocal);
  for(int in = 0; in < nNodesGlobal; in++) {
    if(node.subId[in] == mpi.myId) {
      nodeListLocal[kk++] = in;
    }
  }

  std::vector<int> nNodesLocalVector(mpi.nId);
  std::vector<int> nNodesLocalSum(mpi.nId);

  MPI_Allgather(&nNodesLocal, 1, MPI_INT, &nNodesLocalVector[0], 1, MPI_INT, MPI_COMM_WORLD);

  nNodesLocalSum = nNodesLocalVector;
  for(int i = 1; i < mpi.nId; i++) {
    nNodesLocalSum[i] += nNodesLocalSum[i - 1];
  }

  if(mpi.myId > 0) {
    nodeStart = nNodesLocalSum[mpi.myId - 1];
  }
  nodeEnd = nNodesLocalSum[mpi.myId] - 1;

  printf("nodeStart = %5d \t nodeEnd = %5d \t myId = %5d \n", nodeStart, nodeEnd, mpi.myId);

  std::vector<int> displs(mpi.nId);

  displs[0] = 0;
  for(int i = 0; i < mpi.nId - 1; i++) {
    displs[i + 1] = displs[i] + nNodesLocalVector[i];
  }

  MPI_Allgatherv(&nodeListLocal[0], nNodesLocal, MPI_INT, &node.map[0], &nNodesLocalVector[0], &displs[0], MPI_INT,
                 MPI_COMM_WORLD);
}

void Grid::distributeToLocalDofs()
{
  rowStart = 0;
  for(int in = 0; in < nodeStart; in++) {
    rowStart += node.nDofsOnNodeNew[in];
  }

  rowEnd = 0;
  for(int in = 0; in <= nodeEnd; in++) {
    rowEnd += node.nDofsOnNodeNew[in];
  }

  rowEnd = rowEnd - 1;

  nDofsLocal = 0;
  for(int in = nodeStart; in <= nodeEnd; in++) {
    for(int id = 0; id < node.nDofsOnNodeNew[in]; id++) {
      nDofsLocal++;
    }
  }

  printf("nDofsLocal = %5d/%5d \t rowStart  = %5d \t rowEnd  = %5d \t myId  = %5d \n", nDofsLocal, nDofsGlobal,
         rowStart, rowEnd, mpi.myId);

  int nDofsGlobalCheck;
  MPI_Allreduce(&nDofsLocal, &nDofsGlobalCheck, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(nDofsGlobalCheck != nDofsGlobal)
    std::cout << "Sum of local problem sizes is not equal to global size" << std::endl;
}

void Grid::getNodeNewMap(Dirichlet &dirichletBC)
{
  // initialize variables
  node.nDofsOnNodeNew.resize(node.nNodesGlobal);
  node.isDirichletNew.resize(node.nNodesGlobal);
  node.mapNew.resize(node.nNodesGlobal);
  node.dofsMapNew.resize(node.nNodesGlobal);
  node.dofsBCsMapNew.resize(node.nNodesGlobal);

  int n1;
  for(int in = 0; in < node.nNodesGlobal; in++) {
    n1 = node.map[in];
    node.mapNew[n1] = in;
  }

  for(int in = 0; in < node.nNodesGlobal; in++) {
    node.nDofsOnNodeNew[node.mapNew[in]] = node.nDofsOnNode[in];
  }

  for(int in = 0; in < node.nNodesGlobal; in++) {
    node.isDirichletNew[in].resize(node.nDofsOnNodeNew[in], false);
  }

  int tmp = 0;
  for(int in = 0; in < node.nNodesGlobal; in++) {
    node.dofsMapNew[in].resize(node.nDofsOnNodeNew[in]);
    node.dofsBCsMapNew[in].resize(node.nDofsOnNodeNew[in]);
    for(int id = 0; id < node.nDofsOnNodeNew[in]; id++) {
      node.dofsMapNew[in][id] = tmp;
      node.dofsBCsMapNew[in][id] = tmp;
      tmp++;
    }
  }

  int count;

  for(auto &pair : dirichletBC.velocitySet) {
    std::vector<double> vecTmp;
    n1 = node.mapNew[pair.first];
    count = 0;
    for(auto &value : pair.second) {
      vecTmp.push_back(value);
      node.isDirichletNew[n1][count] = true;
      count++;
    }
  }

  for(int in = 0; in < node.nNodesGlobal; in++) {
    for(int id = 0; id < node.nDofsOnNode[in]; id++) {
      if(node.isDirichletNew[in][id])
        node.dofsBCsMapNew[in][id] = -1;
    }
  }
}

void Grid::getCellNewMap()
{
  for(int ic = 0; ic < nCellsGlobal; ic++) {
    for(int p = 0; p < cell.nNodesInCell; p++) {
      cell(ic).nodeNew[p] = node.mapNew[cell(ic).node[p]];
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
}