/**
 * @file Boundary.cpp
 * @author K.Ueda
 * @date May, 2024
 */

#include "Boundary.h"

/***************************************
 * @brief Initialize dirichlet boundary.
 */
void DirichletBoundary::initialize(Config &conf)
{
  vDirichlet.resize(conf.timeMax);
  vDirichletWall.resize(conf.timeMax);
  pDirichlet.resize(conf.timeMax);
  for(int t = 0; t < conf.timeMax; t++) {
    vDirichlet[t] = conf.vDirichlet;
    vDirichletWall[t] = conf.vDirichletWall;
    pDirichlet[t] = conf.pDirichlet;
  }
}

/***************************************
 * @brief Initialize dirichlet boundary
 *        for adjoint equation.
 */
void DirichletBoundary::initializeAdjoint(Config &conf)
{
  vDirichlet.resize(conf.timeMax);
  vDirichletWall.resize(conf.timeMax);
  pDirichlet.resize(conf.timeMax);
  for(int t = 0; t < conf.timeMax; t++) {
    vDirichlet[t] = conf.vDirichlet;
    vDirichletWall[t] = conf.vDirichletWall;
    pDirichlet[t] = conf.pDirichlet;
  }

  controlBoundaryMap = conf.controlBoundaryMap;
  controlCellMap = conf.controlCellMap;
  controlNodeInCell = conf.controlNodeInCell;
  nControlNodesInCell = conf.controlNodeInCell[0].size();
  nControlCellsGlobal = controlCellMap.size();
  nControlNodesGlobal = controlBoundaryMap.size();
  isBoundaryEdge = conf.isBoundaryEdge;
}

/***************************************************
 * @brief Assign dirichlet Boundary condition's
 *        value and dof location to vector variable.
 */
void DirichletBoundary::assignDirichletBCs(std::vector<std::map<int, std::vector<double>>> &vDirichletNew,
                                           std::vector<std::map<int, double>> &pDirichletNew, Node &node, int &dim,
                                           const int t)
{
  int dofCurrentTmp, dofCurrent;
  int count;

  for(auto &pair : vDirichletNew[0]) {
    dofCurrentTmp = 0;
    dofCurrent = 0;
    count = 0;
    for(int i = 0; i < pair.first; i++)
      dofCurrentTmp += node.nDofsOnNodeNew[i];
    for(auto &value : pair.second) {
      dofCurrent = dofCurrentTmp + count;
      dirichletBCsValueNewInit[dofCurrent] = value;
      count++;
    }
  }

  for(auto &pair : vDirichletNew[t]) {
    dofCurrentTmp = 0;
    dofCurrent = 0;
    count = 0;
    for(int i = 0; i < pair.first; i++)
      dofCurrentTmp += node.nDofsOnNodeNew[i];
    for(auto &value : pair.second) {
      dofCurrent = dofCurrentTmp + count;
      dirichletBCsValueNew[dofCurrent] = value;
      count++;
    }
  }

  for(auto &pair : pDirichletNew[0]) {
    dofCurrentTmp = 0;
    dofCurrent = 0;
    count = 0;
    for(int i = 0; i < pair.first; i++)
      dofCurrentTmp += node.nDofsOnNodeNew[i];
    dofCurrent = dofCurrentTmp + dim;
    dirichletBCsValueNewInit[dofCurrent] = pair.second;
  }

  for(auto &pair : pDirichletNew[t]) {
    dofCurrentTmp = 0;
    dofCurrent = 0;
    count = 0;
    for(int i = 0; i < pair.first; i++)
      dofCurrentTmp += node.nDofsOnNodeNew[i];
    dofCurrent = dofCurrentTmp + dim;
    dirichletBCsValueNew[dofCurrent] = pair.second;
  }
}

/******************************************************
 * @brief Assign constant dirichlet boundary conditions
 *        over all simulation time steps.
 */
void DirichletBoundary::assignConstantDirichletBCs(std::vector<std::map<int, std::vector<double>>> &vDirichletNew,
                                                   std::vector<std::map<int, double>> &pDirichletNew, Node &node,
                                                   int &dim, const int t)
{
  int dofCurrentTmp, dofCurrent;
  int count;

  for(auto &pair : vDirichletNew[0]) {
    dofCurrentTmp = 0;
    dofCurrent = 0;
    count = 0;
    for(int i = 0; i < pair.first; i++)
      dofCurrentTmp += node.nDofsOnNodeNew[i];
    for(auto &value : pair.second) {
      dofCurrent = dofCurrentTmp + count;
      dirichletBCsValueNew[dofCurrent] = value;
      count++;
    }
  }

  for(auto &pair : pDirichletNew[0]) {
    dofCurrentTmp = 0;
    dofCurrent = 0;
    count = 0;
    for(int i = 0; i < pair.first; i++)
      dofCurrentTmp += node.nDofsOnNodeNew[i];
    dofCurrent = dofCurrentTmp + dim;
    dirichletBCsValueNew[dofCurrent] = pair.second;
  }
}

/*******************************************************
 * @brief Assign pulsatile dirichlet boundary condition.
 */
void DirichletBoundary::assignPulsatileBCs(const int t, const double dt, const double T, const int pulseBeginItr,
                                           const int nDofsGlobal)
{
  double timePhase = (t - pulseBeginItr) * dt;
  double pulse = 0.25 * sin((2e0 * PI / T) * timePhase) + 1.0;
  for(int id = 0; id < nDofsGlobal; id++) {
    if(dirichletBCsValueNew[id] > 0)
      dirichletBCsValueNew[id] = dirichletBCsValueNewInit[id] * pulse;
  }
}

/*******************************************
 * @brief Set dirichlet boundary conditions.
 */
void DirichletBoundary::applyDirichletBCs(Cell &cell, PetscSolver &petsc)
{
  std::vector<int> vecTmp;

  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    if(cell(ic).subId == mpi.myId) {
      int nDofsInCell = cell(ic).dofsMap.size();
      PetscScalar FlocalTmp[nDofsInCell];
      PetscScalar KlocalTmp[nDofsInCell * nDofsInCell];

      for(int i = 0; i < nDofsInCell; i++)
        FlocalTmp[i] = 0e0;
      for(int i = 0; i < nDofsInCell * nDofsInCell; i++)
        KlocalTmp[i] = 0e0;

      vecTmp = cell(ic).dofsMap;
      for(int i = 0; i < nDofsInCell; i++) {
        if(cell(ic).dofsBCsMap[i] == -1) {
          KlocalTmp[i + i * nDofsInCell] = 1;
          FlocalTmp[i] = dirichletBCsValueNew[cell(ic).dofsMap[i]];
        }
      }
      MatSetValues(petsc.mtx, nDofsInCell, &vecTmp[0], nDofsInCell, &vecTmp[0], KlocalTmp, INSERT_VALUES);
      VecSetValues(petsc.rhsVec, nDofsInCell, &vecTmp[0], FlocalTmp, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(petsc.mtx, MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(petsc.mtx, MAT_FLUSH_ASSEMBLY);

  VecAssemblyBegin(petsc.rhsVec);
  VecAssemblyEnd(petsc.rhsVec);
}

/*******************************************************
 * @brief Set dirichlet boundary conditions for adjoint.
 */
void DirichletBoundary::applyDirichletBCsAdjoint(Cell &cell, PetscSolver &petsc)
{
  std::vector<int> vecTmp;

  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    if(cell(ic).subId == mpi.myId) {
      int nDofsInCell = cell(ic).dofsMapWall.size();
      PetscScalar FlocalTmp[nDofsInCell];
      PetscScalar KlocalTmp[nDofsInCell * nDofsInCell];

      for(int i = 0; i < nDofsInCell; i++)
        FlocalTmp[i] = 0e0;
      for(int i = 0; i < nDofsInCell * nDofsInCell; i++)
        KlocalTmp[i] = 0e0;

      vecTmp = cell(ic).dofsMapWall;
      for(int i = 0; i < nDofsInCell; i++) {
        if(cell(ic).dofsBCsMapWall[i] == -1) {
          KlocalTmp[i + i * nDofsInCell] = 1;
          FlocalTmp[i] = dirichletBCsValueNew[cell(ic).dofsMapWall[i]];
        }
      }
      MatSetValues(petsc.mtx, nDofsInCell, &vecTmp[0], nDofsInCell, &vecTmp[0], KlocalTmp, INSERT_VALUES);
      VecSetValues(petsc.rhsVec, nDofsInCell, &vecTmp[0], FlocalTmp, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(petsc.mtx, MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(petsc.mtx, MAT_FLUSH_ASSEMBLY);

  VecAssemblyBegin(petsc.rhsVec);
  VecAssemblyEnd(petsc.rhsVec);
}

// add

/***************************************
 * @brief Initialize dirichlet boundary.
 */
void Dirichlet::initialize(Config &conf)
{
  velocitySet = std::move(conf.vDirichlet);
  pressureSet = std::move(conf.pDirichlet);
}

void Dirichlet::getNewArray(std::vector<int> mapNew)
{
  for(const auto &[idx, values] : velocitySet) {
    int newIdx = mapNew[idx];
    velocitySetNew[newIdx] = values;
  }
  for(const auto &[idx, value] : pressureSet) {
    int newIdx = mapNew[idx];
    pressureSetNew[newIdx] = value;
  }
}

void Dirichlet::setValuesZero(int n)
{
  values.resize(n);
  initialValues.resize(n);
  values.fillZero();
  initialValues.fillZero();
}

void Dirichlet::assignBCs(Node &node, const int t)
{
  for(const auto &[idx, vec] : velocitySetNew) {
    for(int d = 0; d < 3; d++) {
      int dof = node.getDof(idx, d);
      values(dof) = vec[d];
      initialValues(dof) = vec[d];
    }
  }
  for(const auto &[idx, value] : pressureSetNew) {
    int dof = node.getDof(idx, 4);
    values(dof) = value;
    initialValues(dof) = value;
  }
}

/*******************************************************
 * @brief Assign pulsatile dirichlet boundary condition.
 */
void Dirichlet::assignPulsatileBCs(const double pulse, const int nDofsGlobal)
{
  for(int id = 0; id < nDofsGlobal; id++) {
    if(values(id) > 0) {
      values(id) = initialValues(id) * pulse;
    }
  }
}

/*******************************************
 * @brief Set dirichlet boundary conditions.
 */
void Dirichlet::applyBCs(Cell &cell, PetscSolver &petsc)
{
  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    if(cell(ic).subId == mpi.myId) {
      std::vector<int> vec = cell(ic).dofsMap;
      int nDofs = vec.size();

      MatrixXd Klocal(nDofs, nDofs);
      VectorXd Flocal(nDofs);

      Klocal.setZero();
      Flocal.setZero();

      for(int i = 0; i < nDofs; i++) {
        int dof = cell(ic).dofsMap[i];
        if(cell(ic).dofsBCsMap[i] == -1) {
          Klocal(i, i) = 1;
          Flocal(i) = values(dof);
        }
      }
      petsc.setMatValue(vec, vec, Klocal);
      petsc.setVecValue(vec, Flocal);
    }
  }

  petsc.flashAssembly();
}

void Dirichlet::updateValues(Array3D<double> &X, const int t)
{
  for(auto &[idx, vec] : velocitySet) {
    for(int d = 0; d < 3; d++) {
      vec[d] = X(t, idx, d);
    }
  }
}