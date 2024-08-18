/**
 * @file Boundary.cpp
 * @author K.Ueda
 * @date May, 2024
 */

#include "Boundary.h"

/***************************************
 * @brief Initialize dirichlet boundary.
 */
void Dirichlet::initialize(Config &conf)
{
  velocitySet = conf.vDirichlet;
  pressureSet = conf.pDirichlet;
}

/*************************************************
 * @brief Erase control nodes from velocitySet.
 *        This is used for solving adjoint matrix.
 */
void Dirichlet::eraseControlNodes(Cell &cell, ControlBoundary &cb)
{
  for(int in=0; in<cb.CBNodeMap.size(); in++){
    velocitySet.erase(cb.CBNodeMap[in]);
  }
}

/***************************************
 * @brief Initialize dirichlet boundary.
 */
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
  values.allocate(n);
  initialValues.allocate(n);
  values.fillZero();
  initialValues.fillZero();
}

void Dirichlet::assignBCs(Node &node)
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

void ControlBoundary::initialize(Config &conf)
{
  // Using move is fine here
  CBNodeMap = std::move(conf.CBNodeMap);
  CBCellMap = std::move(conf.CBCellMap);
  CBNodeMapInCell = std::move(conf.CBNodeMapInCell);
}