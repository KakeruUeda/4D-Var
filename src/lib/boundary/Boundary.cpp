/**
 * @file Boundary.cpp
 * @author K.Ueda
 * @date May, 2024
 */

#include "Boundary.h"

/**
 * @brief Initialize dirichlet boundary.
 */
void Dirichlet::initialize(Config &conf)
{
  velocitySet = conf.vDirichlet;
  pressureSet = conf.pDirichlet;
}

/**
 * @brief Erase control nodes from velocitySet.
 *        This is used for solving adjoint matrix.
 */
void Dirichlet::eraseControlNodes(Cell &cell, ControlBoundary &cb)
{
  for(int in=0; in<cb.CBNodeMap.size(); in++){
    velocitySet.erase(cb.CBNodeMap[in]);
  }
}


void Dirichlet::setValuesZero(int n)
{
  values.allocate(n);
  initialValues.allocate(n);
  values.fillZero();
  initialValues.fillZero();
}

void Dirichlet::assignPulsatileBCs(const double pulse)
{
  for (auto &entry : velocitySet){
    std::vector<double> &velocities = entry.second;

    velocities[0] = velocitySetInit[entry.first][0] * pulse;
    velocities[1] = velocitySetInit[entry.first][1] * pulse;
    velocities[2] = velocitySetInit[entry.first][2] * pulse;
  }
}

/**
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

/**
 * @brief Assign dirichlet boundary conditions to dof array.
 */
void Dirichlet::assignBCs(Node &node)
{
  for(const auto &[idx, vec] : velocitySetNew) {
    for(int d = 0; d < 3; d++) {
      int dof = node.getDof(idx, d);
      values(dof) = vec[d];
    }
  }
  for(const auto &[idx, value] : pressureSetNew) {
    int dof = node.getDof(idx, 4);
    values(dof) = value;
  }
}

/**
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

/**
 * @brief Update dirichlet boundary conditions.
 */
void Dirichlet::updateValues(Array3D<double> &X, const int t)
{
  for(auto &[idx, vec] : velocitySet) {
    if(vec[0] == 0 && vec[1] == 0 && vec[2] == 0){
      continue;
    }
    for(int d = 0; d < 3; d++) {
      vec[d] = X(t, idx, d);
    }
  }
}

/**
 * @brief tmp.
 */
constexpr Dirichlet::AveragedVelocity velArray[32] = 
{
  {0 * 0.02947812, 0.11943478}, {1 * 0.02947812, 0.14436713}, 
  {2 * 0.02947812, 0.17000289}, {3 * 0.02947812, 0.20747005},
  {4 * 0.02947812, 0.2368427}, {5 * 0.02947812, 0.26329655}, 
  {6 * 0.02947812, 0.26288376}, {7 * 0.02947812, 0.25537007},
  {8 * 0.02947812, 0.25249323}, {9 * 0.02947812, 0.24208209}, 
  {10 * 0.02947812, 0.22972814}, {11 * 0.02947812, 0.19955092},
  {12 * 0.02947812, 0.17341562}, {13 * 0.02947812, 0.147503}, 
  {14 * 0.02947812, 0.12755212}, {15 * 0.02947812, 0.11137679},
  {16 * 0.02947812, 0.08759725}, {17 * 0.02947812, 0.07680825},
  {18 * 0.02947812, 0.08648109}, {19 * 0.02947812, 0.09086978},
  {20 * 0.02947812, 0.07851191}, {21 * 0.02947812, 0.06508125}, 
  {22 * 0.02947812, 0.05969055}, {23 * 0.02947812, 0.0566123},
  {24 * 0.02947812, 0.0549937}, {25 * 0.02947812, 0.05120174}, 
  {26 * 0.02947812, 0.04861119}, {27 * 0.02947812, 0.05262863},
  {28 * 0.02947812, 0.05423519}, {29 * 0.02947812, 0.04585582}, 
  {30 * 0.02947812, 0.05972958}, {31 * 0.02947812, 0.06170865}
};

/**
 * @brief tmp.
 */
double Dirichlet::comp_pulse(double timeNow)
{
  std::vector<double> x, y;

  for(int t = 0; t < 32; t++) {
    x.push_back(velArray[t].time);
    y.push_back(velArray[t].value);
  }

  vector<Spline::Coefficients> cf = Spline::compCoefficients(x, y);

  static double maxTime = velArray[31].time; 
  static bool maxTimeUpdated = false; 

  if(timeNow > maxTime && !maxTimeUpdated) {
    maxTime = timeNow; 
    maxTimeUpdated = true; 
  }

  if(timeNow >= 3 * maxTime) {
    if(mpi.myId == 0){
      std::cout << "SimTime has completed three cycles. Exiting..." << std::endl;
    }
    MPI_Finalize();
    std::exit(0); 
  }

  double adjustedTime = std::fmod(timeNow, maxTime);

  double pulse = Spline::evaluate(cf, adjustedTime);
  return pulse / velArray[0].value;
}

/**
 * @brief Initialization.
 */
void ControlBoundary::initialize(Config &conf)
{
  CBNodeMap = std::move(conf.CBNodeMap);
  CBCellMap = std::move(conf.CBCellMap);
  CBNodeMapInCell = std::move(conf.CBNodeMapInCell);
}
