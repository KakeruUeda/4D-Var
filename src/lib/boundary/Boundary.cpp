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
  for(int in = 0; in < cb.CBNodeMap.size(); in++) {
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
  for(auto &entry : velocitySet) {
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
    //if(vec[0] == 0 && vec[1] == 0 && vec[2] == 0) {
    //  continue;
    //}
    for(int d = 0; d < 3; d++) {
      vec[d] = X(t, idx, d);
    }
  }
}

constexpr Dirichlet::Velocity MRI_flowRate_inlet_cycle1[32] = {
    {0 * 0.02947812, 6.43433e-05},
    {1 * 0.02947812, 7.81532e-05},
    {2 * 0.02947812, 9.92647e-05},
    {3 * 0.02947812, 0.000122226},
    {4 * 0.02947812, 0.000146049},
    {5 * 0.02947812, 0.000159925},
    {6 * 0.02947812, 0.000165452},
    {7 * 0.02947812, 0.000168617},
    {8 * 0.02947812, 0.000168412},
    {9 * 0.02947812, 0.000163581},
    {10 * 0.02947812, 0.00015419},
    {11 * 0.02947812, 0.000139544},
    {12 * 0.02947812, 0.000119993},
    {13 * 0.02947812, 0.000102198},
    {14 * 0.02947812, 8.84328e-05},
    {15 * 0.02947812, 6.77136e-05},
    {16 * 0.02947812, 4.18576e-05},
    {17 * 0.02947812, 3.41221e-05},
    {18 * 0.02947812, 4.49835e-05},
    {19 * 0.02947812, 4.48696e-05},
    {20 * 0.02947812, 3.86639e-05},
    {21 * 0.02947812, 2.83121e-05},
    {22 * 0.02947812, 1.23957e-05},
    {23 * 0.02947812, 6.05345e-07},
    {24 * 0.02947812, -1.39074e-06},
    {25 * 0.02947812, -3.35728e-07},
    {26 * 0.02947812, 5.29686e-07},
    {27 * 0.02947812, 1.61686e-06},
    {28 * 0.02947812, 6.05608e-06},
    {29 * 0.02947812, 1.04744e-05},
    {30 * 0.02947812, 1.36829e-05},
    {31 * 0.02947812, 2.07897e-05}
};

constexpr Dirichlet::Velocity MRI_flowRate_inlet_cycle2[32] = {
    {32 * 0.02947812, 6.43433e-05},
    {33 * 0.02947812, 7.81532e-05},
    {34 * 0.02947812, 9.92647e-05},
    {35 * 0.02947812, 0.000122226},
    {36 * 0.02947812, 0.000146049},
    {37 * 0.02947812, 0.000159925},
    {38 * 0.02947812, 0.000165452},
    {39 * 0.02947812, 0.000168617},
    {40 * 0.02947812, 0.000168412},
    {41 * 0.02947812, 0.000163581},
    {42 * 0.02947812, 0.00015419},
    {43 * 0.02947812, 0.000139544},
    {44 * 0.02947812, 0.000119993},
    {45 * 0.02947812, 0.000102198},
    {46 * 0.02947812, 8.84328e-05},
    {47 * 0.02947812, 6.77136e-05},
    {48 * 0.02947812, 4.18576e-05},
    {49 * 0.02947812, 3.41221e-05},
    {50 * 0.02947812, 4.49835e-05},
    {51 * 0.02947812, 4.48696e-05},
    {52 * 0.02947812, 3.86639e-05},
    {53 * 0.02947812, 2.83121e-05},
    {54 * 0.02947812, 1.23957e-05},
    {55 * 0.02947812, 6.05345e-07},
    {56 * 0.02947812, -1.39074e-06},
    {57 * 0.02947812, -3.35728e-07},
    {58 * 0.02947812, 5.29686e-07},
    {59 * 0.02947812, 1.61686e-06},
    {60 * 0.02947812, 6.05608e-06},
    {61 * 0.02947812, 1.04744e-05},
    {62 * 0.02947812, 1.36829e-05},
    {63 * 0.02947812, 2.07897e-05}
};

constexpr Dirichlet::Velocity MRI_flowRate_inlet_cycle3[32] = {
    {64 * 0.02947812, 6.43433e-05},
    {65 * 0.02947812, 7.81532e-05},
    {66 * 0.02947812, 9.92647e-05},
    {67 * 0.02947812, 0.000122226},
    {68 * 0.02947812, 0.000146049},
    {69 * 0.02947812, 0.000159925},
    {70 * 0.02947812, 0.000165452},
    {71 * 0.02947812, 0.000168617},
    {72 * 0.02947812, 0.000168412},
    {73 * 0.02947812, 0.000163581},
    {74 * 0.02947812, 0.00015419},
    {75 * 0.02947812, 0.000139544},
    {76 * 0.02947812, 0.000119993},
    {77 * 0.02947812, 0.000102198},
    {78 * 0.02947812, 8.84328e-05},
    {79 * 0.02947812, 6.77136e-05},
    {80 * 0.02947812, 4.18576e-05},
    {81 * 0.02947812, 3.41221e-05},
    {82 * 0.02947812, 4.49835e-05},
    {83 * 0.02947812, 4.48696e-05},
    {84 * 0.02947812, 3.86639e-05},
    {85 * 0.02947812, 2.83121e-05},
    {86 * 0.02947812, 1.23957e-05},
    {87 * 0.02947812, 6.05345e-07},
    {88 * 0.02947812, -1.39074e-06},
    {89 * 0.02947812, -3.35728e-07},
    {90 * 0.02947812, 5.29686e-07},
    {91 * 0.02947812, 1.61686e-06},
    {92 * 0.02947812, 6.05608e-06},
    {93 * 0.02947812, 1.04744e-05},
    {94 * 0.02947812, 1.36829e-05},
    {95 * 0.02947812, 2.07897e-05}
};

constexpr Dirichlet::Velocity MRI_flowRate_inlet_cycle4[32] = {
    {96 * 0.02947812, 6.43433e-05},
    {97 * 0.02947812, 7.81532e-05},
    {98 * 0.02947812, 9.92647e-05},
    {99 * 0.02947812, 0.000122226},
    {100 * 0.02947812, 0.000146049},
    {101 * 0.02947812, 0.000159925},
    {102 * 0.02947812, 0.000165452},
    {103 * 0.02947812, 0.000168617},
    {104 * 0.02947812, 0.000168412},
    {105 * 0.02947812, 0.000163581},
    {106 * 0.02947812, 0.00015419},
    {107 * 0.02947812, 0.000139544},
    {108 * 0.02947812, 0.000119993},
    {109 * 0.02947812, 0.000102198},
    {110 * 0.02947812, 8.84328e-05},
    {111 * 0.02947812, 6.77136e-05},
    {112 * 0.02947812, 4.18576e-05},
    {113 * 0.02947812, 3.41221e-05},
    {114 * 0.02947812, 4.49835e-05},
    {115 * 0.02947812, 4.48696e-05},
    {116 * 0.02947812, 3.86639e-05},
    {117 * 0.02947812, 2.83121e-05},
    {118 * 0.02947812, 1.23957e-05},
    {119 * 0.02947812, 6.05345e-07},
    {120 * 0.02947812, -1.39074e-06},
    {121 * 0.02947812, -3.35728e-07},
    {122 * 0.02947812, 5.29686e-07},
    {123 * 0.02947812, 1.61686e-06},
    {124 * 0.02947812, 6.05608e-06},
    {125 * 0.02947812, 1.04744e-05},
    {126 * 0.02947812, 1.36829e-05},
    {127 * 0.02947812, 2.07897e-05}
};

double Dirichlet::comp_pulse(double timeNow)
{
  std::vector<double> x, y;

  for(int t = 0; t < 32; t++) {
    x.push_back(MRI_flowRate_inlet_cycle1[t].time);
    y.push_back(MRI_flowRate_inlet_cycle1[t].value);
  }
  for(int t = 0; t < 32; t++) {
    x.push_back(MRI_flowRate_inlet_cycle2[t].time);
    y.push_back(MRI_flowRate_inlet_cycle2[t].value);
  }
  for(int t = 0; t < 32; t++) {
    x.push_back(MRI_flowRate_inlet_cycle3[t].time);
    y.push_back(MRI_flowRate_inlet_cycle3[t].value);
  }
  for(int t = 0; t < 32; t++) {
    x.push_back(MRI_flowRate_inlet_cycle4[t].time);
    y.push_back(MRI_flowRate_inlet_cycle4[t].value);
  }

  vector<Spline::Coefficients> cf = Spline::compCoefficients(x, y);
  double pulse = Spline::evaluate(cf, timeNow);

  // double maxTime = MRI_flowRate_inlet_cycle3[31].time;
  // if(timeNow > maxTime) {
  //   if(mpi.myId == 0) {
  //     std::cout << "SimTime has completed three cycles. Exiting..." << std::endl;
  //   }
  //   MPI_Finalize();
  //   std::exit(0);
  // }

  // static double maxTime = MRI_flowRate_inlet[31].time;
  // static bool maxTimeUpdated = false;

  // if(timeNow > maxTime && !maxTimeUpdated) {
  //   maxTime = timeNow;
  //   maxTimeUpdated = true;
  // }

  // if(timeNow >= 3 * maxTime) {
  //   if(mpi.myId == 0) {
  //     std::cout << "SimTime has completed three cycles. Exiting..." << std::endl;
  //   }
  //   MPI_Finalize();
  //   std::exit(0);
  // }

  // double adjustedTime = std::fmod(timeNow, maxTime);
  // double pulse = Spline::evaluate(cf, adjustedTime);

  return pulse;
}

/**
 * @brief tmp.
 */
double Dirichlet::comp_pulse2(double timeNow, std::vector<std::array<double, 2>> &velArr)
{
  std::vector<double> x, y;

  for(const auto &arr : velArr) {
    x.push_back(arr[0]);
    y.push_back(arr[1]);
  }

  vector<Spline::Coefficients> cf = Spline::compCoefficients(x, y);

  double pulse = Spline::evaluate(cf, timeNow);
  return pulse;
}

/**
 * @brief Initialization.
 */
void ControlBoundary::initialize(Config &conf)
{
  CBNodeMap = std::move(conf.CBNodeMap);
  CBEdgeNodeMap = std::move(conf.CBEdgeNodeMap);
  CBCellMap = std::move(conf.CBCellMap);
  CBNodeMapInCell = std::move(conf.CBNodeMapInCell);
}
