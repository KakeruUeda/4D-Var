/**
 * @file OptimalCondition.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "InverseProblem.h"

/**
 * @brief Compute gradient of control variables.
 */
void InverseProblem::compOptimalCondition()
{
  int nc = inletCB.CBNodeMapInCell[0].size();

  mt2d.nNodesInCell = nc;
  mt2d.N.allocate(nc);
  mt2d.xCurrent.allocate(nc, 2);
  mt2d.dNdr.allocate(nc, 2);
  mt2d.dNdx.allocate(nc, 2);

  mt3d.nNodesInCell = main.grid.cell.nNodesInCell;
  mt3d.N.allocate(main.grid.cell.nNodesInCell);
  mt3d.xCurrent.allocate(main.grid.cell.nNodesInCell, 3);
  mt3d.dNdr.allocate(main.grid.cell.nNodesInCell, 3);
  mt3d.dNdx.allocate(main.grid.cell.nNodesInCell, 3);

  std::vector<std::vector<double>> value1, value2, value3, value4, value5;

  VecTool::resize(value1, main.grid.cell.nNodesInCell, 3);
  VecTool::resize(value2, main.grid.cell.nNodesInCell, 3);
  VecTool::resize(value3, main.grid.cell.nNodesInCell, 3);
  VecTool::resize(value4, main.grid.cell.nNodesInCell, 3);
  VecTool::resize(value5, main.grid.cell.nNodesInCell, 3);

  gradX.fillZero();
  gradX0.fillZero();

  Gauss gauss(2);

  // Optimal condition for inlet boundary condition
  for(int t = 0; t < main.timeMax; t++) {
    for(int ic = 0; ic < inletCB.CBCellMap.size(); ic++) {

      for(int p = 0; p < nc; p++) {
        int in = inletCB.CBNodeMapInCell[ic][p];
        for(int d = 0; d < 2; d++) {
          mt2d.xCurrent(p, d) = main.grid.node.x[in][planeDir[d]];
        }
      }

      for(auto &row : value1) {
        std::fill(row.begin(), row.end(), 0e0);
      }
      for(auto &row : value2) {
        std::fill(row.begin(), row.end(), 0e0);
      }
      for(auto &row : value3) {
        std::fill(row.begin(), row.end(), 0e0);
      }
      for(auto &row : value4) {
        std::fill(row.begin(), row.end(), 0e0);
      }
      for(auto &row : value5) {
        std::fill(row.begin(), row.end(), 0e0);
      }

      for(int i1 = 0; i1 < 2; i1++) {
        for(int i2 = 0; i2 < 2; i2++) {
          mt2d.setShapesInGauss(gauss, i1, i2);
          mt2d.setFactorsInGauss(gauss, i1, i2);
          OptCondX_Term1_inGaussIntegral(value1, nc, ic, t);
          OptCondX_Term2_inGaussIntegral(value2, nc, ic, t);
          OptCondX_Term3_inGaussIntegral(value3, nc, ic, t);
          OptCondX_Term4_inGaussIntegral(value4, nc, ic, t);
          OptCondX_Term5_inGaussIntegral(value5, nc, ic, t);
        }
      }

      for(int p = 0; p < nc; p++) {
        int in = inletCB.CBNodeMapInCell[ic][p];
        for(int d = 0; d < 3; d++) {
          gradX(t, in, d) += bCF * value1[p][d];
          gradX(t, in, d) += bCF * value2[p][d];
          gradX(t, in, d) += bCF * value3[p][d];
          gradX(t, in, d) += bCF * value4[p][d];
          gradX(t, in, d) += value5[p][d];
        }
      }
    }
  }

  VecTool::resize(value1, main.grid.cell.nNodesInCell, 3);
  VecTool::resize(value2, main.grid.cell.nNodesInCell, 3);
  VecTool::resize(value3, main.grid.cell.nNodesInCell, 3);

  /* Optimal condition for initial velicity field */
  for(int ic = 0; ic < main.grid.cell.nCellsGlobal; ic++) {
    for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
      int n = main.grid.cell(ic).node[p];
      for(int d = 0; d < main.dim; d++) {
        mt3d.xCurrent(p, d) = main.grid.node.x[n][d];
      }
    }

    for(auto &row : value1) {
      std::fill(row.begin(), row.end(), 0e0);
    }
    for(auto &row : value2) {
      std::fill(row.begin(), row.end(), 0e0);
    }
    for(auto &row : value3) {
      std::fill(row.begin(), row.end(), 0e0);
    }

    for(int i1 = 0; i1 < 2; i1++) {
      for(int i2 = 0; i2 < 2; i2++) {
        for(int i3 = 0; i3 < 2; i3++) {
          mt3d.setShapesInGauss(gauss, i1, i2, i3);
          mt3d.setFactorsInGauss(gauss, i1, i2, i3);
          OptCondX0_Term1_inGaussIntegral(value1, ic);
          OptCondX0_Term2_inGaussIntegral(value2, ic);
          OptCondX0_Term3_inGaussIntegral(value3, ic);
        }
      }
    }

    for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
      int in = main.grid.cell(ic).node[p];
      for(int d = 0; d < main.dim; d++) {
        gradX0(in, d) += gCF * value1[p][d];
        gradX0(in, d) += gCF * value2[p][d];
        gradX0(in, d) += value3[p][d];
      }
    }
  }

  /*
  for(int icb = 0; icb < inletCB.CBNodeMap.size(); icb++){
    int in = inletCB.CBNodeMap[icb];
    for(int d = 0; d < 3; d++){
      gradX0(in, d) = 0e0;
    }
  }
  */
}

void InverseProblem::OptCondX_Term1_inGaussIntegral(std::vector<std::vector<double>> &value, const int nc, const int ic,
                                                    const int t)
{
  double u[3] = {0e0, 0e0, 0e0};

  for(int d = 0; d < dim; d++) {
    for(int p = 0; p < nc; p++) {
      int n = inletCB.CBNodeMapInCell[ic][p];
      u[d] += mt2d.N(p) * main.vt(t, n, d);
    }
  }

  for(int p = 0; p < nc; p++) {
    for(int d = 0; d < 3; d++) {
      value[p][d] += u[d] * mt2d.N(p) * mt2d.vol;
    }
  }
}

void InverseProblem::OptCondX_Term2_inGaussIntegral(std::vector<std::vector<double>> &value, const int nc, const int ic,
                                                    const int t)
{
  double dudx[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 2; d2++) {
      for(int p = 0; p < nc; p++) {
        int n = inletCB.CBNodeMapInCell[ic][p];
        dudx[d1][d2] += mt2d.dNdx(p, d2) * main.vt(t, n, d1);
      }
    }
  }

  for(int p = 0; p < nc; p++) {
    for(int d1 = 0; d1 < 3; d1++) {
      for(int d2 = 0; d2 < 2; d2++) {
        value[p][d1] += dudx[d1][d2] * mt2d.N(p) * mt2d.vol;
      }
    }
  }
}

void InverseProblem::OptCondX_Term3_inGaussIntegral(std::vector<std::vector<double>> &value, const int nc, const int ic,
                                                    const int t)
{
  double u[3] = {0e0, 0e0, 0e0};
  double ub[3] = {0e0, 0e0, 0e0};
  double dudt[3] = {0e0, 0e0, 0e0};

  for(int p = 0; p < nc; p++) {
    int n = inletCB.CBNodeMapInCell[ic][p];
    for(int d = 0; d < 3; d++) {
      if(t == 0) {
        u[d] += mt2d.N(p) * main.vt(t, n, d);
        ub[d] += mt2d.N(p) * main.v0(n, d);
      } else {
        u[d] += mt2d.N(p) * main.vt(t, n, d);
        ub[d] += mt2d.N(p) * main.vt(t - 1, n, d);
      }
    }
  }

  for(int d = 0; d < 3; d++) {
    dudt[d] = (u[d] - ub[d]) / main.dt;
  }

  for(int p = 0; p < nc; p++) {
    for(int d = 0; d < 3; d++) {
      value[p][d] += dudt[d] * mt2d.N(p) * mt2d.vol;
    }
  }
}

void InverseProblem::OptCondX_Term4_inGaussIntegral(std::vector<std::vector<double>> &value, const int nc, const int ic,
                                                    const int t)
{
  double dudx[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
  double dubdx[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
  double dudxdt[3][2] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 2; d2++) {
      for(int p = 0; p < nc; p++) {
        int n = inletCB.CBNodeMapInCell[ic][p];
        if(t == 0) {
          dudx[d1][d2] += mt2d.dNdx(p, d2) * main.vt(t, n, d1);
          dubdx[d1][d2] += mt2d.dNdx(p, d2) * main.v0(n, d1);
        } else {
          dudx[d1][d2] += mt2d.dNdx(p, d2) * main.vt(t, n, d1);
          dubdx[d1][d2] += mt2d.dNdx(p, d2) * main.vt(t - 1, n, d1);
        }
      }
    }
  }

  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 2; d2++) {
      dudxdt[d1][d2] = (dudx[d1][d2] - dubdx[d1][d2]) / main.dt;
    }
  }

  for(int p = 0; p < nc; p++) {
    for(int d1 = 0; d1 < 3; d1++) {
      for(int d2 = 0; d2 < 2; d2++) {
        value[p][d1] += dudxdt[d1][d2] * mt2d.N(p) * mt2d.vol;
      }
    }
  }
}

void InverseProblem::OptCondX_Term5_inGaussIntegral(std::vector<std::vector<double>> &value, const int nc, const int ic,
                                                    const int t)
{
  double lgp[3] = {0e0, 0e0, 0e0};

  for(int d = 0; d < 3; d++) {
    for(int p = 0; p < nc; p++) {
      int n = inletCB.CBNodeMapInCell[ic][p];
      lgp[d] -= mt2d.N(p) * adjoint.lt(t, n, d);
    }
  }

  for(int p = 0; p < nc; p++) {
    for(int d = 0; d < 3; d++) {
      value[p][d] += lgp[d] * mt2d.N(p) * mt2d.vol;
    }
  }
}

/**
 * @brief Compute value for term1 in optimal condition
 *        on gauss integral point.
 */
void InverseProblem::OptCondX0_Term1_inGaussIntegral(std::vector<std::vector<double>> &value, const int ic)
{
  double u0[3] = {0e0, 0e0, 0e0};

  for(int d = 0; d < 3; d++) {
    for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
      int n = main.grid.cell(ic).node[p];
      u0[d] += mt3d.N(p) * main.v0(n, d);
    }
  }

  for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
    for(int d = 0; d < 3; d++) {
      value[p][d] += u0[d] * mt3d.N(p) * mt3d.vol;
    }
  }
}

/**
 * @brief Compute value for term2 in optimal condition
 *        on gauss integral point.
 */
void InverseProblem::OptCondX0_Term2_inGaussIntegral(std::vector<std::vector<double>> &value, const int ic)
{
  double du0dx[3][3] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 3; d2++) {
      for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
        int n = main.grid.cell(ic).node[p];
        du0dx[d1][d2] += mt3d.dNdx(p, d2) * main.v0(n, d1);
      }
    }
  }

  for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
    for(int d1 = 0; d1 < 3; d1++) {
      for(int d2 = 0; d2 < 3; d2++) {
        value[p][d1] += du0dx[d1][d2] * mt3d.N(p) * mt3d.vol;
      }
    }
  }
}

/**
 * @brief Compute value for term3 in optimal condition
 *        on gauss integral point.
 */
void InverseProblem::OptCondX0_Term3_inGaussIntegral(std::vector<std::vector<double>> &value, const int ic)
{
  setValue(ic);

  double he = main.comp_he(mt3d.xCurrent);
  double f = main.comp_f(main.grid.cell(ic).phi);

  main.tau = main.comp_tau(adjoint.adw, he);
  //main.tau = main.comp_tau2(mt3d.dNdx, adjoint.adw, mt3d.nNodesInCell);

  for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
    // Mass term
    value[p][0] -= main.rho * mt3d.N(p) * adjoint.wk1[0] / main.dt * mt3d.vol;
    value[p][1] -= main.rho * mt3d.N(p) * adjoint.wk1[1] / main.dt * mt3d.vol;
    value[p][2] -= main.rho * mt3d.N(p) * adjoint.wk1[2] / main.dt * mt3d.vol;

    // Diffusion term
    for(int d = 0; d < 3; d++) {
      value[p][0] += 5e-1 * main.mu * mt3d.dNdx(p, d) * adjoint.dwk1dx[0][d] * mt3d.vol;
      value[p][1] += 5e-1 * main.mu * mt3d.dNdx(p, d) * adjoint.dwk1dx[1][d] * mt3d.vol;
      value[p][2] += 5e-1 * main.mu * mt3d.dNdx(p, d) * adjoint.dwk1dx[2][d] * mt3d.vol;
    }

    // Advection term
    value[p][0] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvk1dx[0][0] * adjoint.wk1[0] * mt3d.vol;
    value[p][0] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvkdx[0][0] * adjoint.wk1[0] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][0] += 0.5 * main.rho * mt3d.dNdx(p, d) * adjoint.advk2[d] * adjoint.wk1[0] * mt3d.vol;
    }
    value[p][0] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvk1dx[1][0] * adjoint.wk1[1] * mt3d.vol;
    value[p][0] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvkdx[1][0] * adjoint.wk1[1] * mt3d.vol;
    value[p][0] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvk1dx[2][0] * adjoint.wk1[2] * mt3d.vol;
    value[p][0] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvkdx[2][0] * adjoint.wk1[2] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][0] += 0.5 * (-0.5) * main.rho * mt3d.N(p) * adjoint.dvk2dx[d][0] * adjoint.wk2[d] * mt3d.vol;
      value[p][0] += 0.5 * (-0.5) * main.rho * mt3d.N(p) * adjoint.dvk1dx[d][0] * adjoint.wk2[d] * mt3d.vol;
    }

    value[p][1] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvk1dx[1][1] * adjoint.wk1[1] * mt3d.vol;
    value[p][1] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvkdx[1][1] * adjoint.wk1[1] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][1] += 0.5 * main.rho * mt3d.dNdx(p, d) * adjoint.advk2[d] * adjoint.wk1[1] * mt3d.vol;
    }
    value[p][1] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvk1dx[0][1] * adjoint.wk1[0] * mt3d.vol;
    value[p][1] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvkdx[0][1] * adjoint.wk1[0] * mt3d.vol;
    value[p][1] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvk1dx[2][1] * adjoint.wk1[2] * mt3d.vol;
    value[p][1] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvkdx[2][1] * adjoint.wk1[2] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][1] += 0.5 * (-0.5) * main.rho * mt3d.N(p) * adjoint.dvk2dx[d][1] * adjoint.wk2[d] * mt3d.vol;
      value[p][1] += 0.5 * (-0.5) * main.rho * mt3d.N(p) * adjoint.dvk1dx[d][1] * adjoint.wk2[d] * mt3d.vol;
    }

    value[p][2] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvk1dx[2][2] * adjoint.wk1[2] * mt3d.vol;
    value[p][2] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvkdx[2][2] * adjoint.wk1[2] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][2] += 0.5 * main.rho * mt3d.dNdx(p, d) * adjoint.advk2[d] * adjoint.wk1[2] * mt3d.vol;
    }
    value[p][2] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvk1dx[0][2] * adjoint.wk1[0] * mt3d.vol;
    value[p][2] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvkdx[0][2] * adjoint.wk1[0] * mt3d.vol;
    value[p][2] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvk1dx[1][2] * adjoint.wk1[1] * mt3d.vol;
    value[p][2] += 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dvkdx[1][2] * adjoint.wk1[1] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][2] += 0.5 * (-0.5) * main.rho * mt3d.N(p) * adjoint.dvk2dx[d][2] * adjoint.wk2[d] * mt3d.vol;
      value[p][2] += 0.5 * (-0.5) * main.rho * mt3d.N(p) * adjoint.dvk1dx[d][2] * adjoint.wk2[d] * mt3d.vol;
    }

    // Darcy term
    value[p][0] += 5e-1 * f / main.rho * mt3d.N(p) * adjoint.wk1[0] * mt3d.vol;
    value[p][1] += 5e-1 * f / main.rho * mt3d.N(p) * adjoint.wk1[1] * mt3d.vol;
    value[p][2] += 5e-1 * f / main.rho * mt3d.N(p) * adjoint.wk1[2] * mt3d.vol;

    // SUPG mass term
    for(int d = 0; d < main.dim; d++) {
      value[p][0] += main.tau * 1.5 * main.rho * mt3d.N(p) * (adjoint.vk1[d] - adjoint.vk[d]) * adjoint.dwk1dx[d][0] / main.dt  * mt3d.vol;
      value[p][0] -= main.tau * 0.5 * main.rho * mt3d.N(p) * adjoint.dwk2dx[d][0] * (adjoint.vk2[d] - adjoint.vk1[d]) / main.dt * mt3d.vol;
    }
    value[p][0] -= main.tau * main.rho * mt3d.N(p) * adjoint.advk2[0] * adjoint.dwk1dx[0][0] / main.dt * mt3d.vol;

    for(int d = 0; d < main.dim; d++) {
      value[p][1] += main.tau * 1.5 * main.rho * mt3d.N(p) * (adjoint.vk1[d] - adjoint.vk[d]) * adjoint.dwk1dx[d][1] / main.dt  * mt3d.vol;
      value[p][1] -= main.tau * 0.5 * main.rho * mt3d.N(p) * adjoint.dwk2dx[d][1] * (adjoint.vk2[d] - adjoint.vk1[d]) / main.dt * mt3d.vol;
    }
    value[p][1] -= main.tau * main.rho * mt3d.N(p) * adjoint.advk2[1] * adjoint.dwk1dx[1][1] / main.dt * mt3d.vol;

    for(int d = 0; d < main.dim; d++) {
      value[p][2] += main.tau * 1.5 * main.rho * mt3d.N(p) * (adjoint.vk1[d] - adjoint.vk[d]) * adjoint.dwk1dx[d][2] / main.dt  * mt3d.vol;
      value[p][2] -= main.tau * 0.5 * main.rho * mt3d.N(p) * adjoint.dwk2dx[d][2] * (adjoint.vk2[d] - adjoint.vk1[d]) / main.dt * mt3d.vol;
    }
    value[p][2] -= main.tau * main.rho * mt3d.N(p) * adjoint.advk2[2] * adjoint.dwk1dx[2][2] / main.dt * mt3d.vol;
   
    // SUPG advection term
    std::vector<double> frontAdv2, frontAdv3;
    VecTool::resize(frontAdv2, main.dim);
    VecTool::resize(frontAdv3, main.dim);

    for(int d1 = 0; d1 < main.dim; d1++) {
      for(int d2 = 0; d2 < main.dim; d2++) {
        frontAdv2[d1] += adjoint.advk2[d2] * adjoint.dwk1dx[d1][d2];
        frontAdv3[d1] += adjoint.advk3[d2] * adjoint.dwk2dx[d1][d2];
      }
    }

    value[p][0] += main.tau * main.rho * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[0][0] * mt3d.vol;
    value[p][0] += main.tau * main.rho * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[0][0] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][0] += main.tau * main.rho * frontAdv2[0] * 0.5 * mt3d.dNdx(p, d) * adjoint.advk2[d] * mt3d.vol;
    }
    value[p][0] += main.tau * main.rho * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[1][0] * mt3d.vol;
    value[p][0] += main.tau * main.rho * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[1][0] * mt3d.vol;
    value[p][0] += main.tau * main.rho * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[2][0] * mt3d.vol;
    value[p][0] += main.tau * main.rho * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[2][0] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][0] += main.tau * main.rho * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk2dx[d][0] * mt3d.vol;
      value[p][0] += main.tau * main.rho * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk1dx[d][0] * mt3d.vol;
    }

    value[p][1] += main.tau * main.rho * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[1][1] * mt3d.vol;
    value[p][1] += main.tau * main.rho * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[1][1] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][1] += main.tau * main.rho * frontAdv2[1] * 0.5 * mt3d.dNdx(p, d) * adjoint.advk2[d] * mt3d.vol;
    }
    value[p][1] += main.tau * main.rho * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[0][1] * mt3d.vol;
    value[p][1] += main.tau * main.rho * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[0][1] * mt3d.vol;
    value[p][1] += main.tau * main.rho * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[2][1] * mt3d.vol;
    value[p][1] += main.tau * main.rho * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[2][1] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][1] += main.tau * main.rho * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk2dx[d][1] * mt3d.vol;
      value[p][1] += main.tau * main.rho * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk1dx[d][1] * mt3d.vol;
    }

    value[p][2] += main.tau * main.rho * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[2][2] * mt3d.vol;
    value[p][2] += main.tau * main.rho * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[2][2] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][2] += main.tau * main.rho * frontAdv2[2] * 0.5 * mt3d.dNdx(p, d) * adjoint.advk2[d] * mt3d.vol;
    }
    value[p][2] += main.tau * main.rho * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[0][2] * mt3d.vol;
    value[p][2] += main.tau * main.rho * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[0][2] * mt3d.vol;
    value[p][2] += main.tau * main.rho * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[1][2] * mt3d.vol;
    value[p][2] += main.tau * main.rho * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[1][2] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][2] += main.tau * main.rho * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk2dx[d][2] * mt3d.vol;
      value[p][2] += main.tau * main.rho * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk1dx[d][2] * mt3d.vol;
    }

    std::vector<double> backAdv2L, backAdv3L;
    VecTool::resize(backAdv2L, main.dim);
    VecTool::resize(backAdv3L, main.dim);

    for(int d1 = 0; d1 < main.dim; d1++) {
      for(int d2 = 0; d2 < main.dim; d2++) {
        backAdv2L[d1] += adjoint.advk2[d2] * adjoint.dvk1dx[d1][d2];
        backAdv3L[d1] += adjoint.advk3[d2] * adjoint.dvk2dx[d1][d2];
      }
    }

    std::vector<double> backAdv2R, backAdv3R;
    VecTool::resize(backAdv2R, main.dim);
    VecTool::resize(backAdv3R, main.dim);

    for(int d1 = 0; d1 < main.dim; d1++) {
      for(int d2 = 0; d2 < main.dim; d2++) {
        backAdv2R[d1] += adjoint.advk2[d2] * adjoint.dvkdx[d1][d2];
        backAdv3R[d1] += adjoint.advk3[d2] * adjoint.dvk1dx[d1][d2];
      }
    }

    for(int d = 0; d < main.dim; d++) {
      value[p][0] -= main.tau * 0.5 * 0.5 * main.rho * mt3d.N(p) * adjoint.dwk2dx[d][0] * (backAdv3L[0] + backAdv3R[0]) * mt3d.vol;
      value[p][0] += main.tau * 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dwk1dx[d][0] * (backAdv2L[0] + backAdv2R[0]) * mt3d.vol;
    }

    for(int d = 0; d < main.dim; d++) {
      value[p][1] -= main.tau * 0.5 * 0.5 * main.rho * mt3d.N(p) * adjoint.dwk2dx[d][1] * (backAdv3L[1] + backAdv3R[1]) * mt3d.vol;
      value[p][1] += main.tau * 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dwk1dx[d][1] * (backAdv2L[1] + backAdv2R[1]) * mt3d.vol;
    }

    for(int d = 0; d < main.dim; d++) {
      value[p][2] -= main.tau * 0.5 * 0.5 * main.rho * mt3d.N(p) * adjoint.dwk2dx[d][2] * (backAdv3L[2] + backAdv3R[2]) * mt3d.vol;
      value[p][2] += main.tau * 0.5 * 1.5 * main.rho * mt3d.N(p) * adjoint.dwk1dx[d][2] * (backAdv2L[2] + backAdv2R[2]) * mt3d.vol;
    }

    // SUPG pressure term
    for(int d = 0; d < 3; d++) {
      value[p][0] += main.tau * 1.5 * mt3d.N(p) * adjoint.dpk1dx[d] * adjoint.dwk1dx[d][0] * mt3d.vol;
      value[p][0] -= main.tau * 0.5 * mt3d.N(p) * adjoint.dpk2dx[d] * adjoint.dwk2dx[d][0] * mt3d.vol;
    }
    for(int d = 0; d < 3; d++) {
      value[p][1] += main.tau * 1.5 * mt3d.N(p) * adjoint.dpk1dx[d] * adjoint.dwk1dx[d][1] * mt3d.vol;
      value[p][1] -= main.tau * 0.5 * mt3d.N(p) * adjoint.dpk2dx[d] * adjoint.dwk2dx[d][1] * mt3d.vol;
    }
    for(int d = 0; d < 3; d++) {
      value[p][2] += main.tau * 1.5 * mt3d.N(p) * adjoint.dpk1dx[d] * adjoint.dwk1dx[d][2] * mt3d.vol;
      value[p][2] -= main.tau * 0.5 * mt3d.N(p) * adjoint.dpk2dx[d] * adjoint.dwk2dx[d][2] * mt3d.vol;
    }

    // PSPG mass term
    value[p][0] -= main.tau * mt3d.N(p) * adjoint.dqk1dx[0] / main.dt * mt3d.vol;
    value[p][1] -= main.tau * mt3d.N(p) * adjoint.dqk1dx[1] / main.dt * mt3d.vol;
    value[p][2] -= main.tau * mt3d.N(p) * adjoint.dqk1dx[2] / main.dt * mt3d.vol;

    // PSPG advection term
    value[p][0] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[0][0] * adjoint.dqk1dx[0] * mt3d.vol;
    value[p][0] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[0][0] * adjoint.dqk1dx[0] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][0] += main.tau * 0.5 * mt3d.dNdx(p, d) * adjoint.advk2[d] * adjoint.dqk1dx[0] * mt3d.vol;
    }
    value[p][0] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[1][0] * adjoint.dqk1dx[1] * mt3d.vol;
    value[p][0] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[1][0] * adjoint.dqk1dx[1] * mt3d.vol;
    value[p][0] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[2][0] * adjoint.dqk1dx[2] * mt3d.vol;
    value[p][0] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[2][0] * adjoint.dqk1dx[2] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][0] += main.tau * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk2dx[d][0] * adjoint.dqk2dx[d] * mt3d.vol;
      value[p][0] += main.tau * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk1dx[d][0] * adjoint.dqk2dx[d] * mt3d.vol;
    }

    value[p][1] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[1][1] * adjoint.dqk1dx[1] * mt3d.vol;
    value[p][1] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[1][1] * adjoint.dqk1dx[1] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][1] += main.tau * 0.5 * mt3d.dNdx(p, d) * adjoint.advk2[d] * adjoint.dqk1dx[1] * mt3d.vol;
    }
    value[p][1] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[0][1] * adjoint.dqk1dx[0] * mt3d.vol;
    value[p][1] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[0][1] * adjoint.dqk1dx[0] * mt3d.vol;
    value[p][1] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[2][1] * adjoint.dqk1dx[2] * mt3d.vol;
    value[p][1] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[2][1] * adjoint.dqk1dx[2] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][1] += main.tau * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk2dx[d][1] * adjoint.dqk2dx[d] * mt3d.vol;
      value[p][1] += main.tau * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk1dx[d][1] * adjoint.dqk2dx[d] * mt3d.vol;
    }

    value[p][2] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[2][2] * adjoint.dqk1dx[2] * mt3d.vol;
    value[p][2] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[2][2] * adjoint.dqk1dx[2] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][2] += main.tau * 0.5 * mt3d.dNdx(p, d) * adjoint.advk2[d] * adjoint.dqk1dx[2] * mt3d.vol;
    }
    value[p][2] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[0][2] * adjoint.dqk1dx[0] * mt3d.vol;
    value[p][2] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[0][2] * adjoint.dqk1dx[0] * mt3d.vol;
    value[p][2] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvk1dx[1][2] * adjoint.dqk1dx[1] * mt3d.vol;
    value[p][2] += main.tau * 0.5 * 1.5 * mt3d.N(p) * adjoint.dvkdx[1][2] * adjoint.dqk1dx[1] * mt3d.vol;
    for(int d = 0; d < main.dim; d++) {
      value[p][2] += main.tau * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk2dx[d][2] * adjoint.dqk2dx[d] * mt3d.vol;
      value[p][2] += main.tau * 0.5 * (-0.5) * mt3d.N(p) * adjoint.dvk1dx[d][2] * adjoint.dqk2dx[d] * mt3d.vol;
    }

  }

}

/**
 * @brief Set values needed for matrix assembly
 *        on gauss integral points.
 */
void InverseProblem::setValue(const int ic)
{
  // main var - v
  for(int d = 0; d < 3; d++) {
    adjoint.vk[d] = 0e0;
    adjoint.vk1[d] = 0e0;
    adjoint.vk2[d] = 0e0;
    for(int p = 0; p < adjoint.grid.cell.nNodesInCell; p++) {
      int n = adjoint.grid.cell(ic).node[p];
      adjoint.vk[d] += mt3d.N(p) * main.vt(0, n, d);
      adjoint.vk1[d] += mt3d.N(p) * main.vt(1, n, d);
      adjoint.vk2[d] += mt3d.N(p) * main.vt(2, n, d);
    }
  }

  // main var - dvdx
  for(int d = 0; d < 3; d++) {
    for(int e = 0; e < 3; e++) {
      adjoint.dvkdx[d][e] = 0e0;
      adjoint.dvk1dx[d][e] = 0e0;
      adjoint.dvk2dx[d][e] = 0e0;
      for(int p = 0; p < adjoint.grid.cell.nNodesInCell; p++) {
        int n = adjoint.grid.cell(ic).node[p];
        adjoint.dvkdx[d][e] += mt3d.dNdx(p, e) * main.vt(0, n, d);
        adjoint.dvk1dx[d][e] += mt3d.dNdx(p, e) * main.vt(1, n, d);
        adjoint.dvk2dx[d][e] += mt3d.dNdx(p, e) * main.vt(2, n, d);
      }
    }
  }

  // main var - dpdx
  for(int d = 0; d < 3; d++) {
    adjoint.dpkdx[d] = 0e0;
    adjoint.dpk1dx[d] = 0e0;
    adjoint.dpk2dx[d] = 0e0;
    for(int p = 0; p < adjoint.grid.cell.nNodesInCell; p++) {
      int n = adjoint.grid.cell(ic).node[p];
      adjoint.dpkdx[d] += mt3d.dNdx(p, d) * main.pt(0, n);
      adjoint.dpk1dx[d] += mt3d.dNdx(p, d) * main.pt(1, n);
      adjoint.dpk2dx[d] += mt3d.dNdx(p, d) * main.pt(2, n);
    }
  }

  // main var - adv
  for(int d = 0; d < 3; d++) {
    adjoint.advk1[d] = 0e0;
    adjoint.advk2[d] = 0e0;
    adjoint.advk3[d] = 0e0;
    for(int p = 0; p < adjoint.grid.cell.nNodesInCell; p++) {
      int n = adjoint.grid.cell(ic).node[p];
      adjoint.advk1[d] += mt3d.N(p) * main.v0(n, d);
      adjoint.advk3[d] += mt3d.N(p) * (1.5 * main.vt(1, n, d) - 0.5 * main.vt(0, n, d));
      adjoint.advk2[d] += mt3d.N(p) * (1.5 * main.vt(0, n, d) - 0.5 * main.v0(n, d));
    }
  }

  // lagrange multiplier - w
  for(int d = 0; d < 3; d++) {
    adjoint.wk1[d] = 0e0;
    adjoint.wk2[d] = 0e0;
    for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
      int n = adjoint.grid.cell(ic).node[p];
      adjoint.wk1[d] += mt3d.N(p) * adjoint.wt(0, n, d);
      adjoint.wk2[d] += mt3d.N(p) * adjoint.wt(1, n, d);
    }
  }

  for(int d = 0; d < main.dim; d++) {
    adjoint.adw[d] = 0e0;
    for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
      int n = adjoint.grid.cell(ic).node[p];
      adjoint.adw[d] += mt3d.N(p) * (1.5 * adjoint.wt(1, n, d) - 0.5 * adjoint.wt(2, n, d));
    }
  }

  // lagrange multiplier - dwdx
  for(int d = 0; d < 3; d++) {
    for(int e = 0; e < 3; e++) {
      adjoint.dwk1dx[d][e] = 0e0;
      adjoint.dwk2dx[d][e] = 0e0;
      for(int p = 0; p < main.grid.cell.nNodesInCell; p++) {
        int n = adjoint.grid.cell(ic).node[p];
        adjoint.dwk1dx[d][e] += mt3d.dNdx(p, e) * adjoint.wt(0, n, d);
        adjoint.dwk2dx[d][e] += mt3d.dNdx(p, e) * adjoint.wt(1, n, d);
      }
    }
  }

  // lagrange multiplier - dqdx
  for(int d = 0; d < 3; d++) {
    adjoint.dqk1dx[d] = 0e0;
    adjoint.dqk2dx[d] = 0e0;
    for(int p = 0; p < adjoint.grid.cell.nNodesInCell; p++) {
      int n = adjoint.grid.cell(ic).node[p];
      adjoint.dqk1dx[d] += mt3d.dNdx(p, d) * adjoint.qt(0, n);
      adjoint.dqk2dx[d] += mt3d.dNdx(p, d) * adjoint.qt(1, n);
    }
  }
}
