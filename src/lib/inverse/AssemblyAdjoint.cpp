/**
 * @file AssemblyAdjoint.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "InverseProblem.h"

/*********************************
 * @brief Assemble adjoint system.
 */
void Adjoint::matrixAssemblyAdjoint(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t)
{
  Klocal.setZero();
  Flocal.setZero();

  mt3d.N.allocate(grid.cell.nNodesInCell);
  mt3d.xCurrent.allocate(grid.cell.nNodesInCell, 3);
  mt3d.dNdr.allocate(grid.cell.nNodesInCell, 3);
  mt3d.dNdx.allocate(grid.cell.nNodesInCell, 3);
  mt3d.K.allocate(grid.cell.nNodesInCell, grid.cell.nNodesInCell);

  for(int p = 0; p < grid.cell.nNodesInCell; p++) {
    for(int d = 0; d < main.dim; d++) {
      mt3d.xCurrent(p, d) = grid.node.x[grid.cell(ic).node[p]][d];
    }
  }

  double he = comp_he(mt3d.xCurrent);
  double f = comp_f(grid.cell(ic).phi);

  Gauss g2(2);
  for(int i1 = 0; i1 < 2; i1++) {
    for(int i2 = 0; i2 < 2; i2++) {
      for(int i3 = 0; i3 < 2; i3++) {
        setValuesInGaussIntegral(main, g2, he, i1, i2, i3, ic, t);
        for(int ii = 0; ii < grid.cell.nNodesInCell; ii++) {
          updateRowIndex(grid, ii, ic);
          for(int jj = 0; jj < grid.cell.nNodesInCell; jj++) {
            updateColumnIndex(grid, jj, ic);
            adjointGaussIntegralLHS(main, Klocal, f, ii, jj);
          }
          adjointGaussIntegralRHS(main, Flocal, f, ii);
        }
      }
    }
  }
}

/**************************************
 * @brief Set values in gauss integral.
 */
void Adjoint::setValuesInGaussIntegral(DirectProblem &main, Gauss &g2, const double he, const int i1, const int i2,
                                       const int i3, const int ic, const int t)
{
  ShapeFunction3D::C3D8_N(mt3d.N, g2.point[i1], g2.point[i2], g2.point[i3]);
  ShapeFunction3D::C3D8_dNdr(mt3d.dNdr, g2.point[i1], g2.point[i2], g2.point[i3]);

  MathTools3D::comp_dxdr(mt3d.dxdr, mt3d.dNdr, mt3d.xCurrent, grid.cell.nNodesInCell);
  MathTools3D::comp_dNdx(mt3d.dNdx, mt3d.dNdr, mt3d.dxdr, grid.cell.nNodesInCell);

  mt3d.detJ = MathTools3D::compDeterminant(mt3d.dxdr);
  mt3d.weight = g2.weight[i1] * g2.weight[i2] * g2.weight[i3];

  setValue(main, ic, t);
  tau = comp_tau(advk2, he);
}

/**********************************************************************
 * @brief Compute element sfiffness LHS matrix on gauss integral point.
 */
void Adjoint::adjointGaussIntegralLHS(DirectProblem &main, MatrixXd &Klocal, const double f, const int ii, const int jj)
{
  int n1, n2, n3;
  mt3d.vol = mt3d.detJ * mt3d.weight;

  mt3d.K(ii, jj) = 0e0;
  for(int d = 0; d < 3; d++) {
    mt3d.K(ii, jj) += mt3d.dNdx(ii, d) * mt3d.dNdx(jj, d);
  }

  // Mass term
  Klocal(IU, JU) += mt3d.N(ii) * mt3d.N(jj) / main.dt * mt3d.vol;
  Klocal(IV, JV) += mt3d.N(ii) * mt3d.N(jj) / main.dt * mt3d.vol;
  Klocal(IW, JW) += mt3d.N(ii) * mt3d.N(jj) / main.dt * mt3d.vol;

  // Diffusion term
  for(int d = 0; d < 3; d++) {
    Klocal(IU, JU) += 5e-1 * mt3d.dNdx(ii, d) * mt3d.dNdx(jj, d) / main.Re * mt3d.vol;
    Klocal(IV, JV) += 5e-1 * mt3d.dNdx(ii, d) * mt3d.dNdx(jj, d) / main.Re * mt3d.vol;
    Klocal(IW, JW) += 5e-1 * mt3d.dNdx(ii, d) * mt3d.dNdx(jj, d) / main.Re * mt3d.vol;
  }

  // Advection term
  for(int d = 0; d < main.dim; d++) {
    Klocal(IU, JU) += 5e-1 * mt3d.dNdx(ii, d) * advk1[d] * mt3d.N(jj) * mt3d.vol;
    Klocal(IV, JV) += 5e-1 * mt3d.dNdx(ii, d) * advk1[d] * mt3d.N(jj) * mt3d.vol;
    Klocal(IW, JW) += 5e-1 * mt3d.dNdx(ii, d) * advk1[d] * mt3d.N(jj) * mt3d.vol;
  }

  // Pressure term
  Klocal(IU, JP) -= mt3d.dNdx(ii, 0) * mt3d.N(jj) * mt3d.vol;
  Klocal(IV, JP) -= mt3d.dNdx(ii, 1) * mt3d.N(jj) * mt3d.vol;
  Klocal(IW, JP) -= mt3d.dNdx(ii, 2) * mt3d.N(jj) * mt3d.vol;

  // Continuity term
  Klocal(IP, JU) += mt3d.N(ii) * mt3d.dNdx(jj, 0) * mt3d.vol;
  Klocal(IP, JV) += mt3d.N(ii) * mt3d.dNdx(jj, 1) * mt3d.vol;
  Klocal(IP, JW) += mt3d.N(ii) * mt3d.dNdx(jj, 2) * mt3d.vol;

  // Darcy term
  Klocal(IU, JU) += 5e-1 * f * mt3d.N(ii) * mt3d.N(jj) * mt3d.vol;
  Klocal(IV, JV) += 5e-1 * f * mt3d.N(ii) * mt3d.N(jj) * mt3d.vol;
  Klocal(IW, JW) += 5e-1 * f * mt3d.N(ii) * mt3d.N(jj) * mt3d.vol;

  // SUPG mass term
  Klocal(IU, JU) += tau * mt3d.N(ii) * advk1[0] * mt3d.dNdx(ii, 0) / main.dt * mt3d.vol;
  Klocal(IV, JV) += tau * mt3d.N(ii) * advk1[1] * mt3d.dNdx(ii, 1) / main.dt * mt3d.vol;
  Klocal(IW, JW) += tau * mt3d.N(ii) * advk1[2] * mt3d.dNdx(ii, 2) / main.dt * mt3d.vol;

  // SUPG advection term
  for(int d1 = 0; d1 < main.dim; d1++) {
    for(int d2 = 0; d2 < main.dim; d2++) {
      Klocal(IU, JU) += tau * advk1[d1] * 5e-1 * mt3d.dNdx(ii, d2) * advk1[d2] * mt3d.dNdx(jj, d1) * mt3d.vol;
      Klocal(IV, JV) += tau * advk1[d1] * 5e-1 * mt3d.dNdx(ii, d2) * advk1[d2] * mt3d.dNdx(jj, d1) * mt3d.vol;
      Klocal(IW, JW) += tau * advk1[d1] * 5e-1 * mt3d.dNdx(ii, d2) * advk1[d2] * mt3d.dNdx(jj, d1) * mt3d.vol;
    }
  }

  // SUPG pressure term
  Klocal(IP, IU) += tau * mt3d.dNdx(ii, 0) * advk1[0] * mt3d.dNdx(jj, 0) * mt3d.vol;
  Klocal(IP, IV) += tau * mt3d.dNdx(ii, 1) * advk1[1] * mt3d.dNdx(jj, 1) * mt3d.vol;
  Klocal(IP, IW) += tau * mt3d.dNdx(ii, 2) * advk1[2] * mt3d.dNdx(jj, 2) * mt3d.vol;

  // PSPG mass term
  Klocal(IU, JP) += tau * mt3d.N(ii) * mt3d.dNdx(jj, 0) / main.dt * mt3d.vol;
  Klocal(IV, JP) += tau * mt3d.N(ii) * mt3d.dNdx(jj, 1) / main.dt * mt3d.vol;
  Klocal(IW, JP) += tau * mt3d.N(ii) * mt3d.dNdx(jj, 2) / main.dt * mt3d.vol;

  // PSPG advection term
  for(int d = 0; d < 3; d++) {
    Klocal(IU, JP) += 5e-1 * tau * mt3d.dNdx(ii, d) * advk1[d] * mt3d.dNdx(jj, 0) * mt3d.vol;
    Klocal(IV, JP) += 5e-1 * tau * mt3d.dNdx(ii, d) * advk1[d] * mt3d.dNdx(jj, 1) * mt3d.vol;
    Klocal(IW, JP) += 5e-1 * tau * mt3d.dNdx(ii, d) * advk1[d] * mt3d.dNdx(jj, 2) * mt3d.vol;
  }

  // PSPG pressure term
  Klocal(IP, JP) += tau * mt3d.K(ii, jj) * mt3d.vol;
}

/**********************************************************************
 * @brief Compute element sfiffness RHS matrix on gauss integral point.
 */
void Adjoint::adjointGaussIntegralRHS(DirectProblem &main, VectorXd &Flocal, const double f, const int ii)
{
  int n1, n2, n3;
  mt3d.vol = mt3d.detJ * mt3d.weight;

  // Mass term
  Flocal(IU) += mt3d.N(ii) * wk1[0] / main.dt * mt3d.vol;
  Flocal(IV) += mt3d.N(ii) * wk1[1] / main.dt * mt3d.vol;
  Flocal(IW) += mt3d.N(ii) * wk1[2] / main.dt * mt3d.vol;

  /*
  // Diffusion term
  for(int d=0; d<3; d++){
      if(d == 0){n1 = 2e0; n2 = 1e0; n3 = 1e0;}
      if(d == 1){n1 = 1e0; n2 = 2e0; n3 = 1e0;}
      if(d == 2){n1 = 1e0; n2 = 1e0; n3 = 2e0;}
      Flocal(IU) -= 5e-1 * n1 * mt3d.dNdx(ii, d) * dwk1dx[0][d] / main.Re * mt3d.vol;
      Flocal(IV) -= 5e-1 * n2 * mt3d.dNdx(ii, d) * dwk1dx[1][d] / main.Re * mt3d.vol;
      Flocal(IW) -= 5e-1 * n3 * mt3d.dNdx(ii, d) * dwk1dx[2][d] / main.Re * mt3d.vol;
  }
  Flocal(IU) -= 5e-1 * mt3d.dNdx(ii, 1) * dwk1dx[1][0] / main.Re * mt3d.vol;
  Flocal(IU) -= 5e-1 * mt3d.dNdx(ii, 2) * dwk1dx[2][0] / main.Re * mt3d.vol;
  Flocal(IV) -= 5e-1 * mt3d.dNdx(ii, 0) * dwk1dx[0][1] / main.Re * mt3d.vol;
  Flocal(IV) -= 5e-1 * mt3d.dNdx(ii, 2) * dwk1dx[2][1] / main.Re * mt3d.vol;
  Flocal(IW) -= 5e-1 * mt3d.dNdx(ii, 0) * dwk1dx[0][2] / main.Re * mt3d.vol;
  Flocal(IW) -= 5e-1 * mt3d.dNdx(ii, 1) * dwk1dx[1][2] / main.Re * mt3d.vol;
  */

  // Diffusion term
  for(int d = 0; d < 3; d++) {
    Flocal(IU) -= 5e-1 * mt3d.dNdx(ii, d) * dwk1dx[0][d] / main.Re * mt3d.vol;
    Flocal(IV) -= 5e-1 * mt3d.dNdx(ii, d) * dwk1dx[1][d] / main.Re * mt3d.vol;
    Flocal(IW) -= 5e-1 * mt3d.dNdx(ii, d) * dwk1dx[2][d] / main.Re * mt3d.vol;
  }

  // Advection term
  Flocal(IU) -= 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[0][0] * wk1[0] * mt3d.vol;
  Flocal(IU) -= 0.5 * 1.5 * mt3d.N(ii) * dvkdx[0][0] * wk1[0] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IU) -= 0.5 * mt3d.dNdx(ii, d) * advk2[d] * wk1[0] * mt3d.vol;
  }
  Flocal(IU) -= 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[1][0] * wk1[1] * mt3d.vol;
  Flocal(IU) -= 0.5 * 1.5 * mt3d.N(ii) * dvkdx[1][0] * wk1[1] * mt3d.vol;
  Flocal(IU) -= 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[2][0] * wk1[2] * mt3d.vol;
  Flocal(IU) -= 0.5 * 1.5 * mt3d.N(ii) * dvkdx[2][0] * wk1[2] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IU) -= 0.5 * (-0.5) * mt3d.N(ii) * dvk2dx[d][0] * wk2[d] * mt3d.vol;
    Flocal(IU) -= 0.5 * (-0.5) * mt3d.N(ii) * dvk1dx[d][0] * wk2[d] * mt3d.vol;
  }

  Flocal(IV) -= 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[1][1] * wk1[1] * mt3d.vol;
  Flocal(IV) -= 0.5 * 1.5 * mt3d.N(ii) * dvkdx[1][1] * wk1[1] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IV) -= 0.5 * mt3d.dNdx(ii, d) * advk2[d] * wk1[1] * mt3d.vol;
  }
  Flocal(IV) -= 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[0][1] * wk1[0] * mt3d.vol;
  Flocal(IV) -= 0.5 * 1.5 * mt3d.N(ii) * dvkdx[0][1] * wk1[0] * mt3d.vol;
  Flocal(IV) -= 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[2][1] * wk1[2] * mt3d.vol;
  Flocal(IV) -= 0.5 * 1.5 * mt3d.N(ii) * dvkdx[2][1] * wk1[2] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IV) -= 0.5 * (-0.5) * mt3d.N(ii) * dvk2dx[d][1] * wk2[d] * mt3d.vol;
    Flocal(IV) -= 0.5 * (-0.5) * mt3d.N(ii) * dvk1dx[d][1] * wk2[d] * mt3d.vol;
  }

  Flocal(IW) -= 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[2][2] * wk1[2] * mt3d.vol;
  Flocal(IW) -= 0.5 * 1.5 * mt3d.N(ii) * dvkdx[2][2] * wk1[2] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IW) -= 0.5 * mt3d.dNdx(ii, d) * advk2[d] * wk1[2] * mt3d.vol;
  }
  Flocal(IW) -= 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[0][2] * wk1[0] * mt3d.vol;
  Flocal(IW) -= 0.5 * 1.5 * mt3d.N(ii) * dvkdx[0][2] * wk1[0] * mt3d.vol;
  Flocal(IW) -= 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[1][2] * wk1[1] * mt3d.vol;
  Flocal(IW) -= 0.5 * 1.5 * mt3d.N(ii) * dvkdx[1][2] * wk1[1] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IW) -= 0.5 * (-0.5) * mt3d.N(ii) * dvk2dx[d][2] * wk2[d] * mt3d.vol;
    Flocal(IW) -= 0.5 * (-0.5) * mt3d.N(ii) * dvk1dx[d][2] * wk2[d] * mt3d.vol;
  }

  // Darcy term
  Flocal(IU) -= 5e-1 * f * mt3d.N(ii) * wk1[0] * mt3d.vol;
  Flocal(IV) -= 5e-1 * f * mt3d.N(ii) * wk1[1] * mt3d.vol;
  Flocal(IW) -= 5e-1 * f * mt3d.N(ii) * wk1[2] * mt3d.vol;

  // SUPG mass term
  for(int d = 0; d < main.dim; d++) {
    Flocal(IU) -= tau * 1.5 * mt3d.N(ii) * (vk1[d] - vk[d]) / main.dt * dwk1dx[d][0] * mt3d.vol;
    Flocal(IU) += tau * 0.5 * mt3d.N(ii) * dwk2dx[d][0] * (vk2[d] - vk1[d]) / main.dt * mt3d.vol;
  }
  Flocal(IU) += tau * mt3d.N(ii) * advk2[0] * dwk1dx[0][0] / main.dt * mt3d.vol;

  for(int d = 0; d < main.dim; d++) {
    Flocal(IV) -= tau * 1.5 * mt3d.N(ii) * (vk1[d] - vk[d]) / main.dt * dwk1dx[d][1] * mt3d.vol;
    Flocal(IV) += tau * 0.5 * mt3d.N(ii) * dwk2dx[d][1] * (vk2[d] - vk1[d]) / main.dt * mt3d.vol;
  }
  Flocal(IV) += tau * mt3d.N(ii) * advk2[1] * dwk1dx[1][1] / main.dt * mt3d.vol;

  for(int d = 0; d < main.dim; d++) {
    Flocal(IW) -= tau * 1.5 * mt3d.N(ii) * (vk1[d] - vk[d]) / main.dt * dwk1dx[d][1] * mt3d.vol;
    Flocal(IW) += tau * 0.5 * mt3d.N(ii) * dwk2dx[d][2] * (vk2[d] - vk1[d]) / main.dt * mt3d.vol;
  }
  Flocal(IW) += tau * mt3d.N(ii) * advk2[2] * dwk1dx[2][2] / main.dt * mt3d.vol;

  // SUPG advection term
  std::vector<double> frontAdv2, frontAdv3;
  VecTool::resize(frontAdv2, main.dim);
  VecTool::resize(frontAdv3, main.dim);

  for(int d1 = 0; d1 < main.dim; d1++) {
    for(int d2 = 0; d2 < main.dim; d2++) {
      frontAdv2[d1] += advk2[d2] * dwk1dx[d1][d2];
      frontAdv3[d1] += advk3[d2] * dwk2dx[d1][d2];
    }
  }

  Flocal(IU) -= tau * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[0][0] * mt3d.vol;
  Flocal(IU) -= tau * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[0][0] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IU) -= tau * frontAdv2[0] * 0.5 * mt3d.dNdx(ii, d) * advk2[d] * mt3d.vol;
  }
  Flocal(IU) -= tau * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[1][0] * mt3d.vol;
  Flocal(IU) -= tau * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[1][0] * mt3d.vol;
  Flocal(IU) -= tau * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[2][0] * mt3d.vol;
  Flocal(IU) -= tau * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[2][0] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IU) -= tau * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(ii) * dvk2dx[d][0] * mt3d.vol;
    Flocal(IU) -= tau * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(ii) * dvk1dx[d][0] * mt3d.vol;
  }

  Flocal(IV) -= tau * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[1][1] * mt3d.vol;
  Flocal(IV) -= tau * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[1][1] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IV) -= tau * frontAdv2[1] * 0.5 * mt3d.dNdx(ii, d) * advk2[d] * mt3d.vol;
  }
  Flocal(IV) -= tau * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[0][1] * mt3d.vol;
  Flocal(IV) -= tau * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[0][1] * mt3d.vol;
  Flocal(IV) -= tau * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[2][1] * mt3d.vol;
  Flocal(IV) -= tau * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[2][1] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IV) -= tau * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(ii) * dvk2dx[d][1] * mt3d.vol;
    Flocal(IV) -= tau * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(ii) * dvk1dx[d][1] * mt3d.vol;
  }

  Flocal(IW) -= tau * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[2][2] * mt3d.vol;
  Flocal(IW) -= tau * frontAdv2[2] * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[2][2] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IW) -= tau * frontAdv2[2] * 0.5 * mt3d.dNdx(ii, d) * advk2[d] * mt3d.vol;
  }
  Flocal(IW) -= tau * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[0][2] * mt3d.vol;
  Flocal(IW) -= tau * frontAdv2[0] * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[0][2] * mt3d.vol;
  Flocal(IW) -= tau * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[1][2] * mt3d.vol;
  Flocal(IW) -= tau * frontAdv2[1] * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[1][2] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IW) -= tau * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(ii) * dvk2dx[d][2] * mt3d.vol;
    Flocal(IW) -= tau * frontAdv3[d] * 0.5 * (-0.5) * mt3d.N(ii) * dvk1dx[d][2] * mt3d.vol;
  }

  std::vector<double> backAdv2L, backAdv3L;
  VecTool::resize(backAdv2L, main.dim);
  VecTool::resize(backAdv3L, main.dim);

  for(int d1 = 0; d1 < main.dim; d1++) {
    for(int d2 = 0; d2 < main.dim; d2++) {
      backAdv2L[d1] += advk2[d2] * dvk1dx[d1][d2];
      backAdv3L[d1] += advk3[d2] * dvk2dx[d1][d2];
    }
  }

  std::vector<double> backAdv2R, backAdv3R;
  VecTool::resize(backAdv2R, main.dim);
  VecTool::resize(backAdv3R, main.dim);

  for(int d1 = 0; d1 < main.dim; d1++) {
    for(int d2 = 0; d2 < main.dim; d2++) {
      backAdv2R[d1] += advk2[d2] * dvkdx[d1][d2];
      backAdv3R[d1] += advk3[d2] * dvk1dx[d1][d2];
    }
  }

  for(int d = 0; d < main.dim; d++) {
    Flocal(IU) += tau * 0.5 * mt3d.N(ii) * dwk2dx[d][0] * (backAdv3L[0] + backAdv3R[0]) * mt3d.vol;
    Flocal(IU) -= tau * 1.5 * mt3d.N(ii) * dwk1dx[d][0] * (backAdv2L[0] + backAdv2R[0]) * mt3d.vol;
  }

  for(int d = 0; d < main.dim; d++) {
    Flocal(IV) += tau * 0.5 * mt3d.N(ii) * dwk2dx[d][1] * (backAdv3L[1] + backAdv3R[1]) * mt3d.vol;
    Flocal(IV) -= tau * 1.5 * mt3d.N(ii) * dwk1dx[d][1] * (backAdv2L[1] + backAdv2R[1]) * mt3d.vol;
  }

  for(int d = 0; d < main.dim; d++) {
    Flocal(IW) += tau * 0.5 * mt3d.N(ii) * dwk2dx[d][2] * (backAdv3L[2] + backAdv3R[2]) * mt3d.vol;
    Flocal(IW) -= tau * 1.5 * mt3d.N(ii) * dwk1dx[d][2] * (backAdv2L[2] + backAdv2R[2]) * mt3d.vol;
  }

  // SUPG pressure term
  for(int d = 0; d < 3; d++) {
    Flocal(IU) -= tau * 1.5 * mt3d.N(ii) * dpk1dx[d] * dwk1dx[d][0] * mt3d.vol;
    Flocal(IU) += tau * 0.5 * mt3d.N(ii) * dpk2dx[d] * dwk2dx[d][0] * mt3d.vol;
  }
  for(int d = 0; d < 3; d++) {
    Flocal(IV) -= tau * 1.5 * mt3d.N(ii) * dpk1dx[d] * dwk1dx[d][1] * mt3d.vol;
    Flocal(IV) += tau * 0.5 * mt3d.N(ii) * dpk2dx[d] * dwk2dx[d][1] * mt3d.vol;
  }
  for(int d = 0; d < 3; d++) {
    Flocal(IW) -= tau * 1.5 * mt3d.N(ii) * dpk1dx[d] * dwk1dx[d][2] * mt3d.vol;
    Flocal(IW) += tau * 0.5 * mt3d.N(ii) * dpk2dx[d] * dwk2dx[d][2] * mt3d.vol;
  }

  // PSPG mass term
  Flocal(IU) += tau * mt3d.N(ii) * dqk1dx[0] / main.dt * mt3d.vol;
  Flocal(IV) += tau * mt3d.N(ii) * dqk1dx[1] / main.dt * mt3d.vol;
  Flocal(IW) += tau * mt3d.N(ii) * dqk1dx[2] / main.dt * mt3d.vol;

  // PSPG advection term
  Flocal(IU) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[0][0] * dqk1dx[0] * mt3d.vol;
  Flocal(IU) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[0][0] * dqk1dx[0] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IU) -= tau * 0.5 * mt3d.dNdx(ii, d) * advk2[d] * dqk1dx[0] * mt3d.vol;
  }
  Flocal(IU) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[1][0] * dqk1dx[1] * mt3d.vol;
  Flocal(IU) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[1][0] * dqk1dx[1] * mt3d.vol;
  Flocal(IU) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[2][0] * dqk1dx[2] * mt3d.vol;
  Flocal(IU) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[2][0] * dqk1dx[2] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IU) -= tau * 0.5 * (-0.5) * mt3d.N(ii) * dvk2dx[d][0] * dqk2dx[d] * mt3d.vol;
    Flocal(IU) -= tau * 0.5 * (-0.5) * mt3d.N(ii) * dvk1dx[d][0] * dqk2dx[d] * mt3d.vol;
  }

  Flocal(IV) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[1][1] * dqk1dx[1] * mt3d.vol;
  Flocal(IV) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[1][1] * dqk1dx[1] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IV) -= tau * 0.5 * mt3d.dNdx(ii, d) * advk2[d] * dqk1dx[1] * mt3d.vol;
  }
  Flocal(IV) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[0][1] * dqk1dx[0] * mt3d.vol;
  Flocal(IV) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[0][1] * dqk1dx[0] * mt3d.vol;
  Flocal(IV) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[2][1] * dqk1dx[2] * mt3d.vol;
  Flocal(IV) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[2][1] * dqk1dx[2] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IV) -= tau * 0.5 * (-0.5) * mt3d.N(ii) * dvk2dx[d][1] * dqk2dx[d] * mt3d.vol;
    Flocal(IV) -= tau * 0.5 * (-0.5) * mt3d.N(ii) * dvk1dx[d][1] * dqk2dx[d] * mt3d.vol;
  }

  Flocal(IW) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[2][2] * dqk1dx[2] * mt3d.vol;
  Flocal(IW) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[2][2] * dqk1dx[2] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IW) -= tau * 0.5 * mt3d.dNdx(ii, d) * advk2[d] * dqk1dx[2] * mt3d.vol;
  }
  Flocal(IW) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[0][2] * dqk1dx[0] * mt3d.vol;
  Flocal(IW) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[0][2] * dqk1dx[0] * mt3d.vol;
  Flocal(IW) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvk1dx[1][2] * dqk1dx[1] * mt3d.vol;
  Flocal(IW) -= tau * 0.5 * 1.5 * mt3d.N(ii) * dvkdx[1][2] * dqk1dx[1] * mt3d.vol;
  for(int d = 0; d < main.dim; d++) {
    Flocal(IW) -= tau * 0.5 * (-0.5) * mt3d.N(ii) * dvk2dx[d][2] * dqk2dx[d] * mt3d.vol;
    Flocal(IW) -= tau * 0.5 * (-0.5) * mt3d.N(ii) * dvk1dx[d][2] * dqk2dx[d] * mt3d.vol;
  }
}

/***********************************************
 * @brief Set values needed for matrix assembly
 *        on gauss integral points.
 */
void Adjoint::setValue(DirectProblem &main, const int ic, const int t)
{  
  // main var - v
  for(int d = 0; d < main.dim; d++) {
    vk[d] = 0e0;
    vk1[d] = 0e0;
    vk2[d] = 0e0;
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int n = grid.cell(ic).node[p];
      if(t == timeMax - 1) {
        vk[d] += mt3d.N(p) * main.vt(t, n, d);
        vk1[d] = 0e0;
        vk2[d] = 0e0;
      } else if(t == timeMax - 2) {
        vk[d] += mt3d.N(p) * main.vt(t, n, d);
        vk1[d] += mt3d.N(p) * main.vt(t + 1, n, d);
        vk2[d] = 0e0;
      } else {
        vk[d] += mt3d.N(p) * main.vt(t, n, d);
        vk1[d] += mt3d.N(p) * main.vt(t + 1, n, d);
        vk2[d] += mt3d.N(p) * main.vt(t + 2, n, d);
      }
    }
  }

  // main var - dvdx
  for(int d = 0; d < main.dim; d++) {
    for(int e = 0; e < main.dim; e++) {
      dvkdx[d][e] = 0e0;
      dvk1dx[d][e] = 0e0;
      dvk2dx[d][e] = 0e0;
      for(int p = 0; p < grid.cell.nNodesInCell; p++) {
        int n = grid.cell(ic).node[p];
        if(t == timeMax - 1) {
          dvkdx[d][e] += mt3d.dNdx(p, e) * main.vt(t, n, d);
          dvk1dx[d][e] = 0e0;
          dvk2dx[d][e] = 0e0;
        } else if(t == timeMax - 2) {
          dvkdx[d][e] += mt3d.dNdx(p, e) * main.vt(t, n, d);
          dvk1dx[d][e] += mt3d.dNdx(p, e) * main.vt(t + 1, n, d);
          dvk2dx[d][e] = 0e0;
        } else {
          dvkdx[d][e] += mt3d.dNdx(p, e) * main.vt(t, n, d);
          dvk1dx[d][e] += mt3d.dNdx(p, e) * main.vt(t + 1, n, d);
          dvk2dx[d][e] += mt3d.dNdx(p, e) * main.vt(t + 2, n, d);
        }
      }
    }
  }

  // main var - dpdx
  for(int d = 0; d < main.dim; d++) {
    dpkdx[d] = 0e0;
    dpk1dx[d] = 0e0;
    dpk2dx[d] = 0e0;
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int n = grid.cell(ic).node[p];
      if(t == timeMax - 1) {
        dpkdx[d] += mt3d.dNdx(p, d) * main.pt(t, n);
        dpk1dx[d] = 0e0;
        dpk2dx[d] = 0e0;
      } else if(t == timeMax - 2) {
        dpkdx[d] += mt3d.dNdx(p, d) * main.pt(t, n);
        dpk1dx[d] += mt3d.dNdx(p, d) * main.pt(t + 1, n);
        dpk2dx[d] = 0e0;
      } else {
        dpkdx[d] += mt3d.dNdx(p, d) * main.pt(t, n);
        dpk1dx[d] += mt3d.dNdx(p, d) * main.pt(t + 1, n);
        dpk2dx[d] += mt3d.dNdx(p, d) * main.pt(t + 2, n);
      }
    }
  }

  // main var - adv
  for(int d = 0; d < main.dim; d++) {
    advk1[d] = 0e0;
    advk2[d] = 0e0;
    advk3[d] = 0e0;
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int n = grid.cell(ic).node[p];
      if(t == 0) {
        advk1[d] += mt3d.N(p) * main.v0(n, d);
        advk2[d] += mt3d.N(p) * (1.5 * main.vt(t, n, d) - 0.5 * main.v0(n, d));
        advk3[d] += mt3d.N(p) * (1.5 * main.vt(t + 1, n, d) - 0.5 * main.vt(t, n, d));
      } else if(t == 1) {
        advk1[d] += mt3d.N(p) * (1.5 * main.vt(t - 1, n, d) - 0.5 * main.v0(n, d));
        advk2[d] += mt3d.N(p) * (1.5 * main.vt(t, n, d) - 0.5 * main.vt(t - 1, n, d));
        advk3[d] += mt3d.N(p) * (1.5 * main.vt(t + 1, n, d) - 0.5 *main.vt(t, n, d));
      } else if(t == timeMax - 1) {
        advk1[d] += mt3d.N(p) * (1.5 * main.vt(t - 1, n, d) - 0.5 * main.vt(t - 2, n, d));
        advk2[d] += mt3d.N(p) * (1.5 * main.vt(t, d, d) - 0.5 * main.vt(t - 1, n, d));
        advk3[d] += mt3d.N(p) * main.vt(t, n, d);
      } else {
        advk1[d] += mt3d.N(p) * (1.5 * main.vt(t - 1, n, d) - 0.5 * main.vt(t - 2, n, d));
        advk2[d] += mt3d.N(p) * (1.5 * main.vt(t, n, d) - 0.5 * main.vt(t - 1, n, d));
        advk3[d] += mt3d.N(p) * (1.5 * main.vt(t + 1, n, d) - 0.5 * main.vt(t, n, d));
      }
    }
  }

  // lagrange multiplier - w
  for(int d = 0; d < main.dim; d++) {
    wk1[d] = 0e0;
    wk2[d] = 0e0;
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int n = grid.cell(ic).node[p];
      wk1[d] += mt3d.N(p) * w(n, d);
      wk2[d] += mt3d.N(p) * wPrev(n, d);
    }
  }

  // lagrange multiplier - dwdx
  for(int d = 0; d < main.dim; d++) {
    for(int e = 0; e < main.dim; e++) {
      dwk1dx[d][e] = 0e0;
      dwk2dx[d][e] = 0e0;
      for(int p = 0; p < grid.cell.nNodesInCell; p++) {
        int n = grid.cell(ic).node[p];
        dwk1dx[d][e] += mt3d.dNdx(p, e) * w(n, d);
        dwk2dx[d][e] += mt3d.dNdx(p, e) * wPrev(n, d);
      }
    }
  }

  // lagrange multiplier - dqdx
  for(int d = 0; d < main.dim; d++) {
    dqk1dx[d] = 0e0;
    dqk2dx[d] = 0e0;
    for(int p = 0; p < grid.cell.nNodesInCell; p++) {
      int n = grid.cell(ic).node[p];
      dqk1dx[d] += mt3d.dNdx(p, d) * q(n);
      dqk2dx[d] += mt3d.dNdx(p, d) * qPrev(n);
    }
  }
  
}

/***********************************
 * @brief Compute boundary integral.
 */
void Adjoint::boundaryIntegral(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal, ControlBoundary &cb,
                               const int ic, const int ib)
{
  Klocal.setZero();
  Flocal.setZero();

  int nc = 4;

  mt2d.N.allocate(nc);
  mt2d.xCurrent.allocate(nc, 2);
  mt2d.dNdr.allocate(nc, 2);
  mt2d.dNdx.allocate(nc, 2);

  Gauss g2(2);

  for(int p = 0; p < nc; p++) {
    int n = cb.CBNodeMapInCell[ib][p];
    for(int d = 0; d < 2; d++) {
      mt2d.xCurrent(p, d) = 0e0;
      mt2d.xCurrent(p, d) = main.grid.node.x[n][planeDir[d]];
    }
  }

  for(int i1 = 0; i1 < 2; i1++) {
    for(int i2 = 0; i2 < 2; i2++) {
      ShapeFunction2D::C2D4_N(mt2d.N, g2.point[i1], g2.point[i2]);
      ShapeFunction2D::C2D4_dNdr(mt2d.dNdr, g2.point[i1], g2.point[i2]);
      MathTools2D::comp_dxdr(mt2d.dxdr, mt2d.dNdr, mt2d.xCurrent, nc);
      mt2d.detJ = MathTools2D::compDeterminant(mt2d.dxdr);
      mt2d.weight = g2.weight[i1] * g2.weight[i2];
      for(int ii = 0; ii < nc; ii++) {
        updateRowIndexPlane(grid, ii, ic);
        for(int jj = 0; jj < nc; jj++) {
          updateColumnIndexPlane(grid, jj, ic);
          boundaryInGaussIntegral(Klocal, ii, jj);
        }
      }
    }
  }
}

/*******************************************************************
 * @brief Compute LHS element sfiffness matrix for boundary integral
 *        on gauss integral point.
 */
void Adjoint::boundaryInGaussIntegral(MatrixXd &Klocal, const int ii, const int jj)
{
  mt2d.vol = mt2d.detJ * mt2d.weight;

  Klocal(IU, JLU) += mt2d.N(ii) * mt2d.N(jj) * mt2d.vol;
  Klocal(IV, JLV) += mt2d.N(ii) * mt2d.N(jj) * mt2d.vol;
  Klocal(IW, JLW) += mt2d.N(ii) * mt2d.N(jj) * mt2d.vol;
  Klocal(ILU, JU) += mt2d.N(ii) * mt2d.N(jj) * mt2d.vol;
  Klocal(ILV, JV) += mt2d.N(ii) * mt2d.N(jj) * mt2d.vol;
  Klocal(ILW, JW) += mt2d.N(ii) * mt2d.N(jj) * mt2d.vol;
}

/* Advection term for Optimize Then Discretize
 *
 * for(int d=0; d<main.dim; d++){
 *     Klocal(IU, JU) += 5e-1 * dNdx(ii, d) * vel[d] * N(jj) * mt3d.vol;
 *     Klocal(IV, JV) += 5e-1 * dNdx(ii, d) * vel[d] * N(jj) * mt3d.vol;
 *     Klocal(IW, JW) += 5e-1 * dNdx(ii, d) * vel[d] * N(jj) * mt3d.vol;
 *  }
 *  Klocal(IU, JV) += 5e-1 * N(ii) * dvdx[1][0] * N(jj) * mt3d.vol;
 *  Klocal(IU, JW) += 5e-1 * N(ii) * dvdx[2][0] * N(jj) * mt3d.vol;
 *  Klocal(IV, JU) += 5e-1 * N(ii) * dvdx[0][1] * N(jj) * mt3d.vol;
 *  Klocal(IV, JW) += 5e-1 * N(ii) * dvdx[2][1] * N(jj) * mt3d.vol;
 *  Klocal(IW, JU) += 5e-1 * N(ii) * dvdx[0][2] * N(jj) * mt3d.vol;
 *  Klocal(IW, JV) += 5e-1 * N(ii) * dvdx[1][2] * N(jj) * mt3d.vol;
 */