/**
 * @file AssemblyUSNS.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "DirectProblem.h"

/**
 * @brief Assemble system.
 */
void DirectProblem::matrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t)
{
  MathTools3D tools(grid.cell.nNodesInCell);

  Klocal.setZero();
  Flocal.setZero();

  for(int p = 0; p < grid.cell.nNodesInCell; p++) {
    int n = grid.cell(ic).node[p];
    for(int d = 0; d < 3; d++) {
      tools.xCurrent(p, d) = 0e0;
      tools.xCurrent(p, d) = grid.node.x[n][d];
    }
  }

  double he = comp_he(tools.xCurrent);
  double f = comp_f(grid.cell(ic).phi);

  Gauss g2(2);

  for(int i1 = 0; i1 < 2; i1++) {
    for(int i2 = 0; i2 < 2; i2++) {
      for(int i3 = 0; i3 < 2; i3++) {
        settingInGauss(tools, g2, he, i1, i2, i3, ic, t);
        setLocalValuesInGauss(Klocal, Flocal, tools, f, ic);
      }
    }
  }
}

/**
 * @brief Assemble local matrix and vector.
 */
void DirectProblem::setLocalValuesInGauss(MatrixXd &Klocal, VectorXd &Flocal, MathTools3D &tools, double f,
                                                 const int ic)
{
  for(int ii = 0; ii < grid.cell.nNodesInCell; ii++) {
    updateRowIndex(grid, ii, ic);
    for(int jj = 0; jj < grid.cell.nNodesInCell; jj++) {
      updateColumnIndex(grid, jj, ic);
      usnsGaussIntegralLHS(Klocal, tools, f, ii, jj);
    }
    usnsGaussIntegralRHS(Flocal, tools, f, ii);
  }
}

/**
 * @brief Set Values in gauss integral.
 */
void DirectProblem::settingInGauss(MathTools3D &tools, Gauss &g2, const double he, const int i1, const int i2,
                                             const int i3, const int ic, const int t)
{
  ShapeFunction3D::C3D8_N(tools.N, g2.point[i1], g2.point[i2], g2.point[i3]);
  ShapeFunction3D::C3D8_dNdr(tools.dNdr, g2.point[i1], g2.point[i2], g2.point[i3]);

  MathTools3D::comp_dxdr(tools.dxdr, tools.dNdr, tools.xCurrent, grid.cell.nNodesInCell);
  MathTools3D::comp_dNdx(tools.dNdx, tools.dNdr, tools.dxdr, grid.cell.nNodesInCell);

  tools.detJ = MathTools3D::compDeterminant(tools.dxdr);
  tools.weight = g2.weight[i1] * g2.weight[i2] * g2.weight[i3];

  setVelocityValue(tools, ic, t);
  tau = comp_tau(adv_gp, he);
  //tau = comp_tau2(tools.dNdx, adv_gp, tools.nNodesInCell);
}

/**
 * @brief Compute element sfiffness LHS matrix on gauss integral point.
 */
void DirectProblem::usnsGaussIntegralLHS(MatrixXd &Klocal, MathTools3D &tools, const double f, const int ii,
                                         const int jj)
{
  tools.vol = tools.detJ * tools.weight;

  tools.K(ii, jj) = 0e0;
  for(int d = 0; d < 3; d++) {
    tools.K(ii, jj) += tools.dNdx(ii, d) * tools.dNdx(jj, d);
  }

  // Mass term
  Klocal(IU, JU) += rho * tools.N(ii) * tools.N(jj) / dt * tools.vol;
  Klocal(IV, JV) += rho * tools.N(ii) * tools.N(jj) / dt * tools.vol;
  Klocal(IW, JW) += rho * tools.N(ii) * tools.N(jj) / dt * tools.vol;

  // Diffusion term
  for(int d = 0; d < 3; d++) {
    Klocal(IU, JU) += 5e-1 * mu * tools.dNdx(ii, d) * tools.dNdx(jj, d) * tools.vol;
    Klocal(IV, JV) += 5e-1 * mu * tools.dNdx(ii, d) * tools.dNdx(jj, d) * tools.vol;
    Klocal(IW, JW) += 5e-1 * mu * tools.dNdx(ii, d) * tools.dNdx(jj, d) * tools.vol;
  }

  // Advection term
  for(int d = 0; d < 3; d++) {
    Klocal(IU, JU) += 5e-1 * rho * tools.N(ii) * adv_gp[d] * tools.dNdx(jj, d) * tools.vol;
    Klocal(IV, JV) += 5e-1 * rho * tools.N(ii) * adv_gp[d] * tools.dNdx(jj, d) * tools.vol;
    Klocal(IW, JW) += 5e-1 * rho * tools.N(ii) * adv_gp[d] * tools.dNdx(jj, d) * tools.vol;
  }

  // Pressure term
  Klocal(IU, JP) -= tools.N(jj) * tools.dNdx(ii, 0) * tools.vol;
  Klocal(IV, JP) -= tools.N(jj) * tools.dNdx(ii, 1) * tools.vol;
  Klocal(IW, JP) -= tools.N(jj) * tools.dNdx(ii, 2) * tools.vol;

  // Continuity term
  Klocal(IP, JU) += tools.N(ii) * tools.dNdx(jj, 0) * tools.vol;
  Klocal(IP, JV) += tools.N(ii) * tools.dNdx(jj, 1) * tools.vol;
  Klocal(IP, JW) += tools.N(ii) * tools.dNdx(jj, 2) * tools.vol;

  // Darcy term
  Klocal(IU, JU) += 5e-1 * f / rho * tools.N(ii) * tools.N(jj) * tools.vol;
  Klocal(IV, JV) += 5e-1 * f / rho * tools.N(ii) * tools.N(jj) * tools.vol;
  Klocal(IW, JW) += 5e-1 * f / rho * tools.N(ii) * tools.N(jj) * tools.vol;

  // SUPG　mass term
  for(int d = 0; d < 3; d++) {
    Klocal(IU, JU) += tau * rho * tools.dNdx(ii, d) * adv_gp[d] * tools.N(jj) / dt * tools.vol;
    Klocal(IV, JV) += tau * rho * tools.dNdx(ii, d) * adv_gp[d] * tools.N(jj) / dt * tools.vol;
    Klocal(IW, JW) += tau * rho * tools.dNdx(ii, d) * adv_gp[d] * tools.N(jj) / dt * tools.vol;
  }

  // SUPG advection term
  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 3; d2++) {
      Klocal(IU, JU) += 5e-1 * tau * rho * adv_gp[d2] * tools.dNdx(ii, d2) * adv_gp[d1] * tools.dNdx(jj, d1) * tools.vol;
      Klocal(IV, JV) += 5e-1 * tau * rho * adv_gp[d2] * tools.dNdx(ii, d2) * adv_gp[d1] * tools.dNdx(jj, d1) * tools.vol;
      Klocal(IW, JW) += 5e-1 * tau * rho * adv_gp[d2] * tools.dNdx(ii, d2) * adv_gp[d1] * tools.dNdx(jj, d1) * tools.vol;
    }
  }

  // SUPG pressure term
  for(int d = 0; d < 3; d++) {
    Klocal(IU, JP) += tau * tools.dNdx(ii, d) * adv_gp[d] * tools.dNdx(jj, 0) * tools.vol;
    Klocal(IV, JP) += tau * tools.dNdx(ii, d) * adv_gp[d] * tools.dNdx(jj, 1) * tools.vol;
    Klocal(IW, JP) += tau * tools.dNdx(ii, d) * adv_gp[d] * tools.dNdx(jj, 2) * tools.vol;
  }

  // PSPG mass term
  Klocal(IP, JU) += tau * tools.dNdx(ii, 0) * tools.N(jj) / dt * tools.vol;
  Klocal(IP, JV) += tau * tools.dNdx(ii, 1) * tools.N(jj) / dt * tools.vol;
  Klocal(IP, JW) += tau * tools.dNdx(ii, 2) * tools.N(jj) / dt * tools.vol;

  // PSPG asvection term
  for(int d = 0; d < 3; d++) {
    Klocal(IP, JU) += 5e-1 * tau * tools.dNdx(ii, 0) * adv_gp[d] * tools.dNdx(jj, d) * tools.vol;
    Klocal(IP, JV) += 5e-1 * tau * tools.dNdx(ii, 1) * adv_gp[d] * tools.dNdx(jj, d) * tools.vol;
    Klocal(IP, JW) += 5e-1 * tau * tools.dNdx(ii, 2) * adv_gp[d] * tools.dNdx(jj, d) * tools.vol;
  }

  // PSPG pressure term
  Klocal(IP, JP) += tau / rho * tools.K(ii, jj) * tools.vol;
}

/**
 * @brief Compute element sfiffness RHS vector on gauss integral point.
 */
void DirectProblem::usnsGaussIntegralRHS(VectorXd &Flocal, MathTools3D &tools, const double f, const int ii)
{
  int n1, n2, n3;
  tools.vol = tools.detJ * tools.weight;

  // Mass term
  Flocal(IU) += rho * tools.N(ii) * v_gp[0] / dt * tools.vol;
  Flocal(IV) += rho * tools.N(ii) * v_gp[1] / dt * tools.vol;
  Flocal(IW) += rho * tools.N(ii) * v_gp[2] / dt * tools.vol;

  // Diffusion term
  for(int d = 0; d < 3; d++) {
    Flocal(IU) -= 5e-1 * mu * tools.dNdx(ii, d) * dvdx_gp[0][d] * tools.vol;
    Flocal(IV) -= 5e-1 * mu * tools.dNdx(ii, d) * dvdx_gp[1][d] * tools.vol;
    Flocal(IW) -= 5e-1 * mu * tools.dNdx(ii, d) * dvdx_gp[2][d] * tools.vol;
  }

  // Advection term
  for(int d = 0; d < 3; d++) {
    Flocal(IU) -= 5e-1 * rho * tools.N(ii) * adv_gp[d] * dvdx_gp[0][d] * tools.vol;
    Flocal(IV) -= 5e-1 * rho * tools.N(ii) * adv_gp[d] * dvdx_gp[1][d] * tools.vol;
    Flocal(IW) -= 5e-1 * rho * tools.N(ii) * adv_gp[d] * dvdx_gp[2][d] * tools.vol;
  }

  // Darcy term
  Flocal(IU) -= 5e-1 * f / rho * tools.N(ii) * v_gp[0] * tools.vol;
  Flocal(IV) -= 5e-1 * f / rho * tools.N(ii) * v_gp[1] * tools.vol;
  Flocal(IW) -= 5e-1 * f / rho * tools.N(ii) * v_gp[2] * tools.vol;

  // SUPG mass term
  for(int d = 0; d < 3; d++) {
    Flocal(IU) += tau * rho * tools.dNdx(ii, d) * adv_gp[d] * v_gp[0] / dt * tools.vol;
    Flocal(IV) += tau * rho * tools.dNdx(ii, d) * adv_gp[d] * v_gp[1] / dt * tools.vol;
    Flocal(IW) += tau * rho * tools.dNdx(ii, d) * adv_gp[d] * v_gp[2] / dt * tools.vol;
  }

  // SUPG advection term
  for(int d1 = 0; d1 < 3; d1++) {
    for(int d2 = 0; d2 < 3; d2++) {
      Flocal(IU) -= 5e-1 * tau * rho * v_gp[d2] * tools.dNdx(ii, d2) * adv_gp[d1] * dvdx_gp[0][d1] * tools.vol;
      Flocal(IV) -= 5e-1 * tau * rho * v_gp[d2] * tools.dNdx(ii, d2) * adv_gp[d1] * dvdx_gp[1][d1] * tools.vol;
      Flocal(IW) -= 5e-1 * tau * rho * v_gp[d2] * tools.dNdx(ii, d2) * adv_gp[d1] * dvdx_gp[2][d1] * tools.vol;
    }
  }

  // PSPG mass term
  Flocal(IP) += tau * tools.dNdx(ii, 0) * v_gp[0] / dt * tools.vol;
  Flocal(IP) += tau * tools.dNdx(ii, 1) * v_gp[1] / dt * tools.vol;
  Flocal(IP) += tau * tools.dNdx(ii, 2) * v_gp[2] / dt * tools.vol;

  // PSPG advection term
  for(int d = 0; d < 3; d++) {
    Flocal(IP) -= 5e-1 * tau * tools.dNdx(ii, 0) * adv_gp[d] * dvdx_gp[0][d] * tools.vol;
    Flocal(IP) -= 5e-1 * tau * tools.dNdx(ii, 1) * adv_gp[d] * dvdx_gp[1][d] * tools.vol;
    Flocal(IP) -= 5e-1 * tau * tools.dNdx(ii, 2) * adv_gp[d] * dvdx_gp[2][d] * tools.vol;
  }
}

void DirectProblem::setVelocityValue(MathTools3D &tools, const int ic, const int t)
{
  compVelOnGaussPoint(grid, tools, v_gp, ic);
  compAdVelOnGaussPoint(grid, tools, adv_gp, ic);
  compVelDerivOnGaussPoint(grid, tools, dvdx_gp, ic);
}

/*
   // Diffusion term
  for(int d=0; d<3; d++){
      if(d == 0){n1 = 2; n2 = 1; n3 = 1;}
      if(d == 1){n1 = 1; n2 = 2; n3 = 1;}
      if(d == 2){n1 = 1; n2 = 1; n3 = 2;}
      Klocal(IU, JU) += 5e-1 * n1 * tools.dNdx(ii, d) * tools.dNdx(jj, d) / Re * tools.vol;
      Klocal(IV, JV) += 5e-1 * n2 * tools.dNdx(ii, d) * tools.dNdx(jj, d) / Re * tools.vol;
      Klocal(IW, JW) += 5e-1 * n3 * tools.dNdx(ii, d) * tools.dNdx(jj, d) / Re * tools.vol;
  }
  Klocal(IU, JV) += 5e-1 * tools.dNdx(ii, 1) * tools.dNdx(jj, 0) / Re * tools.vol;
  Klocal(IU, JW) += 5e-1 * tools.dNdx(ii, 2) * tools.dNdx(jj, 0) / Re * tools.vol;
  Klocal(IV, JU) += 5e-1 * tools.dNdx(ii, 0) * tools.dNdx(jj, 1) / Re * tools.vol;
  Klocal(IV, JW) += 5e-1 * tools.dNdx(ii, 2) * tools.dNdx(jj, 1) / Re * tools.vol;
  Klocal(IW, JU) += 5e-1 * tools.dNdx(ii, 0) * tools.dNdx(jj, 2) / Re * tools.vol;
  Klocal(IW, JV) += 5e-1 * tools.dNdx(ii, 1) * tools.dNdx(jj, 2) / Re * tools.vol;
*/
