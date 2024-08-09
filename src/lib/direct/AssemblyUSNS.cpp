/**
 * @file AssemblyUSNS.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "DirectProblem.h"

/**************************
 * @brief Assemble system.
 */
void DirectProblem::matrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, MathTools3D &tools, const int ic, const int t)
{
  for (int p = 0; p < grid.cell.nNodesInCell; p++)
  {
    for (int d = 0; d < dim; d++)
    {
      tools.xCurrent(p, d) = 0e0;
      tools.xCurrent(p, d) = grid.node.x[grid.cell(ic).node[p]][d];
    }
  }

  double he = fabs(tools.xCurrent(1, 0) - tools.xCurrent(0, 0));
  double f = resistance * alpha * (1e0 - grid.cell(ic).phi) / (alpha + grid.cell(ic).phi);

  Gauss g2(2);
  for (int i1 = 0; i1 < 2; i1++)
  {
    for (int i2 = 0; i2 < 2; i2++)
    {
      for (int i3 = 0; i3 < 2; i3++)
      {
        setValuesInGaussIntegral(tools, g2, he, i1, i2, i3, ic, t);
        for (int ii = 0; ii < grid.cell.nNodesInCell; ii++)
        {
          updateRowIndex(ii, ic);
          for (int jj = 0; jj < grid.cell.nNodesInCell; jj++)
          {
            updateColumnIndex(jj, ic);
            usnsGaussIntegralLHS(Klocal, tools, f, ii, jj);
          }
          usnsGaussIntegralRHS(Flocal, tools, f, ii);
        }
      }
    }
  }
}

/**************************************
 * @brief Set Values in gauss integral.
 */
void DirectProblem::setValuesInGaussIntegral(MathTools3D &tools, Gauss &g2, const double he,
                                             const int i1, const int i2, const int i3,
                                             const int ic, const int t)
{
  double dxdr[3][3];
  ShapeFunction3D::C3D8_N(tools.N, g2.point[i1], g2.point[i2], g2.point[i3]);
  ShapeFunction3D::C3D8_dNdr(tools.dNdr, g2.point[i1], g2.point[i2], g2.point[i3]);

  MathTools3D::comp_dxdr(dxdr, tools.dNdr, tools.xCurrent, grid.cell.nNodesInCell);
  MathTools3D::comp_dNdx(tools.dNdx, tools.dNdr, dxdr, grid.cell.nNodesInCell);

  tools.detJ = MathCommon::compDeterminant_3x3(dxdr);
  tools.weight = g2.weight[i1] * g2.weight[i2] * g2.weight[i3];

  setVelocityValue(tools, ic, t);
  tau = MathCommon::comp_tau(advgp, he, Re, dt);
}

/**********************************************************************
 * @brief Compute element sfiffness LHS matrix on gauss integral point.
 */
void DirectProblem::usnsGaussIntegralLHS(MatrixXd &Klocal, MathTools3D &tools, const double f, const int ii, const int jj)
{
  int n1, n2, n3;
  tools.vol = tools.detJ * tools.weight;

  tools.K(ii, jj) = 0e0;
  for (int d = 0; d < 3; d++)
  {
    tools.K(ii, jj) += tools.dNdx(ii, d) * tools.dNdx(jj, d);
  }

  // Mass term
  Klocal(IU, JU) += tools.N(ii) * tools.N(jj) / dt * tools.vol;
  Klocal(IV, JV) += tools.N(ii) * tools.N(jj) / dt * tools.vol;
  Klocal(IW, JW) += tools.N(ii) * tools.N(jj) / dt * tools.vol;

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

  // Diffusion term
  for (int d = 0; d < 3; d++)
  {
    Klocal(IU, JU) += 5e-1 * tools.dNdx(ii, d) * tools.dNdx(jj, d) / Re * tools.vol;
    Klocal(IV, JV) += 5e-1 * tools.dNdx(ii, d) * tools.dNdx(jj, d) / Re * tools.vol;
    Klocal(IW, JW) += 5e-1 * tools.dNdx(ii, d) * tools.dNdx(jj, d) / Re * tools.vol;
  }

  // Advection term
  for (int d = 0; d < 3; d++)
  {
    Klocal(IU, JU) += 5e-1 * tools.N(jj) * advgp[d] * tools.dNdx(jj, d) * tools.vol;
    Klocal(IV, JV) += 5e-1 * tools.N(jj) * advgp[d] * tools.dNdx(jj, d) * tools.vol;
    Klocal(IW, JW) += 5e-1 * tools.N(jj) * advgp[d] * tools.dNdx(jj, d) * tools.vol;
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
  Klocal(IU, JU) += 5e-1 * f * tools.N(ii) * tools.N(jj) * tools.vol;
  Klocal(IV, JV) += 5e-1 * f * tools.N(ii) * tools.N(jj) * tools.vol;
  Klocal(IW, JW) += 5e-1 * f * tools.N(ii) * tools.N(jj) * tools.vol;

  // SUPG　mass term
  for (int d = 0; d < 3; d++)
  {
    Klocal(IU, JU) += tau * tools.dNdx(ii, d) * advgp[d] * tools.N(ii) / dt * tools.vol;
    Klocal(IV, JV) += tau * tools.dNdx(ii, d) * advgp[d] * tools.N(ii) / dt * tools.vol;
    Klocal(IW, JW) += tau * tools.dNdx(ii, d) * advgp[d] * tools.N(ii) / dt * tools.vol;
  }

  // SUPG advection term
  for (int d1 = 0; d1 < 3; d1++)
  {
    for (int d2 = 0; d2 < 3; d2++)
    {
      Klocal(IU, JU) += 5e-1 * tau * advgp[d2] * tools.dNdx(ii, d2) * advgp[d1] * tools.dNdx(jj, d1) * tools.vol;
      Klocal(IV, JV) += 5e-1 * tau * advgp[d2] * tools.dNdx(ii, d2) * advgp[d1] * tools.dNdx(jj, d1) * tools.vol;
      Klocal(IW, JW) += 5e-1 * tau * advgp[d2] * tools.dNdx(ii, d2) * advgp[d1] * tools.dNdx(jj, d1) * tools.vol;
    }
  }

  // SUPG pressure term
  for (int d = 0; d < 3; d++)
  {
    Klocal(IU, JP) += tau * tools.dNdx(ii, d) * advgp[d] * tools.dNdx(jj, 0) * tools.vol;
    Klocal(IV, JP) += tau * tools.dNdx(ii, d) * advgp[d] * tools.dNdx(ii, 1) * tools.vol;
    Klocal(IW, JP) += tau * tools.dNdx(ii, d) * advgp[d] * tools.dNdx(ii, 2) * tools.vol;
  }

  // PSPG mass term
  Klocal(IP, JU) += tau * tools.dNdx(ii, 0) * tools.N(jj) / dt * tools.vol;
  Klocal(IP, JV) += tau * tools.dNdx(ii, 1) * tools.N(jj) / dt * tools.vol;
  Klocal(IP, JW) += tau * tools.dNdx(ii, 2) * tools.N(jj) / dt * tools.vol;

  // PSPG asvection term
  for (int d = 0; d < 3; d++)
  {
    Klocal(IP, JU) += 5e-1 * tau * tools.dNdx(ii, 0) * advgp[d] * tools.dNdx(jj, d) * tools.vol;
    Klocal(IP, JV) += 5e-1 * tau * tools.dNdx(ii, 1) * advgp[d] * tools.dNdx(jj, d) * tools.vol;
    Klocal(IP, JW) += 5e-1 * tau * tools.dNdx(ii, 2) * advgp[d] * tools.dNdx(jj, d) * tools.vol;
  }

  // PSPG pressure term
  Klocal(IP, JP) += tau * tools.K(ii, jj) * tools.vol;
}

/**********************************************************************
 * @brief Compute element sfiffness RHS vector on gauss integral point.
 */
void DirectProblem::usnsGaussIntegralRHS(VectorXd &Flocal, MathTools3D &tools, const double f, const int ii)
{
  int n1, n2, n3;
  tools.vol = tools.detJ * tools.weight;

  // Mass term
  Flocal(IU) += tools.N(ii) * vgp[0] / dt * tools.vol;
  Flocal(IV) += tools.N(ii) * vgp[1] / dt * tools.vol;
  Flocal(IW) += tools.N(ii) * vgp[2] / dt * tools.vol;

  /*
  // Diffusion term
  for(int d=0; d<3; d++){
      if(d == 0){n1 = 2; n2 = 1; n3 = 1;}
      if(d == 1){n1 = 1; n2 = 2; n3 = 1;}
      if(d == 2){n1 = 1; n2 = 1; n3 = 2;}
      Flocal(IU) -= 5e-1 * n1 * tools.dNdx(ii, d) * dvgpdx[0][d] / Re * tools.vol;
      Flocal(IV) -= 5e-1 * n2 * tools.dNdx(ii, d) * dvgpdx[1][d] / Re * tools.vol;
      Flocal(IW) -= 5e-1 * n3 * tools.dNdx(ii, d) * dvgpdx[2][d] / Re * tools.vol;
  }
  Flocal(IU) -= 5e-1 * tools.dNdx(ii, 1) * dvgpdx[1][0] / Re * tools.vol;
  Flocal(IU) -= 5e-1 * tools.dNdx(ii, 2) * dvgpdx[2][0] / Re * tools.vol;
  Flocal(IV) -= 5e-1 * tools.dNdx(ii, 0) * dvgpdx[0][1] / Re * tools.vol;
  Flocal(IV) -= 5e-1 * tools.dNdx(ii, 2) * dvgpdx[2][1] / Re * tools.vol;
  Flocal(IW) -= 5e-1 * tools.dNdx(ii, 0) * dvgpdx[0][2] / Re * tools.vol;
  Flocal(IW) -= 5e-1 * tools.dNdx(ii, 1) * dvgpdx[1][2] / Re * tools.vol;
  */

  // Diffusion term
  for (int d = 0; d < 3; d++)
  {
    Flocal(IU) -= 5e-1 * tools.dNdx(ii, d) * dvgpdx[0][d] / Re * tools.vol;
    Flocal(IV) -= 5e-1 * tools.dNdx(ii, d) * dvgpdx[1][d] / Re * tools.vol;
    Flocal(IW) -= 5e-1 * tools.dNdx(ii, d) * dvgpdx[2][d] / Re * tools.vol;
  }

  // Advection term
  for (int d = 0; d < 3; d++)
  {
    Flocal(IU) -= 5e-1 * tools.N(ii) * advgp[d] * dvgpdx[0][d] * tools.vol;
    Flocal(IV) -= 5e-1 * tools.N(ii) * advgp[d] * dvgpdx[1][d] * tools.vol;
    Flocal(IW) -= 5e-1 * tools.N(ii) * advgp[d] * dvgpdx[2][d] * tools.vol;
  }

  // Darcy term
  Flocal(IU) -= 5e-1 * f * tools.N(ii) * vgp[0] * tools.vol;
  Flocal(IV) -= 5e-1 * f * tools.N(ii) * vgp[1] * tools.vol;
  Flocal(IW) -= 5e-1 * f * tools.N(ii) * vgp[2] * tools.vol;

  // SUPG mass term
  for (int d = 0; d < 3; d++)
  {
    Flocal(IU) += tau * tools.dNdx(ii, d) * advgp[d] * vgp[0] / dt * tools.vol;
    Flocal(IV) += tau * tools.dNdx(ii, d) * advgp[d] * vgp[1] / dt * tools.vol;
    Flocal(IW) += tau * tools.dNdx(ii, d) * advgp[d] * vgp[2] / dt * tools.vol;
  }

  // SUPG advection term
  for (int d1 = 0; d1 < 3; d1++)
  {
    for (int d2 = 0; d2 < 3; d2++)
    {
      Flocal(IU) -= 5e-1 * tau * vgp[d2] * tools.dNdx(ii, d2) * advgp[d1] * dvgpdx[0][d1] * tools.vol;
      Flocal(IV) -= 5e-1 * tau * vgp[d2] * tools.dNdx(ii, d2) * advgp[d1] * dvgpdx[1][d1] * tools.vol;
      Flocal(IW) -= 5e-1 * tau * vgp[d2] * tools.dNdx(ii, d2) * advgp[d1] * dvgpdx[2][d1] * tools.vol;
    }
  }

  // PSPG mass term
  Flocal(IP) += tau * tools.dNdx(ii, 0) * vgp[0] / dt * tools.vol;
  Flocal(IP) += tau * tools.dNdx(ii, 1) * vgp[1] / dt * tools.vol;
  Flocal(IP) += tau * tools.dNdx(ii, 2) * vgp[2] / dt * tools.vol;

  // PSPG advection term
  for (int d = 0; d < 3; d++)
  {
    Flocal(IP) -= 5e-1 * tau * tools.dNdx(ii, 0) * advgp[d] * dvgpdx[0][d] * tools.vol;
    Flocal(IP) -= 5e-1 * tau * tools.dNdx(ii, 1) * advgp[d] * dvgpdx[1][d] * tools.vol;
    Flocal(IP) -= 5e-1 * tau * tools.dNdx(ii, 2) * advgp[d] * dvgpdx[2][d] * tools.vol;
  }
}

void DirectProblem::setVelocityValue(MathTools3D &tools, const int ic, const int t)
{
  for (int d = 0; d < dim; d++)
  {
    advgp[d] = 0e0;
    vgp[d] = 0e0;
    for (int p = 0; p < grid.cell.nNodesInCell; p++)
    {
      if (t == 0)
      {
        advgp[d] += tools.N(p) * grid.node.v[grid.cell(ic).node[p]][d];
      }
      else
      {
        advgp[d] += tools.N(p) * (1.5 * grid.node.v[grid.cell(ic).node[p]][d] - 0.5 * grid.node.vPrev[grid.cell(ic).node[p]][d]);
      }
      vgp[d] += tools.N(p) * grid.node.v[grid.cell(ic).node[p]][d];
    }
  }
  for (int d = 0; d < dim; d++)
  {
    for (int e = 0; e < dim; e++)
    {
      dvgpdx[d][e] = 0e0;
      for (int p = 0; p < grid.cell.nNodesInCell; p++)
      {
        dvgpdx[d][e] += tools.dNdx(p, e) * grid.node.v[grid.cell(ic).node[p]][d];
      }
    }
  }
}
