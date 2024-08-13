/**
 * @file MathCommon.cpp
 * @ref  git@github.com:oubiomechlab/voxelFEMfluid.git
 */

#include "MathTool.h"

MathTools2D::MathTools2D(const int nNodesInCell)
{
  N.resize(nNodesInCell);
  xCurrent.resize(nNodesInCell, 2);
  dNdr.resize(nNodesInCell, 2);
  dNdx.resize(nNodesInCell, 2);
  K.resize(nNodesInCell, nNodesInCell);
}

MathTools3D::MathTools3D(const int nNodesInCell)
{
  N.resize(nNodesInCell);
  xCurrent.resize(nNodesInCell, 3);
  dNdr.resize(nNodesInCell, 3);
  dNdx.resize(nNodesInCell, 3);
  K.resize(nNodesInCell, nNodesInCell);
}

void MathTools2D::setZero()
{
  N.fillZero();
  xCurrent.fillZero();
  dNdr.fillZero();
  dNdx.fillZero();
  K.fillZero();
}

void MathTools3D::setZero()
{
  N.fillZero();
  xCurrent.fillZero();
  dNdr.fillZero();
  dNdx.fillZero();
  K.fillZero();
}

void MathTools2D::compInverseMatrix(double (&inv_a)[2][2], const double (&a)[2][2])
{
  double det;
  det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

  inv_a[0][0] = a[1][1] / det;
  inv_a[0][1] = -a[0][1] / det;
  inv_a[1][0] = -a[1][0] / det;
  inv_a[1][1] = a[0][0] / det;
}

void MathTools3D::compInverseMatrix(double (&inv_a)[3][3], const double (&a)[3][3])
{
  double det;
  det = compDeterminant(a);

  inv_a[0][0] = a[1][1] * a[2][2] - a[1][2] * a[2][1];
  inv_a[0][1] = a[0][2] * a[2][1] - a[0][1] * a[2][2];
  inv_a[0][2] = a[0][1] * a[1][2] - a[0][2] * a[1][1];
  inv_a[1][0] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
  inv_a[1][1] = a[0][0] * a[2][2] - a[0][2] * a[2][0];
  inv_a[1][2] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
  inv_a[2][0] = a[1][0] * a[2][1] - a[1][1] * a[2][0];
  inv_a[2][1] = a[0][1] * a[2][0] - a[0][0] * a[2][1];
  inv_a[2][2] = a[0][0] * a[1][1] - a[0][1] * a[1][0];

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++)
      inv_a[i][j] = inv_a[i][j] / det;
  }
}

double MathTools2D::compDeterminant(const double (&a)[2][2])
{
  double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
  return det;
}

double MathTools3D::compDeterminant(const double (&a)[3][3])
{
  double det = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2] -
               a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
  return det;
}

void MathTools2D::comp_dxdr(double (&dxdr)[2][2], Array2D<double> &dNdr, Array2D<double> &x1, const int &numOfNodeInElm)
{
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      dxdr[i][j] = 0e0;
      for(int p = 0; p < numOfNodeInElm; p++) {
        dxdr[i][j] += dNdr(p, j) * x1(p, i);
      }
    }
  }
}

void MathTools2D::comp_dNdx(Array2D<double> &dNdx, Array2D<double> &dNdr, const double (&dxdr)[2][2],
                            const int &numOfNodeInElm)
{
  double drdx[2][2];
  MathTools2D::compInverseMatrix(drdx, dxdr);

  for(int p = 0; p < numOfNodeInElm; p++) {
    for(int i = 0; i < 2; i++) {
      dNdx(p, i) = 0e0;
      for(int j = 0; j < 2; j++) {
        dNdx(p, i) += dNdr(p, j) * drdx[j][i];
      }
    }
  }
}

void MathTools3D::comp_dxdr(double (&dxdr)[3][3], Array2D<double> &dNdr, Array2D<double> &x1, const int &numOfNodeInElm)
{
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      dxdr[i][j] = 0e0;
      for(int p = 0; p < numOfNodeInElm; p++) {
        dxdr[i][j] += dNdr(p, j) * x1(p, i);
      }
    }
  }
}

void MathTools3D::comp_dNdx(Array2D<double> &dNdx, Array2D<double> &dNdr, const double (&dxdr)[3][3],
                            const int &numOfNodeInElm)
{
  double drdx[3][3];
  MathTools3D::compInverseMatrix(drdx, dxdr);

  for(int p = 0; p < numOfNodeInElm; p++) {
    for(int i = 0; i < 3; i++) {
      dNdx(p, i) = 0e0;
      for(int j = 0; j < 3; j++) {
        dNdx(p, i) += dNdr(p, j) * drdx[j][i];
      }
    }
  }
}

double MathTools2D::comp_tau(std::vector<double> &vel, const double &he, const double &Re, const double &dt)
{
  double tau = 0e0;
  double velMag = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

  double term1 = (2e0 / dt) * (2e0 / dt);
  double term2 = (2e0 * velMag / he) * (2e0 * velMag / he);
  double term3 = (4e0 / (Re * he * he)) * (4e0 / (Re * he * he));

  return tau = pow(term1 + term2 + term3, -5e-1);
}

