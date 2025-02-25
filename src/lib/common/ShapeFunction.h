/**
 * @file ShapeFunction.h
 * @ref git@github.com:oubiomechlab/voxelFEMfluid.git
 */

#ifndef SHAPEFUNCTION_H
#define SHAPEFUNCTION_H

#include "Array.h"
#include <iostream>
#include <vector>
using namespace std;

class ShapeFunction1D
{
public:
  static void P2_N(std::vector<double> &N, const double &g1)
  {
    N[0] = 5e-1 * (1e0 - g1);
    N[1] = 5e-1 * (1e0 + g1);
  }
  static void P2_dNdr(std::vector<double> &dNdr, const double &g1)
  {
    dNdr[0] = -5e-1;
    dNdr[1] = 5e-1;
  }

  static void P3_N(std::vector<double> &N, const double &g1)
  {
    N[0] = -5e-1 * g1 * (1e0 - g1);
    N[1] = 5e-1 * g1 * (1e0 + g1);
    N[2] = 1e0 - g1 * g1;
  }
  static void P3_dNdr(std::vector<double> &dNdr, const double &g1)
  {
    dNdr[0] = g1 - 5e-1;
    dNdr[1] = g1 + 5e-1;
    dNdr[2] = -2e0 * g1;
  }

  static void P2_N(Array1D<double> &N, const double &g1)
  {
    N(0) = 5e-1 * (1e0 - g1);
    N(1) = 5e-1 * (1e0 + g1);
  }
  static void P2_dNdr(Array1D<double> &dNdr, const double &g1)
  {
    dNdr(0) = -5e-1;
    dNdr(1) = 5e-1;
  }

  static void P3_N(Array1D<double> &N, const double &g1)
  {
    N(0) = -5e-1 * g1 * (1e0 - g1);
    N(1) = 5e-1 * g1 * (1e0 + g1);
    N(2) = 1e0 - g1 * g1;
  }
  static void P3_dNdr(Array1D<double> &dNdr, const double &g1)
  {
    dNdr(0) = g1 - 5e-1;
    dNdr(1) = g1 + 5e-1;
    dNdr(2) = -2e0 * g1;
  }
};

class ShapeFunction2D
{
public:
  static void C2D3_N(std::vector<double> &N, const double &L1, const double &L2, const double &L3)
  {
    N[0] = L1;
    N[1] = L2;
    N[2] = L3;
  }

  static void C2D4_N(std::vector<double> &N, const double &g1, const double &g2)
  {
    N[0] = 2.5e-1 * (1e+0 - g1) * (1e+0 - g2);
    N[1] = 2.5e-1 * (1e+0 + g1) * (1e+0 - g2);
    N[2] = 2.5e-1 * (1e+0 + g1) * (1e+0 + g2);
    N[3] = 2.5e-1 * (1e+0 - g1) * (1e+0 + g2);
  }

  static void C2D6_N(std::vector<double> &N, const double &L1, const double &L2, const double &L3)
  {
    N[0] = L1 * (2e0 * L1 - 1e0);
    N[1] = L2 * (2e0 * L2 - 1e0);
    N[2] = L3 * (2e0 * L3 - 1e0);
    N[3] = 4e0 * L1 * L2;
    N[4] = 4e0 * L2 * L3;
    N[5] = 4e0 * L1 * L3;
  }

  static void C2D8_N(std::vector<double> &N, const double &g1, const double &g2)
  {
    N[0] = 2.5e-1 * (1e+0 - g1) * (1e+0 - g2) * (-1e+0 - g1 - g2);
    N[1] = 2.5e-1 * (1e+0 + g1) * (1e+0 - g2) * (-1e+0 + g1 - g2);
    N[2] = 2.5e-1 * (1e+0 + g1) * (1e+0 + g2) * (-1e+0 + g1 + g2);
    N[3] = 2.5e-1 * (1e+0 - g1) * (1e+0 + g2) * (-1e+0 - g1 + g2);
    N[4] = 5e-1 * (1e+0 - g1 * g1) * (1e+0 - g2);
    N[5] = 5e-1 * (1e+0 + g1) * (1e+0 - g2 * g2);
    N[6] = 5e-1 * (1e+0 - g1 * g1) * (1e+0 + g2);
    N[7] = 5e-1 * (1e+0 - g1) * (1e+0 - g2 * g2);
  }

  static void C2D3_dNdr(std::vector<std::vector<double>> &dNdr, const double &L1, const double &L2, const double &L3)
  {
    dNdr[0][0] = -1e0;
    dNdr[0][1] = -1e0;
    dNdr[1][0] = 1e0;
    dNdr[1][1] = 0e0;
    dNdr[2][0] = 0e0;
    dNdr[2][1] = 1e0;
  }

  static void C2D4_dNdr(std::vector<std::vector<double>> &dNdr, const double &g1, const double &g2)
  {
    dNdr[0][0] = -2.5e-1 * (1e+0 - g2);
    dNdr[0][1] = -2.5e-1 * (1e+0 - g1);
    dNdr[1][0] = 2.5e-1 * (1e+0 - g2);
    dNdr[1][1] = -2.5e-1 * (1e+0 + g1);
    dNdr[2][0] = 2.5e-1 * (1e+0 + g2);
    dNdr[2][1] = 2.5e-1 * (1e+0 + g1);
    dNdr[3][0] = -2.5e-1 * (1e+0 + g2);
    dNdr[3][1] = 2.5e-1 * (1e+0 - g1);
  }

  static void C2D6_dNdr(std::vector<std::vector<double>> &dNdr, const double &L1, const double &L2, const double &L3)
  {
    dNdr[0][0] = -4e0 * L1 + 1e0;
    dNdr[0][1] = -4e0 * L1 + 1e0;
    dNdr[1][0] = 4e0 * L2 - 1e0;
    dNdr[1][1] = 0e0;
    dNdr[2][0] = 0e0;
    dNdr[2][1] = 4e0 * L3 - 1e0;
    dNdr[3][0] = 4e0 * (L1 - L2);
    dNdr[3][1] = -4e0 * L2;
    dNdr[4][0] = 4e0 * L3;
    dNdr[4][1] = 4e0 * L2;
    dNdr[5][0] = -4e0 * L3;
    dNdr[5][1] = 4e0 * (L1 - L3);
  }

  static void C2D3_N(Array1D<double> &N, const double &L1, const double &L2, const double &L3)
  {
    N(0) = L1;
    N(1) = L2;
    N(2) = L3;
  }

  static void C2D4_N(Array1D<double> &N, const double &g1, const double &g2)
  {
    N(0) = 2.5e-1 * (1e+0 - g1) * (1e+0 - g2);
    N(1) = 2.5e-1 * (1e+0 + g1) * (1e+0 - g2);
    N(2) = 2.5e-1 * (1e+0 + g1) * (1e+0 + g2);
    N(3) = 2.5e-1 * (1e+0 - g1) * (1e+0 + g2);
  }

  static void C2D6_N(Array1D<double> &N, const double &L1, const double &L2, const double &L3)
  {
    N(0) = L1 * (2e0 * L1 - 1e0);
    N(1) = L2 * (2e0 * L2 - 1e0);
    N(2) = L3 * (2e0 * L3 - 1e0);
    N(3) = 4e0 * L1 * L2;
    N(4) = 4e0 * L2 * L3;
    N(5) = 4e0 * L1 * L3;
  }

  static void C2D8_N(Array1D<double> &N, const double &g1, const double &g2)
  {
    N(0) = 2.5e-1 * (1e+0 - g1) * (1e+0 - g2) * (-1e+0 - g1 - g2);
    N(1) = 2.5e-1 * (1e+0 + g1) * (1e+0 - g2) * (-1e+0 + g1 - g2);
    N(2) = 2.5e-1 * (1e+0 + g1) * (1e+0 + g2) * (-1e+0 + g1 + g2);
    N(3) = 2.5e-1 * (1e+0 - g1) * (1e+0 + g2) * (-1e+0 - g1 + g2);
    N(4) = 5e-1 * (1e+0 - g1 * g1) * (1e+0 - g2);
    N(5) = 5e-1 * (1e+0 + g1) * (1e+0 - g2 * g2);
    N(6) = 5e-1 * (1e+0 - g1 * g1) * (1e+0 + g2);
    N(7) = 5e-1 * (1e+0 - g1) * (1e+0 - g2 * g2);
  }

  static void C2D3_dNdr(Array2D<double> &dNdr, const double &L1, const double &L2, const double &L3)
  {
    dNdr(0, 0) = -1e0;
    dNdr(0, 1) = -1e0;
    dNdr(1, 0) = 1e0;
    dNdr(1, 1) = 0e0;
    dNdr(2, 0) = 0e0;
    dNdr(2, 1) = 1e0;
  }

  static void C2D4_dNdr(Array2D<double> &dNdr, const double &g1, const double &g2)
  {
    dNdr(0, 0) = -2.5e-1 * (1e+0 - g2);
    dNdr(0, 1) = -2.5e-1 * (1e+0 - g1);
    dNdr(1, 0) = 2.5e-1 * (1e+0 - g2);
    dNdr(1, 1) = -2.5e-1 * (1e+0 + g1);
    dNdr(2, 0) = 2.5e-1 * (1e+0 + g2);
    dNdr(2, 1) = 2.5e-1 * (1e+0 + g1);
    dNdr(3, 0) = -2.5e-1 * (1e+0 + g2);
    dNdr(3, 1) = 2.5e-1 * (1e+0 - g1);
  }

  static void C2D6_dNdr(Array2D<double> &dNdr, const double &L1, const double &L2, const double &L3)
  {
    dNdr(0, 0) = -4e0 * L1 + 1e0;
    dNdr(0, 1) = -4e0 * L1 + 1e0;
    dNdr(1, 0) = 4e0 * L2 - 1e0;
    dNdr(1, 1) = 0e0;
    dNdr(2, 0) = 0e0;
    dNdr(2, 1) = 4e0 * L3 - 1e0;
    dNdr(3, 0) = 4e0 * (L1 - L2);
    dNdr(3, 1) = -4e0 * L2;
    dNdr(4, 0) = 4e0 * L3;
    dNdr(4, 1) = 4e0 * L2;
    dNdr(5, 0) = -4e0 * L3;
    dNdr(5, 1) = 4e0 * (L1 - L3);
  }
};

class ShapeFunction3D
{
public:
  static void C3D8_N(std::vector<double> &N, const double &g1, const double &g2, const double &g3)
  {
    N[0] = 1.25e-1 * (1e0 - g1) * (1e0 - g2) * (1e0 - g3);
    N[1] = 1.25e-1 * (1e0 + g1) * (1e0 - g2) * (1e0 - g3);
    N[2] = 1.25e-1 * (1e0 + g1) * (1e0 + g2) * (1e0 - g3);
    N[3] = 1.25e-1 * (1e0 - g1) * (1e0 + g2) * (1e0 - g3);
    N[4] = 1.25e-1 * (1e0 - g1) * (1e0 - g2) * (1e0 + g3);
    N[5] = 1.25e-1 * (1e0 + g1) * (1e0 - g2) * (1e0 + g3);
    N[6] = 1.25e-1 * (1e0 + g1) * (1e0 + g2) * (1e0 + g3);
    N[7] = 1.25e-1 * (1e0 - g1) * (1e0 + g2) * (1e0 + g3);
  }

  static void C3D8_dNdr(std::vector<std::vector<double>> &dNdr, const double &g1, const double &g2, const double &g3)
  {
    dNdr[0][0] = -1.25e-1 * (1e0 - g2) * (1e0 - g3);
    dNdr[0][1] = -1.25e-1 * (1e0 - g1) * (1e0 - g3);
    dNdr[0][2] = -1.25e-1 * (1e0 - g1) * (1e0 - g2);
    dNdr[1][0] = 1.25e-1 * (1e0 - g2) * (1e0 - g3);
    dNdr[1][1] = -1.25e-1 * (1e0 + g1) * (1e0 - g3);
    dNdr[1][2] = -1.25e-1 * (1e0 + g1) * (1e0 - g2);
    dNdr[2][0] = 1.25e-1 * (1e0 + g2) * (1e0 - g3);
    dNdr[2][1] = 1.25e-1 * (1e0 + g1) * (1e0 - g3);
    dNdr[2][2] = -1.25e-1 * (1e0 + g1) * (1e0 + g2);
    dNdr[3][0] = -1.25e-1 * (1e0 + g2) * (1e0 - g3);
    dNdr[3][1] = 1.25e-1 * (1e0 - g1) * (1e0 - g3);
    dNdr[3][2] = -1.25e-1 * (1e0 - g1) * (1e0 + g2);
    dNdr[4][0] = -1.25e-1 * (1e0 - g2) * (1e0 + g3);
    dNdr[4][1] = -1.25e-1 * (1e0 - g1) * (1e0 + g3);
    dNdr[4][2] = 1.25e-1 * (1e0 - g1) * (1e0 - g2);
    dNdr[5][0] = 1.25e-1 * (1e0 - g2) * (1e0 + g3);
    dNdr[5][1] = -1.25e-1 * (1e0 + g1) * (1e0 + g3);
    dNdr[5][2] = 1.25e-1 * (1e0 + g1) * (1e0 - g2);
    dNdr[6][0] = 1.25e-1 * (1e0 + g2) * (1e0 + g3);
    dNdr[6][1] = 1.25e-1 * (1e0 + g1) * (1e0 + g3);
    dNdr[6][2] = 1.25e-1 * (1e0 + g1) * (1e0 + g2);
    dNdr[7][0] = -1.25e-1 * (1e0 + g2) * (1e0 + g3);
    dNdr[7][1] = 1.25e-1 * (1e0 - g1) * (1e0 + g3);
    dNdr[7][2] = 1.25e-1 * (1e0 - g1) * (1e0 + g2);
  }

  static void C3D8_N(Array1D<double> &N, const double &g1, const double &g2, const double &g3)
  {
    N(0) = 1.25e-1 * (1e0 - g1) * (1e0 - g2) * (1e0 - g3);
    N(1) = 1.25e-1 * (1e0 + g1) * (1e0 - g2) * (1e0 - g3);
    N(2) = 1.25e-1 * (1e0 + g1) * (1e0 + g2) * (1e0 - g3);
    N(3) = 1.25e-1 * (1e0 - g1) * (1e0 + g2) * (1e0 - g3);
    N(4) = 1.25e-1 * (1e0 - g1) * (1e0 - g2) * (1e0 + g3);
    N(5) = 1.25e-1 * (1e0 + g1) * (1e0 - g2) * (1e0 + g3);
    N(6) = 1.25e-1 * (1e0 + g1) * (1e0 + g2) * (1e0 + g3);
    N(7) = 1.25e-1 * (1e0 - g1) * (1e0 + g2) * (1e0 + g3);
  }

  static void C3D8_dNdr(Array2D<double> &dNdr, const double &g1, const double &g2, const double &g3)
  {
    dNdr(0, 0) = -1.25e-1 * (1e0 - g2) * (1e0 - g3);
    dNdr(0, 1) = -1.25e-1 * (1e0 - g1) * (1e0 - g3);
    dNdr(0, 2) = -1.25e-1 * (1e0 - g1) * (1e0 - g2);
    dNdr(1, 0) = 1.25e-1 * (1e0 - g2) * (1e0 - g3);
    dNdr(1, 1) = -1.25e-1 * (1e0 + g1) * (1e0 - g3);
    dNdr(1, 2) = -1.25e-1 * (1e0 + g1) * (1e0 - g2);
    dNdr(2, 0) = 1.25e-1 * (1e0 + g2) * (1e0 - g3);
    dNdr(2, 1) = 1.25e-1 * (1e0 + g1) * (1e0 - g3);
    dNdr(2, 2) = -1.25e-1 * (1e0 + g1) * (1e0 + g2);
    dNdr(3, 0) = -1.25e-1 * (1e0 + g2) * (1e0 - g3);
    dNdr(3, 1) = 1.25e-1 * (1e0 - g1) * (1e0 - g3);
    dNdr(3, 2) = -1.25e-1 * (1e0 - g1) * (1e0 + g2);
    dNdr(4, 0) = -1.25e-1 * (1e0 - g2) * (1e0 + g3);
    dNdr(4, 1) = -1.25e-1 * (1e0 - g1) * (1e0 + g3);
    dNdr(4, 2) = 1.25e-1 * (1e0 - g1) * (1e0 - g2);
    dNdr(5, 0) = 1.25e-1 * (1e0 - g2) * (1e0 + g3);
    dNdr(5, 1) = -1.25e-1 * (1e0 + g1) * (1e0 + g3);
    dNdr(5, 2) = 1.25e-1 * (1e0 + g1) * (1e0 - g2);
    dNdr(6, 0) = 1.25e-1 * (1e0 + g2) * (1e0 + g3);
    dNdr(6, 1) = 1.25e-1 * (1e0 + g1) * (1e0 + g3);
    dNdr(6, 2) = 1.25e-1 * (1e0 + g1) * (1e0 + g2);
    dNdr(7, 0) = -1.25e-1 * (1e0 + g2) * (1e0 + g3);
    dNdr(7, 1) = 1.25e-1 * (1e0 - g1) * (1e0 + g3);
    dNdr(7, 2) = 1.25e-1 * (1e0 - g1) * (1e0 + g2);
  }
};

#endif