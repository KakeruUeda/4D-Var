/**
 * @file MathTool.h
 * @ref  git@github.com:oubiomechlab/voxelFEMfluid.git
 */

#ifndef MATHTOOL_H
#define MATHTOOL_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "Array.h"

using namespace std;

class MathTools2D
{
public:
  double detJ, weight, vol;
  double dxdr[2][2];
  Array2D<double> xCurrent;
  Array1D<double> N;
  Array2D<double> dNdr;
  Array2D<double> dNdx;
  Array2D<double> K;

  MathTools2D(const int nNodesInCell)
  {
    N.resize(nNodesInCell);
    xCurrent.resize(nNodesInCell, 2);
    dNdr.resize(nNodesInCell, 2);
    dNdx.resize(nNodesInCell, 2);
    K.resize(nNodesInCell, 2);
  }

  static void compInverseMatrix(double (&inv_a)[2][2], const double (&a)[2][2]);
  static double compDeterminant(const double (&a)[2][2]);

  static void comp_dxdr(double (&dxdr)[2][2], Array2D<double> &dNdr, Array2D<double> &x1, const int &nNodesInCell);
  static void comp_dNdx(Array2D<double> &dNdx, Array2D<double> &dNdr, const double (&dxdr)[2][2], const int &nNodesInCell);
  static double comp_tau(std::vector<double> &vel, const double &he, const double &Re, const double &dt);
};


class MathTools3D
{
public:
  double detJ, weight, vol;
  double dxdr[3][3];
  Array2D<double> xCurrent;
  Array1D<double> N;
  Array2D<double> dNdr;
  Array2D<double> dNdx;
  Array2D<double> K;

  MathTools3D(const int nNodesInCell)
  {
    N.resize(nNodesInCell);
    xCurrent.resize(nNodesInCell, 3);
    dNdr.resize(nNodesInCell, 3);
    dNdx.resize(nNodesInCell, 3);
    K.resize(nNodesInCell, 3);
  }

  static void compInverseMatrix(double (&inv_a)[3][3], const double (&a)[3][3]);
  static double compDeterminant(const double (&a)[3][3]);

  static void comp_dxdr(double (&dxdr)[3][3], Array2D<double>&dNdr, Array2D<double> &x1, const int &nNodesInCell);
  static void comp_dNdx(Array2D<double> &dNdx, Array2D<double> &dNdr, const double (&dxdr)[3][3], const int &nNodesInCell);
  static double comp_tau(std::vector<double> &vel, const double &he, const double &Re, const double &dt);
};

#endif