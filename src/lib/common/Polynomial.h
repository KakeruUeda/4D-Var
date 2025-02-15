/**
 * @file Polynomial.h
 * @author K.Ueda
 * @date January, 2025
 */

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
using namespace Eigen;

class Polynomial3D
{
public:
  // Generate the design matrix for a 3rd-degree polynomial surface fit
  static MatrixXd generateDesignMatrix(const std::vector<double> &x, const std::vector<double> &y);

  // Perform 3D polynomial fitting
  static VectorXd fitPolynomial3D(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z);
  // Interpolate values on a new 2D grid
//   static void interpolateGrid(const VectorXd &coefficients, const std::vector<double> &mask, std::vector<double> &x_new,
//                               std::vector<double> &y_new, std::vector<double> &z_new, int new_nx, int new_ny);
};

#endif