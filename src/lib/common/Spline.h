/**
 * @file Spline.h
 * @ref https://qiita.com/khidaka/items/84610cd890ecb8443d96
 * @author K.Ueda
 * @date July, 2024
 */

#ifndef SPLINE_H
#define SPLINE_H

#include "petscksp.h"
#include "petscmat.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
using namespace Eigen;

class Spline
{
public:
  struct Coefficients
  {
    double a, b, c, d, x;
  };
  static std::vector<Coefficients> compCoefficients(const std::vector<double> &x, const std::vector<double> &y);
  static double evaluate(const std::vector<Coefficients> &coefficients, double x);
};

class Spline2D
{
public:
  using Coefficients = Spline::Coefficients;

  static std::vector<std::vector<Coefficients>> computeCoefficients(const std::vector<double> &x,
                                                                    const std::vector<double> &y,
                                                                    const std::vector<std::vector<double>> &z);
  static double evaluate(const std::vector<std::vector<Coefficients>> &coeffs_y_fixed, const std::vector<double> &x,
                         const std::vector<double> &y, double x_query, double y_query);
};

class Spline3D
{
public:
  using Coefficients = Spline::Coefficients;

  static std::vector<std::vector<std::vector<Coefficients>>>
  computeCoefficients(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                      const std::vector<std::vector<std::vector<double>>> &w);

  static double evaluate(const std::vector<std::vector<std::vector<Coefficients>>> &coeffs_z_fixed,
                         const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                         double x_query, double y_query, double z_query);
};

#endif