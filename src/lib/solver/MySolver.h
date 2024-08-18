/**
 * @file Initialize.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#ifndef MYSOLVER_H
#define MYSOLVER_H

#include <functional>
#include <iostream>
#include <vector>

// Central Difference Method
class CDM
{
public:
  static double xDerivative(const std::function<double(int, int, int)> &f, int i, int j, int k, double dx)
  {
    return (f(i + 1, j, k) - f(i - 1, j, k)) / (2.0 * dx);
  }

  static double yDerivative(const std::function<double(int, int, int)> &f, int i, int j, int k, double dy)
  {
    return (f(i, j + 1, k) - f(i, j - 1, k)) / (2.0 * dy);
  }

  static double zDerivative(const std::function<double(int, int, int)> &f, int i, int j, int k, double dz)
  {
    return (f(i, j, k + 1) - f(i, j, k - 1)) / (2.0 * dz);
  }
};

#endif  // MYSOLVER_H
