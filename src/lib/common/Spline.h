/**
 * @file ShapeFunction.h
 * @ref https://qiita.com/khidaka/items/84610cd890ecb8443d96
 * @author K.Ueda
 * @date July, 2024
*/

#ifndef SPLINE_H
#define SPLINE_H

#include <iostream>
#include <vector>
#include "petscksp.h"
#include "petscmat.h"
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;


class Spline
{
public:
    struct Coefficients{
        double a, b, c, d, x;
    };
    static std::vector<Coefficients> compCoefficients(const std::vector<double>& x, const std::vector<double>& y);
    static double evaluate(const std::vector<Coefficients>& coefficients, double x);
};

#endif