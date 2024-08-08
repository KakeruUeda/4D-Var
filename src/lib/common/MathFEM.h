/**
 * @file MathFEM.h
 * @author K.Ueda
 * @date Jun, 2024
*/

#ifndef MATHFEM_H
#define MATHFEM_H

#include <iostream>
#include <vector>
#include "Define.h"
#include "MathCommon.h"
#include "ShapeFunction.h"
#include "Gauss.h"

class MathFEM
{   
public:
    static void comp_dxdr2D(double (&dxdr)[2][2], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &nNodesInCell);
    static void comp_dNdx2D(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[2][2], const int &nNodesInCell);
    static void comp_dxdr(double (&dxdr)[3][3], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &nNodesInCell);
    static void comp_dNdx(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[3][3], const int &nNodesInCell);
    static double comp_tau(std::vector<double> &vel, const double &he, const double &Re, const double &dt);
};

#endif
