/**
 * @file   Function.h
 * @author K.Ueda
 * @date   July, 2024
*/

#ifndef FUNCTION_H
#define FUNCTION_H

#include <iostream>
#include <vector>
#include "Tool.h"

class Function
{  
    public:
        Function(const int nNodesInCell, const int dim)
        {
            VecTool::resize(dxdr, dim, dim);
            VecTool::resize(xCurrent, nNodesInCell, dim);
            VecTool::resize(N, nNodesInCell);
            VecTool::resize(dNdr, nNodesInCell, dim);
            VecTool::resize(dNdx, nNodesInCell, dim);
            VecTool::resize(K, nNodesInCell, nNodesInCell);
        }
        double detJ, weight, vol;
        std::vector<std::vector<double>> dxdr;
        std::vector<std::vector<double>> xCurrent;
        std::vector<double> N;
        std::vector<std::vector<double>> dNdr;
        std::vector<std::vector<double>> dNdx;
        std::vector<std::vector<double>> K;
};

#endif