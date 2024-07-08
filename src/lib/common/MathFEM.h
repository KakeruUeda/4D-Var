#ifndef MATHFEM_H
#define MATHFEM_H

#include <iostream>
#include <vector>
#include "define.h"
#include "MathCommon.h"
#include "ShapeFunction.h"
#include "Gauss.h"

using namespace std;

class MathFEM
{   
    public:
        int IU, IV, IW, IP;
        int ILU, ILV, ILW;
        int JU, JV, JW, JP;
        int JLU, JLV, JLW;

        int ngp;

        double weight, detJ;
        double tau;

        double dxdr2D[2][2];
        double dxdr[3][3];

        std::vector<std::vector<int>> dofStart;
        std::vector<std::vector<int>> dofStartPlane;

        std::vector<double> vgp;
        std::vector<double> advgp;
        std::vector<std::vector<double>> dvgpdx;

        std::vector<double> N;
        std::vector<std::vector<double>> dNdr;
        std::vector<std::vector<double>> dNdx;
        std::vector<std::vector<double>> xCurrent;

        std::vector<double> N2D;
        std::vector<std::vector<double>> dNdr2D;
        std::vector<std::vector<double>> dNdx2D;
        std::vector<std::vector<double>> xCurrent2D;
 
        std::vector<std::vector<double>> K;

        void updateRowIndex(const int &ii, const int &ic);
        void updateColumnIndex(const int &jj, const int &ic);
        void updateRowIndexPlane(const int &ii, const int &ic);
        void updateColumnIndexPlane(const int &jj, const int &ic);

        static void comp_dxdr2D(double (&dxdr)[2][2], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &nNodesInCell);
        static void comp_dNdx2D(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[2][2], const int &nNodesInCell);
        static void comp_dxdr(double (&dxdr)[3][3], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &nNodesInCell);
        static void comp_dNdx(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[3][3], const int &nNodesInCell);
        static double comp_tau(std::vector<double> &vel, const double &he, const double &Re, const double &dt);
};

#endif
