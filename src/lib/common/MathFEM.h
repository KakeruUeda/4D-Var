#include <iostream>
#include <vector>
#include "define.h"
#include "MathCommon.h"

using namespace std;

class MathFEM
{
    public:
        static void calc_dxdr2D(double (&dxdr)[2][2], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &nNodesInCell);
        static void calc_dNdx2D(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[2][2], const int &nNodesInCell);
        static void calc_dxdr(double (&dxdr)[3][3],std::vector<std::vector<double>> &dNdr,std::vector<std::vector<double>> &x1, const int &nNodesInCell);
        static void calc_dNdx(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[3][3], const int &nNodesInCell);
        static double calc_tau(const double (&vel)[3], double he, double Re, double dt);
};
