/**
 * @file MathCommon.h
 * @ref  git@github.com:oubiomechlab/voxelFEMfluid.git
 */

#ifndef MATHCOMMON_H
#define MATHCOMMON_H

#include <string>
#include <iostream>
#include <cmath>

using namespace std;

class MathCommon{
    public:
        static void compostinverseMatrix_2x2(double (&inv_a)[2][2], const double (&a)[2][2]);
        static void compostinverseMatrix_3x3(double (&inv_a)[3][3], const double (&a)[3][3]);
        static double compDeterminant_2x2(const double (&a)[2][2]);
        static double compDeterminant_3x3(const double (&a)[3][3]);
};

#endif