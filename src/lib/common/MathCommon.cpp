/**
 * @file MathCommon.cpp
 * @ref  git@github.com:oubiomechlab/voxelFEMfluid.git
 */

#include "MathCommon.h"

void MathCommon::compInverseMatrix_2x2(double (&inv_a)[2][2], const double (&a)[2][2])
{   
    double det;
    det = a[0][0]*a[1][1]-a[0][1]*a[1][0];

    inv_a[0][0] =  a[1][1]/det;
    inv_a[0][1] = -a[0][1]/det;
    inv_a[1][0] = -a[1][0]/det;
    inv_a[1][1] =  a[0][0]/det;
}

void MathCommon::compInverseMatrix_3x3(double (&inv_a)[3][3], const double (&a)[3][3])
{
    double det;
    det = compDeterminant_3x3(a);

    inv_a[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
    inv_a[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
    inv_a[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
    inv_a[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
    inv_a[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
    inv_a[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
    inv_a[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
    inv_a[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
    inv_a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++) inv_a[i][j] = inv_a[i][j] / det;
    }
}

double MathCommon::compDeterminant_2x2(const double (&a)[2][2])
{
    double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    return det;
}

double MathCommon::compDeterminant_3x3(const double (&a)[3][3])
{
    double det  = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2]
                - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
    return det;
}