#include <omp.h>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

class MathCommon{
    public:
        static void calcInverseMatrix_2x2(double (&inv_a)[2][2],const double (&a)[2][2]);
        static void calcInverseMatrix_3x3(double (&inv_a)[3][3],const double (&a)[3][3]);
        static double calcDeterminant_3x3(const double (&a)[3][3]);
};