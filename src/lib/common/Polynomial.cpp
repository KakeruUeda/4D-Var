#include "Polynomial.h"

MatrixXd Polynomial3D::generateDesignMatrix(const std::vector<double> &x, const std::vector<double> &y)
{
  int n = x.size();
  MatrixXd designMatrix(n, 10);  // 10 terms for a 3rd-degree polynomial in two variables

  for(int i = 0; i < n; ++i) {
    double xi = x[i];
    double yi = y[i];

    designMatrix(i, 0) = 1.0;
    designMatrix(i, 1) = xi;
    designMatrix(i, 2) = yi;
    designMatrix(i, 3) = xi * xi;
    designMatrix(i, 4) = xi * yi;
    designMatrix(i, 5) = yi * yi;
    designMatrix(i, 6) = xi * xi * xi;
    designMatrix(i, 7) = xi * xi * yi;
    designMatrix(i, 8) = xi * yi * yi;
    designMatrix(i, 9) = yi * yi * yi;
  }

  return designMatrix;
}

VectorXd Polynomial3D::fitPolynomial3D(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z)
{
  MatrixXd designMatrix = generateDesignMatrix(x, y);
  VectorXd zVec = VectorXd::Map(z.data(), z.size());

  VectorXd coefficients = (designMatrix.transpose() * designMatrix).ldlt().solve(designMatrix.transpose() * zVec);

  return coefficients;
}

// void Polynomial3D::interpolateGrid(const VectorXd &coefficients, const std::vector<double> &mask,
//                                           std::vector<double> &x, std::vector<double> &y, std::vector<double> &z_new,
//                                           int new_nx, int new_ny)
// {

//   for(int j = 0; j < new_ny; ++j) {
//     for(int i = 0; i < new_nx; ++i) {
//       double x = i * dx_new + dx_new / 2;
//       double y = j * dy_new + dy_new / 2;

//       int mask_i = static_cast<int>(xi / dx);
//       int mask_j = static_cast<int>(yi / dz);

//       if(mask_i < nx && mask_j < nz) {
//         int mask_index = (mask_j * nx * ny) + mask_i;

//         if(mask[mask_index] != 0) {
//           double zi = coefficients[0] + coefficients[1] * xi + coefficients[2] * yi + coefficients[3] * xi * xi +
//                       coefficients[4] * xi * yi + coefficients[5] * yi * yi + coefficients[6] * xi * xi * xi +
//                       coefficients[7] * xi * xi * yi + coefficients[8] * xi * yi * yi + coefficients[9] * yi * yi * yi;

//           z_new.push_back(zi);
//         }
//       }
//     }
//   }
// }