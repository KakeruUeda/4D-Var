#include "Spline.h"

/**
 * @brief Compute coefficients of the spline
 * @ref   https://qiita.com/khidaka/items/84610cd890ecb8443d96
 */
std::vector<Spline::Coefficients> Spline::compCoefficients(const std::vector<double> &x, const std::vector<double> &y)
{
  int n = x.size() - 1;
  std::vector<double> h(n);
  for(int i = 0; i < n; ++i) {
    h[i] = x[i + 1] - x[i];
  }

  MatrixXd A = MatrixXd::Zero(4 * n, 4 * n);
  VectorXd b = VectorXd::Zero(4 * n);

  int row = 0;

  for(int i = 0; i < n; ++i) {
    A(row, 4 * i) = 1;
    b(row) = y[i];
    row++;
    A(row, 4 * i) = 1;
    A(row, 4 * i + 1) = h[i];
    A(row, 4 * i + 2) = h[i] * h[i];
    A(row, 4 * i + 3) = h[i] * h[i] * h[i];
    b(row) = y[i + 1];
    row++;
  }

  for(int i = 1; i < n; ++i) {
    A(row, 4 * (i - 1) + 1) = 1;
    A(row, 4 * (i - 1) + 2) = 2 * h[i - 1];
    A(row, 4 * (i - 1) + 3) = 3 * h[i - 1] * h[i - 1];
    A(row, 4 * i + 1) = -1;
    row++;
  }

  for(int i = 1; i < n; ++i) {
    A(row, 4 * (i - 1) + 2) = 2;
    A(row, 4 * (i - 1) + 3) = 6 * h[i - 1];
    A(row, 4 * i + 2) = -2;
    row++;
  }

  A(row, 2) = 2;
  row++;
  A(row, 4 * (n - 1) + 2) = 2;
  A(row, 4 * (n - 1) + 3) = 6 * h[n - 1];
  row++;

  VectorXd coeffs = A.lu().solve(b);

  std::vector<Spline::Coefficients> coefficients(n);
  for(int i = 0; i < n; ++i) {
    coefficients[i].a = coeffs(4 * i);
    coefficients[i].b = coeffs(4 * i + 1);
    coefficients[i].c = coeffs(4 * i + 2);
    coefficients[i].d = coeffs(4 * i + 3);
    coefficients[i].x = x[i];
  }

  return coefficients;
}

/**
 * @brief Evaluate the spline
 * @ref   https://qiita.com/khidaka/items/84610cd890ecb8443d96
 */
double Spline::evaluate(const std::vector<Coefficients> &coefficients, double x)
{
  if(coefficients.empty()) {
    PetscPrintf(PETSC_COMM_WORLD, "Spline coefficients are empty.\n");
  }

  const Coefficients *spline = nullptr;
  for(const auto &s : coefficients) {
    if(x >= s.x) {
      spline = &s;
    } else {
      break;
    }
  }

  if(!spline) {
    PetscPrintf(PETSC_COMM_WORLD, "x is out of the range of the spline.\n");
  }

  double dx = x - spline->x;
  return spline->a + spline->b * dx + spline->c * dx * dx + spline->d * dx * dx * dx;
}

/**
 * @brief 2D spline interpolation
 */
std::vector<std::vector<Spline2D::Coefficients>>
Spline2D::computeCoefficients(const std::vector<double> &x, const std::vector<double> &y,
                              const std::vector<std::vector<double>> &z)
{
  int ny = y.size();
  int nx = x.size();

  std::vector<std::vector<Coefficients>> coeffs_y_fixed(ny);
  for(int i = 0; i < ny; ++i) {
    std::vector<double> z_i(nx);
    for(int j = 0; j < nx; ++j) {
      z_i[j] = z[i][j];
    }
    coeffs_y_fixed[i] = Spline::compCoefficients(x, z_i);
  }

  return coeffs_y_fixed;
}

double Spline2D::evaluate(const std::vector<std::vector<Coefficients>> &coeffs_y_fixed, const std::vector<double> &x,
                          const std::vector<double> &y, double x_query, double y_query)
{
  std::vector<double> z_interpolated(y.size());
  for(int i = 0; i < y.size(); ++i) {
    z_interpolated[i] = Spline::evaluate(coeffs_y_fixed[i], x_query);
  }

  std::vector<Coefficients> coeffs_x_fixed = Spline::compCoefficients(y, z_interpolated);
  return Spline::evaluate(coeffs_x_fixed, y_query);
}

std::vector<std::vector<std::vector<Spline::Coefficients>>>
Spline3D::computeCoefficients(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                              const std::vector<std::vector<std::vector<double>>> &w)
{
  int nz = z.size();
  std::vector<std::vector<std::vector<Coefficients>>> coeffs_z_fixed(nz);

  for(int k = 0; k < nz; ++k) {
    std::vector<std::vector<double>> w_k(y.size(), std::vector<double>(x.size()));
    for(int i = 0; i < y.size(); ++i) {
      for(int j = 0; j < x.size(); ++j) {
        w_k[i][j] = w[k][i][j];
      }
    }
    coeffs_z_fixed[k] = Spline2D::computeCoefficients(x, y, w_k);
  }

  return coeffs_z_fixed;
}

double Spline3D::evaluate(const std::vector<std::vector<std::vector<Coefficients>>> &coeffs_z_fixed,
                          const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                          double x_query, double y_query, double z_query)
{
  std::vector<double> w_interpolated(z.size());
  for(int k = 0; k < z.size(); ++k) {
    w_interpolated[k] = Spline2D::evaluate(coeffs_z_fixed[k], x, y, x_query, y_query);
  }

  std::vector<Coefficients> coeffs_z = Spline::compCoefficients(z, w_interpolated);
  return Spline::evaluate(coeffs_z, z_query);
}
