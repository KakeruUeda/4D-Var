#include "MathFEM.h"

void MathFEM::calc_dxdr2D(double (&dxdr)[2][2], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &numOfNodeInElm)
{
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      dxdr[i][j] = 0e0;
      for(int p=0;p<numOfNodeInElm;p++){
        dxdr[i][j] += dNdr[p][j] * x1[p][i];
      }
    }
  }
}

void MathFEM::calc_dNdx2D(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[2][2], const int &numOfNodeInElm)
{
  double drdx[2][2];
  MathCommon::calcInverseMatrix_2x2(drdx, dxdr);

  for(int p=0; p<numOfNodeInElm; p++){
      for(int i=0; i<2; i++){
          dNdx[p][i] = 0e0;
          for(int j=0; j<2; j++) dNdx[p][i] += dNdr[p][j] * drdx[j][i];
      }
  }
}

void MathFEM::calc_dxdr(double (&dxdr)[3][3], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &numOfNodeInElm)
{
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      dxdr[i][j] = 0e0;
      for(int p=0; p<numOfNodeInElm; p++){
        dxdr[i][j] += dNdr[p][j] * x1[p][i];
      }
    }
  }
}

void MathFEM::calc_dNdx(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[3][3], const int &numOfNodeInElm)
{
  double drdx[3][3];
  MathCommon::calcInverseMatrix_3x3(drdx, dxdr);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0; i<3; i++){
      dNdx[p][i] = 0e0;
      for(int j=0; j<3; j++) dNdx[p][i] += dNdr[p][j] * drdx[j][i];
    }
  }
}

double MathFEM::calc_tau(const double (&vel)[3], double he, double Re, double dt)
{
    double tau = 0e0;
    double velMag = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
    
    double term1 = (2e0 / dt) * (2e0 / dt);
    double term2 = (2e0 * velMag / he) * (2e0 * velMag / he);
    double term3 = (4e0 / (Re * he * he)) * (4e0 / (Re * he * he));
  
    return tau = pow(term1 + term2 + term3, -5e-1);
}
