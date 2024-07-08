#include "MathFEM.h"

void MathFEM::updateRowIndex(const int &ii, const int &ic)
{
    IU = dofStart[ic][ii]; IV = IU + 1;  IW = IU + 2;
    IP = IU + 3; ILU = IU + 4; ILV = IU + 5; ILW = IU + 6;
}

void MathFEM::updateColumnIndex(const int &jj, const int &ic)
{
    JU = dofStart[ic][jj]; JV = JU + 1;  JW = JU + 2;
    JP = JU + 3; JLU = JU + 4; JLV = JU + 5; JLW = JU + 6;
}

void MathFEM::updateRowIndexPlane(const int &ii, const int &ic)
{
    IU = dofStartPlane[ic][ii]; IV = IU + 1; IW = IU + 2;
    IP = IU + 3; ILU = IU + 4; ILV = IU + 5; ILW = IU + 6;
}

void MathFEM::updateColumnIndexPlane(const int &jj, const int &ic)
{
    JU = dofStartPlane[ic][jj];  JV = JU + 1;  JW = JU + 2;
    JP = JU + 3; JLU = JU + 4; JLV = JU + 5; JLW = JU + 6;
}

void MathFEM::MathFEM::comp_dxdr2D(double (&dxdr)[2][2], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &numOfNodeInElm)
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

void MathFEM::MathFEM::comp_dNdx2D(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[2][2], const int &numOfNodeInElm)
{
    double drdx[2][2];
    MathCommon::compInverseMatrix_2x2(drdx, dxdr);

    for(int p=0; p<numOfNodeInElm; p++){
        for(int i=0; i<2; i++){
            dNdx[p][i] = 0e0;
            for(int j=0; j<2; j++){
                dNdx[p][i] += dNdr[p][j] * drdx[j][i];
            }
        }
    }
}

void MathFEM::MathFEM::comp_dxdr(double (&dxdr)[3][3], std::vector<std::vector<double>> &dNdr, std::vector<std::vector<double>> &x1, const int &numOfNodeInElm)
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

void MathFEM::comp_dNdx(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> &dNdr, const double (&dxdr)[3][3], const int &numOfNodeInElm)
{
    double drdx[3][3];
    MathCommon::compInverseMatrix_3x3(drdx, dxdr);

    for(int p=0;p<numOfNodeInElm;p++){
        for(int i=0; i<3; i++){
            dNdx[p][i] = 0e0;
            for(int j=0; j<3; j++){
                dNdx[p][i] += dNdr[p][j] * drdx[j][i];
            }
        }
    }
}

double MathFEM::comp_tau(std::vector<double> &vel, const double &he, const double &Re, const double &dt)
{
    double tau = 0e0;
    double velMag = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
    
    double term1 = (2e0 / dt) * (2e0 / dt);
    double term2 = (2e0 * velMag / he) * (2e0 * velMag / he);
    double term3 = (4e0 / (Re * he * he)) * (4e0 / (Re * he * he));
  
    return tau = pow(term1 + term2 + term3, -5e-1);
}
