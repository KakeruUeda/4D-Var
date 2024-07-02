#include "DirectProblem.h"

void DirectProblem::matrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, 
                                       const int ic, const int t)
{
    int n1, n2, n3;
    int IU,IV,IW,IP;
    int JU,JV,JW,JP;
    int nGaussPoint = 2;

    std::vector<double> N;
    std::vector<std::vector<double>> xCurrent;
    std::vector<std::vector<double>> dNdr;
    std::vector<std::vector<double>> dNdx;
    std::vector<std::vector<double>> K;

    VecTool::resize(xCurrent, grid.cell.nNodesInCell, dim);
    VecTool::resize(N, grid.cell.nNodesInCell);
    VecTool::resize(dNdr, grid.cell.nNodesInCell, dim);
    VecTool::resize(dNdx, grid.cell.nNodesInCell, dim);
    VecTool::resize(K, grid.cell.nNodesInCell, grid.cell.nNodesInCell);
    
    for(int p=0; p<grid.cell.nNodesInCell; p++)
        for(int d=0; d<dim; d++)
            xCurrent[p][d] = grid.node.x[grid.cell(ic).node[p]][d];
            
    double he = fabs(xCurrent[1][0] - xCurrent[0][0]);
    double f = resistance * alpha * (1e0 - grid.cell(ic).phi) / (alpha + grid.cell(ic).phi);
    double detJ, weight;
    Gauss gauss(nGaussPoint);

    double dxdr[3][3] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    for(int i1=0; i1<nGaussPoint; i1++){
        for(int i2=0; i2<nGaussPoint; i2++){
            for(int i3=0; i3<nGaussPoint; i3++){
                ShapeFunction3D::C3D8_N(N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                ShapeFunction3D::C3D8_dNdr(dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                MathFEM::calc_dxdr(dxdr, dNdr, xCurrent, grid.cell.nNodesInCell);
                detJ = MathCommon::calcDeterminant_3x3(dxdr);
                weight = gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];

                MathFEM::calc_dNdx(dNdx, dNdr, dxdr, grid.cell.nNodesInCell);
                MathFEM::calc_dNdx(dNdx, dNdr, dxdr, grid.cell.nNodesInCell);

                double vel[3] = {0e0, 0e0, 0e0};
                double advel[3] {0e0, 0e0, 0e0};
                double dvdx[3][3] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

                setVelocityValue(vel, advel, dvdx, N, dNdx, ic, t);
                double tau = MathFEM::calc_tau(advel, he, Re, dt);

                for(int ii=0; ii<grid.cell.nNodesInCell; ii++){  
                    IU = 4 * ii;
                    IV = IU + 1;
                    IW = IU + 2;
                    IP = IU + 3;
                    for(int jj=0; jj<grid.cell.nNodesInCell; jj++){  
                        JU = 4 * jj;
                        JV = JU + 1;
                        JW = JU + 2;
                        JP = JU + 3;
                        K[ii][jj] = 0e0;

                        for(int k=0;k<3;k++){
                            K[ii][jj] += dNdx[ii][k] *dNdx[jj][k];
                        }
                        // MASS TERM
                        Klocal(IU, JU) += N[ii] * N[jj] / dt * detJ * weight;
                        Klocal(IV, JV) += N[ii] * N[jj] / dt * detJ * weight;
                        Klocal(IW, JW) += N[ii] * N[jj] / dt * detJ * weight;
                         // DIFFUSION TERM
                        for(int mm=0; mm<3; mm++){
                            if(mm == 0){n1 = 2e0; n2 = 1e0; n3 = 1e0;}
                            if(mm == 1){n1 = 1e0; n2 = 2e0; n3 = 1e0;}
                            if(mm == 2){n1 = 1e0; n2 = 1e0; n3 = 2e0;}
                            Klocal(IU, JU) += 5e-1 * n1 * dNdx[ii][mm] * dNdx[jj][mm] / Re * detJ * weight;
                            Klocal(IV, JV) += 5e-1 * n2 * dNdx[ii][mm] * dNdx[jj][mm] / Re * detJ * weight;
                            Klocal(IW, JW) += 5e-1 * n3 * dNdx[ii][mm] * dNdx[jj][mm] / Re * detJ * weight;
                        }
                        Klocal(IU, JV) += 5e-1 * dNdx[ii][1] * dNdx[jj][0] / Re * detJ * weight;
                        Klocal(IU, JW) += 5e-1 * dNdx[ii][2] * dNdx[jj][0] / Re * detJ * weight;
                        Klocal(IV, JU) += 5e-1 * dNdx[ii][0] * dNdx[jj][1] / Re * detJ * weight;
                        Klocal(IV, JW) += 5e-1 * dNdx[ii][2] * dNdx[jj][1] / Re * detJ * weight;
                        Klocal(IW, JU) += 5e-1 * dNdx[ii][0] * dNdx[jj][2] / Re * detJ * weight;
                        Klocal(IW, JV) += 5e-1 * dNdx[ii][1] * dNdx[jj][2] / Re * detJ * weight;

                        // ADVECTION TERM
                        for(int mm=0; mm<3; mm++){
                            Klocal(IU, JU) += 5e-1 * N[ii] * advel[mm] * dNdx[jj][mm] * detJ * weight;
                            Klocal(IV, JV) += 5e-1 * N[ii] * advel[mm] * dNdx[jj][mm] * detJ * weight;
                            Klocal(IW, JW) += 5e-1 * N[ii] * advel[mm] * dNdx[jj][mm] * detJ * weight;
                        }

                        // PRESSURE TERM
                        Klocal(IU, JP) -= N[jj] * dNdx[ii][0] * detJ * weight;
                        Klocal(IV, JP) -= N[jj] * dNdx[ii][1] * detJ * weight;
                        Klocal(IW, JP) -= N[jj] * dNdx[ii][2] * detJ * weight;

                        // CONTINUITY TERM
                        Klocal(IP, JU) += N[ii] * dNdx[jj][0] * detJ * weight;
                        Klocal(IP, JV) += N[ii] * dNdx[jj][1] * detJ * weight;
                        Klocal(IP, JW) += N[ii] * dNdx[jj][2] * detJ * weight;

                        // DARCY TERM
                        Klocal(IU, JU) += 5e-1 * f * N[ii] * N[jj] * detJ * weight;
                        Klocal(IV, JV) += 5e-1 * f * N[ii] * N[jj] * detJ * weight;
                        Klocal(IW, JW) += 5e-1 * f * N[ii] * N[jj] * detJ * weight;   

                        // SUPG TERM 
                        // MASS TERM 
                        for(int mm=0; mm<3; mm++){
                            Klocal(IU, JU) += tau * dNdx[ii][mm] * advel[mm] * N[jj] / dt * detJ * weight;
                            Klocal(IV, JV) += tau * dNdx[ii][mm] * advel[mm] * N[jj] / dt * detJ * weight;
                            Klocal(IW, JW) += tau * dNdx[ii][mm] * advel[mm] * N[jj] / dt * detJ * weight;
                        }

                        // ADVECTION TERM
                        for(int mm=0; mm<3; mm++){
                            for(int nn=0; nn<3; nn++){
                                Klocal(IU, JU) += 5e-1 * tau * advel[nn] * dNdx[ii][nn] * advel[mm] * dNdx[jj][mm] * detJ * weight;
                                Klocal(IV, JV) += 5e-1 * tau * advel[nn] * dNdx[ii][nn] * advel[mm] * dNdx[jj][mm] * detJ * weight;
                                Klocal(IW, JW) += 5e-1 * tau * advel[nn] * dNdx[ii][nn] * advel[mm] * dNdx[jj][mm] * detJ * weight;
                            }
                        }

                        // PRESSURE TERM
                        for(int mm=0; mm<3; mm++){
                            Klocal(IU, JP) += tau * dNdx[ii][mm] * advel[mm] * dNdx[jj][0] * detJ * weight;
                            Klocal(IV, JP) += tau * dNdx[ii][mm] * advel[mm] * dNdx[jj][1] * detJ * weight;
                            Klocal(IW, JP) += tau * dNdx[ii][mm] * advel[mm] * dNdx[jj][2] * detJ * weight;
                        }

                        // PSPG TERM
                        // MASS TERM 
                        Klocal(IP, JU) += tau * dNdx[ii][0] * N[jj] / dt * detJ * weight;
                        Klocal(IP, JV) += tau * dNdx[ii][1] * N[jj] / dt * detJ * weight;
                        Klocal(IP, JW) += tau * dNdx[ii][2] * N[jj] / dt * detJ * weight;

                        // ADVECTION TERM
                        for(int mm=0; mm<3; mm++){
                            Klocal(IP, JU) += 5e-1 * tau * dNdx[ii][0] * advel[mm] * dNdx[jj][mm] * detJ * weight;
                            Klocal(IP, JV) += 5e-1 * tau * dNdx[ii][1] * advel[mm] * dNdx[jj][mm] * detJ * weight;
                            Klocal(IP, JW) += 5e-1 * tau * dNdx[ii][2] * advel[mm] * dNdx[jj][mm] * detJ * weight;
                        }
                        // PRESSURE TERM
                        Klocal(IP, JP) += tau * K[ii][jj] * detJ * weight;
                    }

                    // MASS TERM
                    Flocal(IU) += N[ii] * vel[0] / dt * detJ * weight;
                    Flocal(IV) += N[ii] * vel[1] / dt * detJ * weight;
                    Flocal(IW) += N[ii] * vel[2] / dt * detJ * weight;

                    // DIFFUSION TERM
                    for(int mm=0; mm<3; mm++){
                        if(mm == 0){n1 = 2e0; n2 = 1e0; n3 = 1e0;}
                        if(mm == 1){n1 = 1e0; n2 = 2e0; n3 = 1e0;}
                        if(mm == 2){n1 = 1e0; n2 = 1e0; n3 = 2e0;}
                        Flocal(IU) -= 5e-1 * n1 * dNdx[ii][mm] * dvdx[0][mm] / Re * detJ * weight;
                        Flocal(IV) -= 5e-1 * n2 * dNdx[ii][mm] * dvdx[1][mm] / Re * detJ * weight;
                        Flocal(IW) -= 5e-1 * n3 * dNdx[ii][mm] * dvdx[2][mm] / Re * detJ * weight;
                    }
                    Flocal(IU) -= 5e-1 * dNdx[ii][1] * dvdx[1][0] / Re * detJ * weight;
                    Flocal(IU) -= 5e-1 * dNdx[ii][2] * dvdx[2][0] / Re * detJ * weight;
                    Flocal(IV) -= 5e-1 * dNdx[ii][0] * dvdx[0][1] / Re * detJ * weight;
                    Flocal(IV) -= 5e-1 * dNdx[ii][2] * dvdx[2][1] / Re * detJ * weight;
                    Flocal(IW) -= 5e-1 * dNdx[ii][0] * dvdx[0][2] / Re * detJ * weight;
                    Flocal(IW) -= 5e-1 * dNdx[ii][1] * dvdx[1][2] / Re * detJ * weight;

                    // ADVECTION TERM
                    for(int mm=0;mm<3;mm++){
                        Flocal(IU) -= 5e-1 * N[ii] * advel[mm] * dvdx[0][mm] * detJ * weight;
                        Flocal(IV) -= 5e-1 * N[ii] * advel[mm] * dvdx[1][mm] * detJ * weight;
                        Flocal(IW) -= 5e-1 * N[ii] * advel[mm] * dvdx[2][mm] * detJ * weight;
                    }
                    // DARCY TERM
                    Flocal(IU) -= 5e-1 * f * N[ii] * vel[0] * detJ * weight;
                    Flocal(IV) -= 5e-1 * f * N[ii] * vel[1] * detJ * weight;
                    Flocal(IW) -= 5e-1 * f * N[ii] * vel[2] * detJ * weight;  

                    // SUPG TERM
                    // MASS TERM
                    for(int mm=0; mm<3; mm++){
                        Flocal(IU) += tau * dNdx[ii][mm] * advel[mm] * vel[0] / dt * detJ * weight;
                        Flocal(IV) += tau * dNdx[ii][mm] * advel[mm] * vel[1] / dt * detJ * weight;
                        Flocal(IW) += tau * dNdx[ii][mm] * advel[mm] * vel[2] / dt * detJ * weight;
                    }

                    // ADVECTION TERM
                    for(int mm=0; mm<3; mm++){
                        for(int nn=0; nn<3; nn++){
                            Flocal(IU) -= 5e-1 * tau * vel[nn] * dNdx[ii][nn] * advel[mm] * dvdx[0][mm] * detJ * weight;
                            Flocal(IV) -= 5e-1 * tau * vel[nn] * dNdx[ii][nn] * advel[mm] * dvdx[1][mm] * detJ * weight;
                            Flocal(IW) -= 5e-1 * tau * vel[nn] * dNdx[ii][nn] * advel[mm] * dvdx[2][mm] * detJ * weight;
                        }
                    }

                    // PSPG TERM
                    // MASS TERM 
                    Flocal(IP) += tau * dNdx[ii][0] * vel[0] / dt * detJ * weight;
                    Flocal(IP) += tau * dNdx[ii][1] * vel[1] / dt * detJ * weight;
                    Flocal(IP) += tau * dNdx[ii][2] * vel[2] / dt * detJ * weight;

                    // ADVECTION TERM
                    for(int mm=0; mm<3; mm++){
                        Flocal(IP) -= 5e-1 * tau * dNdx[ii][0] * advel[mm] * dvdx[0][mm] * detJ * weight;
                        Flocal(IP) -= 5e-1 * tau * dNdx[ii][1] * advel[mm] * dvdx[1][mm] * detJ * weight;
                        Flocal(IP) -= 5e-1 * tau * dNdx[ii][2] * advel[mm] * dvdx[2][mm] * detJ * weight;
                    }
                }
            }
        }
    }
}


void DirectProblem::setVelocityValue(double (&vel)[3], double (&advel)[3], double (&dvdx)[3][3],
                                     std::vector<double> &N, std::vector<std::vector<double>> &dNdx, 
                                     const int ic, const int t)
{
    for(int d=0; d<dim; d++){
        for(int p=0; p<grid.cell.nNodesInCell; p++){
            advel[d] += N[p] * (1.5 * grid.node.v[grid.cell(ic).node[p]][d] - 0.5 * grid.node.vPrev[grid.cell(ic).node[p]][d]);
            vel[d] += N[p] * grid.node.v[grid.cell(ic).node[p]][d];
        }
    }
    for(int d=0; d<dim; d++){
        for(int e=0; e<dim; e++){
            for(int p=0; p<grid.cell.nNodesInCell; p++){
                dvdx[d][e] += dNdx[p][e] * grid.node.v[grid.cell(ic).node[p]][d];
            }
        }
    }
}

