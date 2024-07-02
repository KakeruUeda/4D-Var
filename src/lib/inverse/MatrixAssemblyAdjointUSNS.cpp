#include "InverseProblem.h"

void Adjoint::matrixAssemblyAdjointUSNS(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal, 
                                        const int ic, const int t)
{
    int n1, n2, n3;
    int IU, IV, IW, IP;
    int JU, JV, JW, JP;
    int ILU, ILV, ILW;
    int nGaussPoint = 2;
    
    std::vector<std::vector<double>> xCurrent;
    std::vector<double> N;
    std::vector<std::vector<double>> dNdr;
    std::vector<std::vector<double>> dNdx;
    std::vector<std::vector<double>> K;
    
    VecTool::resize(xCurrent, grid.cell.nNodesInCell, main.dim);
    VecTool::resize(N, grid.cell.nNodesInCell);
    VecTool::resize(dNdr, grid.cell.nNodesInCell, main.dim);
    VecTool::resize(dNdx, grid.cell.nNodesInCell, main.dim);
    VecTool::resize(K, grid.cell.nNodesInCell, grid.cell.nNodesInCell);
    
    for(int p=0; p<grid.cell.nNodesInCell; p++){
        for(int d=0; d<main.dim; d++){
            xCurrent[p][d] = grid.node.x[grid.cell(ic).node[p]][d];
        }
    }
            
    int s[grid.cell.nNodesInCell];
    for(int p=0; p<grid.cell.nNodesInCell; p++){
        s[p] = 0;
    }       
    for(int p=1; p<grid.cell.nNodesInCell; p++){
        s[p] = s[p-1] + grid.node.nDofsOnNode[grid.cell(ic).node[p-1]]; 
    }

    double he = fabs(xCurrent[1][0] - xCurrent[0][0]);
    double f = main.resistance * main.alpha * (1e0 - grid.cell(ic).phi) / (main.alpha + grid.cell(ic).phi);
    double detJ, weight;
    Gauss gauss(nGaussPoint);

    double dxdr[3][3] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    for(int i1=0; i1<nGaussPoint; i1++){
        for(int i2=0; i2<nGaussPoint; i2++){
            for(int i3=0; i3<nGaussPoint; i3++){
                ShapeFunction3D::C3D8_N(N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                ShapeFunction3D::C3D8_dNdr(dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                MathFEM::calc_dxdr(dxdr, dNdr, xCurrent, grid.cell.nNodesInCell);
                MathFEM::calc_dNdx(dNdx, dNdr, dxdr, grid.cell.nNodesInCell);
                detJ = MathCommon::calcDeterminant_3x3(dxdr);
                weight = gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];

                double vel[3];
                double dvdx[3][3];
                double wgp[3];
                double dwgpdx[3][3];

                for(int d=0; d<main.dim; d++){
                    vel[d] = 0e0;
                    for(int p=0; p<grid.cell.nNodesInCell; p++){
                        vel[d] += N[p] * main.grid.node.vt[t][grid.cell(ic).node[p]][d];
                    }
                }

                for(int d=0; d<main.dim; d++){
                    for(int e=0; e<main.dim; e++){
                        dvdx[d][e] = 0e0;
                        for(int p=0; p<grid.cell.nNodesInCell; p++){
                            dvdx[d][e] += dNdx[p][e] * main.grid.node.vt[t][grid.cell(ic).node[p]][d];
                        }
                    }
                }    
                for(int d=0; d<main.dim; d++){
                    wgp[d] = 0e0;
                    for(int p=0; p<grid.cell.nNodesInCell; p++){
                        wgp[d] += N[p] * grid.node.w[grid.cell(ic).node[p]][d];
                    }
                }  
                for(int d=0; d<main.dim; d++){
                    for(int e=0; e<main.dim; e++){
                        dwgpdx[d][e] = 0e0;
                        for(int p=0; p<grid.cell.nNodesInCell; p++){
                            dwgpdx[d][e] += dNdx[p][e] * grid.node.w[grid.cell(ic).node[p]][d];
                        }
                    }
                }    
                
                vector<double> tmp(grid.cell.nNodesInCell, 0e0);
                for(int p=0; p<grid.cell.nNodesInCell; p++){
                    for(int d=0; d<main.dim; d++){
                        tmp[p] += dNdx[p][d] * vel[d];
                    } 
                }

                //setVelocityValue(vel, advel, dvdx, N, dNdx, ic, t);
                double tau = MathFEM::calc_tau(wgp, he, main.Re, main.dt);

                for(int ii=0; ii<grid.cell.nNodesInCell; ii++){  
                    IU = s[ii];
                    IV = IU + 1;
                    IW = IU + 2;
                    IP = IU + 3;
                    for(int jj=0; jj<grid.cell.nNodesInCell; jj++){  
                        JU = s[jj];
                        JV = JU + 1;
                        JW = JU + 2;
                        JP = JU + 3;

                        K[ii][jj] = 0e0;
                        for(int k=0;k<3;k++){
                            K[ii][jj] += dNdx[ii][k] * dNdx[jj][k];
                        }

                        // MASS TERM
                        Klocal(IU, JU) += N[ii] * N[jj] / main.dt * detJ * weight;
                        Klocal(IV, JV) += N[ii] * N[jj] / main.dt * detJ * weight;
                        Klocal(IW, JW) += N[ii] * N[jj] / main.dt * detJ * weight;

                         // DIFFUSION TERM
                        for(int mm=0; mm<3; mm++){
                            if(mm == 0){n1 = 2e0; n2 = 1e0; n3 = 1e0;}
                            if(mm == 1){n1 = 1e0; n2 = 2e0; n3 = 1e0;}
                            if(mm == 2){n1 = 1e0; n2 = 1e0; n3 = 2e0;}
                            Klocal(IU, JU) += 5e-1 * n1 * dNdx[ii][mm] * dNdx[jj][mm] / main.Re * detJ * weight;
                            Klocal(IV, JV) += 5e-1 * n2 * dNdx[ii][mm] * dNdx[jj][mm] / main.Re * detJ * weight;
                            Klocal(IW, JW) += 5e-1 * n3 * dNdx[ii][mm] * dNdx[jj][mm] / main.Re * detJ * weight;
                        }
                        Klocal(IU, JV) += 5e-1 * dNdx[ii][1] * dNdx[jj][0] / main.Re * detJ * weight;
                        Klocal(IU, JW) += 5e-1 * dNdx[ii][2] * dNdx[jj][0] / main.Re * detJ * weight;
                        Klocal(IV, JU) += 5e-1 * dNdx[ii][0] * dNdx[jj][1] / main.Re * detJ * weight;
                        Klocal(IV, JW) += 5e-1 * dNdx[ii][2] * dNdx[jj][1] / main.Re * detJ * weight;
                        Klocal(IW, JU) += 5e-1 * dNdx[ii][0] * dNdx[jj][2] / main.Re * detJ * weight;
                        Klocal(IW, JV) += 5e-1 * dNdx[ii][1] * dNdx[jj][2] / main.Re * detJ * weight;
                        
                        // ADVECTION TERM 
                        Klocal(IU, JU) += 5e-1 * N[ii] * N[jj] * dvdx[0][0] * detJ * weight;
                        Klocal(IU, JV) += 5e-1 * N[ii] * N[jj] * dvdx[0][1] * detJ * weight;
                        Klocal(IU, JW) += 5e-1 * N[ii] * N[jj] * dvdx[0][2] * detJ * weight;
                        Klocal(IV, JU) += 5e-1 * N[ii] * N[jj] * dvdx[1][0] * detJ * weight;
                        Klocal(IV, JV) += 5e-1 * N[ii] * N[jj] * dvdx[1][1] * detJ * weight;
                        Klocal(IV, JW) += 5e-1 * N[ii] * N[jj] * dvdx[1][2] * detJ * weight;
                        Klocal(IW, JU) += 5e-1 * N[ii] * N[jj] * dvdx[2][0] * detJ * weight;
                        Klocal(IW, JV) += 5e-1 * N[ii] * N[jj] * dvdx[2][1] * detJ * weight;
                        Klocal(IW, JW) += 5e-1 * N[ii] * N[jj] * dvdx[2][2] * detJ * weight;
    
                        Klocal(IU, JU) += 5e-1 * tmp[ii] * N[jj] * detJ * weight;
                        Klocal(IV, JV) += 5e-1 * tmp[ii] * N[jj] * detJ * weight;
                        Klocal(IW, JW) += 5e-1 * tmp[ii] * N[jj] * detJ * weight;
        
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
                            Klocal(IU, JU) += tau * dNdx[ii][mm] * vel[mm] * N[jj] / main.dt * detJ * weight;
                            Klocal(IV, JV) += tau * dNdx[ii][mm] * vel[mm] * N[jj] / main.dt * detJ * weight;
                            Klocal(IW, JW) += tau * dNdx[ii][mm] * vel[mm] * N[jj] / main.dt * detJ * weight;
                        }
        
                        // ADVECTION TERM
                        for(int mm=0; mm<3; mm++){
                            for(int nn=0; nn<3; nn++){
                                Klocal(IU, JU) += 5e-1 * tau * vel[nn] * dNdx[ii][nn] * vel[mm] * dNdx[jj][mm] * detJ * weight;
                                Klocal(IV, JV) += 5e-1 * tau * vel[nn] * dNdx[ii][nn] * vel[mm] * dNdx[jj][mm] * detJ * weight;
                                Klocal(IW, JW) += 5e-1 * tau * vel[nn] * dNdx[ii][nn] * vel[mm] * dNdx[jj][mm] * detJ * weight;
                            }
                        }
       
                        // PRESSURE TERM
                        for(int mm=0; mm<3; mm++){
                            Klocal(IU, JP) += tau * dNdx[ii][mm] * vel[mm] * dNdx[jj][0] * detJ * weight;
                            Klocal(IV, JP) += tau * dNdx[ii][mm] * vel[mm] * dNdx[jj][1] * detJ * weight;
                            Klocal(IW, JP) += tau * dNdx[ii][mm] * vel[mm] * dNdx[jj][2] * detJ * weight;
                        }

                        // PSPG TERM
                        // PRESSURE TERM
                        Klocal(IP, JP) += tau * K[ii][jj] * detJ * weight;
                    }

                    // MASS TERM
                    Flocal(IU) += N[ii] * wgp[0] / main.dt * detJ * weight;
                    Flocal(IV) += N[ii] * wgp[1] / main.dt * detJ * weight;
                    Flocal(IW) += N[ii] * wgp[2] / main.dt * detJ * weight;

                    // DIFFUSION TERM
                    for(int d=0; d<3; d++){
                        if(d == 0){n1 = 2e0; n2 = 1e0; n3 = 1e0;}
                        if(d == 1){n1 = 1e0; n2 = 2e0; n3 = 1e0;}
                        if(d == 2){n1 = 1e0; n2 = 1e0; n3 = 2e0;}
                        Flocal(IU) -= 5e-1 * n1 * dNdx[ii][d] * dwgpdx[0][d] / main.Re * detJ * weight;
                        Flocal(IV) -= 5e-1 * n2 * dNdx[ii][d] * dwgpdx[1][d] / main.Re * detJ * weight;
                        Flocal(IW) -= 5e-1 * n3 * dNdx[ii][d] * dwgpdx[2][d] / main.Re * detJ * weight;
                    }
                    Flocal(IU) -= 5e-1 * dNdx[ii][1] * dwgpdx[1][0] / main.Re * detJ * weight;
                    Flocal(IU) -= 5e-1 * dNdx[ii][2] * dwgpdx[2][0] / main.Re * detJ * weight;
                    Flocal(IV) -= 5e-1 * dNdx[ii][0] * dwgpdx[0][1] / main.Re * detJ * weight;
                    Flocal(IV) -= 5e-1 * dNdx[ii][2] * dwgpdx[2][1] / main.Re * detJ * weight;
                    Flocal(IW) -= 5e-1 * dNdx[ii][0] * dwgpdx[0][2] / main.Re * detJ * weight;
                    Flocal(IW) -= 5e-1 * dNdx[ii][1] * dwgpdx[1][2] / main.Re * detJ * weight;

                    // ADVECTION TERM
                    for(int d=0; d<3; d++){
                        Flocal(IU) -= 5e-1 * N[ii] * wgp[d] * dvdx[0][d] * detJ * weight;
                        Flocal(IV) -= 5e-1 * N[ii] * wgp[d] * dvdx[1][d] * detJ * weight;
                        Flocal(IW) -= 5e-1 * N[ii] * wgp[d] * dvdx[2][d] * detJ * weight;
                    }
                    Klocal(IU, JU) -= 5e-1 * tmp[ii] * wgp[0] * detJ * weight;
                    Klocal(IV, JV) -= 5e-1 * tmp[ii] * wgp[1] * detJ * weight;
                    Klocal(IW, JW) -= 5e-1 * tmp[ii] * wgp[2] * detJ * weight;

                    // DARCY TERM
                    Flocal(IU) -= 5e-1 * f * N[ii] * wgp[0] * detJ * weight;
                    Flocal(IV) -= 5e-1 * f * N[ii] * wgp[1] * detJ * weight;
                    Flocal(IW) -= 5e-1 * f * N[ii] * wgp[2] * detJ * weight;  

                    // SUPG TERM
                    // MASS TERM
                    for(int mm=0; mm<3; mm++){
                        Flocal(IU) += tau * dNdx[ii][mm] * vel[mm] * wgp[0] / main.dt * detJ * weight;
                        Flocal(IV) += tau * dNdx[ii][mm] * vel[mm] * wgp[1] / main.dt * detJ * weight;
                        Flocal(IW) += tau * dNdx[ii][mm] * vel[mm] * wgp[2] / main.dt * detJ * weight;
                    }
         
                    // ADVECTION TERM
                    for(int mm=0; mm<3; mm++){
                        for(int nn=0; nn<3; nn++){
                            Flocal(IU) -= 5e-1 * tau * wgp[nn] * dNdx[ii][nn] * vel[mm] * dvdx[0][mm] * detJ * weight;
                            Flocal(IV) -= 5e-1 * tau * wgp[nn] * dNdx[ii][nn] * vel[mm] * dvdx[1][mm] * detJ * weight;
                            Flocal(IW) -= 5e-1 * tau * wgp[nn] * dNdx[ii][nn] * vel[mm] * dvdx[2][mm] * detJ * weight;
                        }
                    }
                }
            }
        }
    }

}

