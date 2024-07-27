/**
 * @file AssemblyAdjoint.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "InverseProblem.h"

/*********************************
 * @brief Assemble adjoint system.
 */
void Adjoint::matrixAssemblyAdjoint(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal, 
                                    Function &func, const int ic, const int t)
{   
    Gauss g2(2);
    double dxdr[3][3];
    
    for(int p=0; p<grid.cell.nNodesInCell; p++){
        for(int d=0; d<main.dim; d++){
            func.xCurrent[p][d] = grid.node.x[grid.cell(ic).node[p]][d];
        }
    }
    double he = fabs(func.xCurrent[1][0] - func.xCurrent[0][0]);
    double f = main.resistance * main.alpha * (1e0 - grid.cell(ic).phi) 
                                     / (main.alpha + grid.cell(ic).phi);

    for(int i1=0; i1<2; i1++){
        for(int i2=0; i2<2; i2++){
            for(int i3=0; i3<2; i3++){
                ShapeFunction3D::C3D8_N(func.N, g2.point[i1], g2.point[i2], g2.point[i3]);
                ShapeFunction3D::C3D8_dNdr(func.dNdr, g2.point[i1], g2.point[i2], g2.point[i3]);

                MathFEM::comp_dxdr(dxdr, func.dNdr, func.xCurrent, grid.cell.nNodesInCell);
                MathFEM::comp_dNdx(func.dNdx, func.dNdr, dxdr, grid.cell.nNodesInCell);

                func.detJ = MathCommon::compDeterminant_3x3(dxdr);
                func.weight = g2.weight[i1] * g2.weight[i2] * g2.weight[i3];

                setValue(main, func, ic, t);
                tau = MathFEM::comp_tau(wk1, he, main.Re, main.dt);

                for(int ii=0; ii<grid.cell.nNodesInCell; ii++){  
                    updateRowIndex(ii, ic);
                    for(int jj=0; jj<grid.cell.nNodesInCell; jj++){  
                        updateColumnIndex(jj, ic);
                        adjointGaussIntegralLHS(main, Klocal, func, f, ii, jj);
                    }
                    adjointGaussIntegralRHS(main, Flocal, func, f, ii);                    
                }
            }
        }
    }

}

/**********************************************************************
 * @brief Compute element sfiffness LHS matrix on gauss integral point. 
 */
void Adjoint::adjointGaussIntegralLHS(DirectProblem &main, MatrixXd &Klocal, Function &func, 
                                      const double f, const int ii, const int jj)
{
    int n1, n2, n3;
    func.vol = func.detJ * func.weight;

    func.K[ii][jj] = 0e0;
    for(int k=0;k<3;k++){
        func.K[ii][jj] += func.dNdx[ii][k] * func.dNdx[jj][k];
    }

    // Mass term
    Klocal(IU, JU) += func.N[ii] * func.N[jj] / main.dt * func.vol;
    Klocal(IV, JV) += func.N[ii] * func.N[jj] / main.dt * func.vol;
    Klocal(IW, JW) += func.N[ii] * func.N[jj] / main.dt * func.vol;

     // Diffusion term
    for(int mm=0; mm<3; mm++){
        if(mm == 0){n1 = 2e0; n2 = 1e0; n3 = 1e0;}
        if(mm == 1){n1 = 1e0; n2 = 2e0; n3 = 1e0;}
        if(mm == 2){n1 = 1e0; n2 = 1e0; n3 = 2e0;}
        Klocal(IU, JU) += 5e-1 * n1 * func.dNdx[ii][mm] * func.dNdx[jj][mm] / main.Re * func.vol;
        Klocal(IV, JV) += 5e-1 * n2 * func.dNdx[ii][mm] * func.dNdx[jj][mm] / main.Re * func.vol;
        Klocal(IW, JW) += 5e-1 * n3 * func.dNdx[ii][mm] * func.dNdx[jj][mm] / main.Re * func.vol;
    }
    Klocal(IU, JV) += 5e-1 * func.dNdx[ii][1] * func.dNdx[jj][0] / main.Re * func.vol;
    Klocal(IU, JW) += 5e-1 * func.dNdx[ii][2] * func.dNdx[jj][0] / main.Re * func.vol;
    Klocal(IV, JU) += 5e-1 * func.dNdx[ii][0] * func.dNdx[jj][1] / main.Re * func.vol;
    Klocal(IV, JW) += 5e-1 * func.dNdx[ii][2] * func.dNdx[jj][1] / main.Re * func.vol;
    Klocal(IW, JU) += 5e-1 * func.dNdx[ii][0] * func.dNdx[jj][2] / main.Re * func.vol;
    Klocal(IW, JV) += 5e-1 * func.dNdx[ii][1] * func.dNdx[jj][2] / main.Re * func.vol;

    // Advection term
    for(int d=0; d<main.dim; d++){
        Klocal(IU, JU) += 5e-1 * func.dNdx[ii][d] * advk1[d] * func.N[jj] * func.vol;
        Klocal(IV, JV) += 5e-1 * func.dNdx[ii][d] * advk1[d] * func.N[jj] * func.vol;
        Klocal(IW, JW) += 5e-1 * func.dNdx[ii][d] * advk1[d] * func.N[jj] * func.vol;
    }

    // Pressure term
    Klocal(IU, JP) -= func.dNdx[ii][0] * func.N[jj] * func.vol;
    Klocal(IV, JP) -= func.dNdx[ii][1] * func.N[jj] * func.vol;
    Klocal(IW, JP) -= func.dNdx[ii][2] * func.N[jj] * func.vol;  
    
    // Continuity term
    Klocal(IP, JU) += func.N[ii] * func.dNdx[jj][0] * func.vol;
    Klocal(IP, JV) += func.N[ii] * func.dNdx[jj][1] * func.vol;
    Klocal(IP, JW) += func.N[ii] * func.dNdx[jj][2] * func.vol;
    
    // Darcy term
    Klocal(IU, JU) += 5e-1 * f * func.N[ii] * func.N[jj] * func.vol;
    Klocal(IV, JV) += 5e-1 * f * func.N[ii] * func.N[jj] * func.vol;
    Klocal(IW, JW) += 5e-1 * f * func.N[ii] * func.N[jj] * func.vol;   

    /*
    // PSPG advection term
    for(int d=0; d<3; d++){
        Klocal(IP, JU) += 5e-1 * tau * dNdx[ii][0] * vel[d] * dNdx[jj][d] * func.vol;
        Klocal(IP, JV) += 5e-1 * tau * dNdx[ii][1] * vel[d] * dNdx[jj][d] * func.vol;
        Klocal(IP, JW) += 5e-1 * tau * dNdx[ii][2] * vel[d] * dNdx[jj][d] * func.vol;
    }
    */

    // PSPG pressure term
    Klocal(IP, JP) += tau * func.K[ii][jj] * func.vol;
}

/**********************************************************************
 * @brief Compute element sfiffness RHS matrix on gauss integral point. 
 */
void Adjoint::adjointGaussIntegralRHS(DirectProblem &main, VectorXd &Flocal, 
                                      Function &func, const double f, const int ii)
{
    int n1, n2, n3;
    func.vol = func.detJ * func.weight;

    // Mass term
    Flocal(IU) += func.N[ii] * wk1[0] / main.dt * func.vol;
    Flocal(IV) += func.N[ii] * wk1[1] / main.dt * func.vol;
    Flocal(IW) += func.N[ii] * wk1[2] / main.dt * func.vol;
                
    // Diffusion term
    for(int d=0; d<3; d++){
        if(d == 0){n1 = 2e0; n2 = 1e0; n3 = 1e0;}
        if(d == 1){n1 = 1e0; n2 = 2e0; n3 = 1e0;}
        if(d == 2){n1 = 1e0; n2 = 1e0; n3 = 2e0;}
        Flocal(IU) -= 5e-1 * n1 * func.dNdx[ii][d] * dwk1dx[0][d] / main.Re * func.vol;
        Flocal(IV) -= 5e-1 * n2 * func.dNdx[ii][d] * dwk1dx[1][d] / main.Re * func.vol;
        Flocal(IW) -= 5e-1 * n3 * func.dNdx[ii][d] * dwk1dx[2][d] / main.Re * func.vol;
    }
    Flocal(IU) -= 5e-1 * func.dNdx[ii][1] * dwk1dx[1][0] / main.Re * func.vol;
    Flocal(IU) -= 5e-1 * func.dNdx[ii][2] * dwk1dx[2][0] / main.Re * func.vol;
    Flocal(IV) -= 5e-1 * func.dNdx[ii][0] * dwk1dx[0][1] / main.Re * func.vol;
    Flocal(IV) -= 5e-1 * func.dNdx[ii][2] * dwk1dx[2][1] / main.Re * func.vol;
    Flocal(IW) -= 5e-1 * func.dNdx[ii][0] * dwk1dx[0][2] / main.Re * func.vol;
    Flocal(IW) -= 5e-1 * func.dNdx[ii][1] * dwk1dx[1][2] / main.Re * func.vol;
                
    // Advection term
    Flocal(IU) -= 0.5 * 1.5 * func.N[ii] * dvk1dx[0][0] * wk1[0] * func.vol;
    Flocal(IU) -= 0.5 * 1.5 * func.N[ii] * dvkdx[0][0] * wk1[0] * func.vol;
    for(int d=0; d<main.dim; d++){
        Flocal(IU) -= 0.5 * func.dNdx[ii][d] * advk2[d] * wk1[0] * func.vol;
    }
    Flocal(IU) -= 0.5 * 1.5 * func.N[ii] * dvk1dx[1][0] * wk1[1] * func.vol;
    Flocal(IU) -= 0.5 * 1.5 * func.N[ii] * dvkdx[1][0] * wk1[1] * func.vol;
    Flocal(IU) -= 0.5 * 1.5 * func.N[ii] * dvk1dx[2][0] * wk1[2] * func.vol;
    Flocal(IU) -= 0.5 * 1.5 * func.N[ii] * dvkdx[2][0] * wk1[2] * func.vol;
    for(int d=0; d<main.dim; d++){
        Flocal(IU) -= 0.5 * (-0.5) * func.N[ii] * dvk2dx[d][0] * wk2[d] * func.vol;
        Flocal(IU) -= 0.5 * (-0.5) * func.N[ii] * dvk1dx[d][0] * wk2[d] * func.vol;
    }
                
    Flocal(IV) -= 0.5 * 1.5 * func.N[ii] * dvk1dx[1][1] * wk1[1] * func.vol;
    Flocal(IV) -= 0.5 * 1.5 * func.N[ii] * dvkdx[1][1] * wk1[1] * func.vol;
    for(int d=0; d<main.dim; d++){
        Flocal(IV) -= 0.5 * func.dNdx[ii][d] * advk2[d] * wk1[1] * func.vol;
    }
    Flocal(IV) -= 0.5 * 1.5 * func.N[ii] * dvk1dx[0][1] * wk1[0] * func.vol;
    Flocal(IV) -= 0.5 * 1.5 * func.N[ii] * dvkdx[0][1] * wk1[0] * func.vol;
    Flocal(IV) -= 0.5 * 1.5 * func.N[ii] * dvk1dx[2][1] * wk1[2] * func.vol;
    Flocal(IV) -= 0.5 * 1.5 * func.N[ii] * dvkdx[2][1] * wk1[2] * func.vol;
    for(int d=0; d<main.dim; d++){
        Flocal(IV) -= 0.5 * (-0.5) * func.N[ii] * dvk2dx[d][1] * wk2[d] * func.vol;
        Flocal(IV) -= 0.5 * (-0.5) * func.N[ii] * dvk1dx[d][1] * wk2[d] * func.vol;
    }
                
    Flocal(IW) -= 0.5 * 1.5 * func.N[ii] * dvk1dx[2][2] * wk1[2] * func.vol;
    Flocal(IW) -= 0.5 * 1.5 * func.N[ii] * dvkdx[2][2] * wk1[2] * func.vol;
    for(int d=0; d<main.dim; d++){
        Flocal(IW) -= 0.5 * func.dNdx[ii][d] * advk2[d] * wk1[2] * func.vol;
    }
    Flocal(IW) -= 0.5 * 1.5 * func.N[ii] * dvk1dx[0][2] * wk1[0] * func.vol;
    Flocal(IW) -= 0.5 * 1.5 * func.N[ii] * dvkdx[0][2] * wk1[0] * func.vol;
    Flocal(IW) -= 0.5 * 1.5 * func.N[ii] * dvk1dx[1][2] * wk1[1] * func.vol;
    Flocal(IW) -= 0.5 * 1.5 * func.N[ii] * dvkdx[1][2] * wk1[1] * func.vol;   
    for(int d=0; d<main.dim; d++){
        Flocal(IW) -= 0.5 * (-0.5) * func.N[ii] * dvk2dx[d][2] * wk2[d] * func.vol;
        Flocal(IW) -= 0.5 * (-0.5) * func.N[ii] * dvk1dx[d][2] * wk2[d] * func.vol;
    }
    
    // Darcy term
    Flocal(IU) -= 5e-1 * f * func.N[ii] * wk1[0] * func.vol;
    Flocal(IV) -= 5e-1 * f * func.N[ii] * wk1[1] * func.vol;
    Flocal(IW) -= 5e-1 * f * func.N[ii] * wk1[2] * func.vol;
                
    /*
    // PSPG Advection term
    for(int mm=0; mm<3; mm++){
        Flocal(IP) -= 5e-1 * tau * dNdx[ii][0] * vel[mm] * dwgpdx[0][mm] * func.vol;
        Flocal(IP) -= 5e-1 * tau * dNdx[ii][1] * vel[mm] * dwgpdx[1][mm] * func.vol;
        Flocal(IP) -= 5e-1 * tau * dNdx[ii][2] * vel[mm] * dwgpdx[2][mm] * func.vol;
    }
    */
}

/***********************************************
 * @brief Set values needed for matrix assembly 
 *        on gauss integral points.
 */
void Adjoint::setValue(DirectProblem &main, Function &func, const int ic, const int t)
{
    // main var - v
    for(int d=0; d<main.dim; d++){
        vk[d] = 0e0;
        vk1[d] = 0e0;
        vk2[d] = 0e0;
        for(int p=0; p<grid.cell.nNodesInCell; p++){
            int n = grid.cell(ic).node[p];
            if(t == timeMax-1){
                vk[d]  += func.N[p] * main.grid.node.vt[t][n][d];
                vk1[d] = 0e0;
                vk2[d] = 0e0;
            }else if(t == timeMax-2){
                vk[d]  += func.N[p] * main.grid.node.vt[t][n][d];
                vk1[d] += func.N[p] * main.grid.node.vt[t+1][n][d];
                vk2[d] = 0e0;
            }else{
                vk[d]  += func.N[p] * main.grid.node.vt[t][n][d];
                vk1[d] += func.N[p] * main.grid.node.vt[t+1][n][d];
                vk2[d] += func.N[p] * main.grid.node.vt[t+2][n][d];
            }
        }
    }

    // main var - dvdx
    for(int d=0; d<main.dim; d++){
        for(int e=0; e<main.dim; e++){
            dvkdx[d][e] = 0e0;
            dvk1dx[d][e] = 0e0;
            dvk2dx[d][e] = 0e0;
            for(int p=0; p<grid.cell.nNodesInCell; p++){
                int n = grid.cell(ic).node[p];
                if(t == timeMax-1){
                    dvkdx[d][e]  += func.dNdx[p][e] * main.grid.node.vt[t][n][d];
                    dvk1dx[d][e] = 0e0;
                    dvk2dx[d][e] = 0e0;
                }else if(t == timeMax-2){
                    dvkdx[d][e]  += func.dNdx[p][e] * main.grid.node.vt[t][n][d];
                    dvk1dx[d][e] += func.dNdx[p][e] * main.grid.node.vt[t+1][n][d];
                    dvk2dx[d][e] = 0e0;
                }else{
                    dvkdx[d][e]  += func.dNdx[p][e] * main.grid.node.vt[t][n][d];
                    dvk1dx[d][e] += func.dNdx[p][e] * main.grid.node.vt[t+1][n][d];
                    dvk2dx[d][e] += func.dNdx[p][e] * main.grid.node.vt[t+2][n][d];
                }
            }
        }
    }

    // main var - adv
    for(int d=0; d<main.dim; d++){
        advk1[d] = 0e0; 
        advk2[d] = 0e0;
        for(int p=0; p<grid.cell.nNodesInCell; p++){
            int n = grid.cell(ic).node[p];
            if(t == 0){
                advk1[d] += func.N[p] * main.grid.node.v0[n][d];
                advk2[d] += func. N[p] * (1.5 * main.grid.node.vt[t][n][d]
                          - 0.5 * main.grid.node.v0[n][d]);
            }else if(t == 1){
                advk1[d] += func.N[p] * (1.5 * main.grid.node.vt[t-1][n][d]
                          - 0.5 * main.grid.node.v0[n][d]);
                advk2[d] += func.N[p] * (1.5 * main.grid.node.vt[t][n][d] 
                          - 0.5 * main.grid.node.vt[t-1][n][d]);
            }else{
                advk1[d] += func.N[p] * (1.5 * main.grid.node.vt[t-1][n][d] 
                                 - 0.5 * main.grid.node.vt[t-2][n][d]);
                advk2[d] += func.N[p] * (1.5 * main.grid.node.vt[t][n][d] 
                                 - 0.5 * main.grid.node.vt[t-1][n][d]);
            }
        }
    }

    // lagrange multiplier - w
    for(int d=0; d<main.dim; d++){
        wk1[d] = 0e0;
        wk2[d] = 0e0;
        for(int p=0; p<grid.cell.nNodesInCell; p++){
            int n = grid.cell(ic).node[p];
            wk1[d] += func.N[p] * grid.node.w[n][d];
            wk2[d] += func.N[p] * grid.node.wPrev[n][d];
        }
    }

    // lagrange multiplier - dwdx
    for(int d=0; d<main.dim; d++){
        for(int e=0; e<main.dim; e++){
            dwk1dx[d][e] = 0e0;
            dwk2dx[d][e] = 0e0;
            for(int p=0; p<grid.cell.nNodesInCell; p++){
                int n = grid.cell(ic).node[p];
                dwk1dx[d][e] += func.dNdx[p][e] * grid.node.w[n][d];
                dwk2dx[d][e] += func.dNdx[p][e] * grid.node.wPrev[n][d];
            }
        }
    }

}

/***********************************
 * @brief Compute boundary integral.
 */
void Adjoint::boundaryIntegral(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal,
                               Function &func, const int ic, const int ib)
{   
    int nc = grid.dirichlet.nControlNodesInCell;
    
    Gauss g2(2);
    double dxdr[2][2];
    
    for(int p=0; p<nc; p++){
        int n = grid.dirichlet.controlNodeInCell[ib][p];
        for(int d=0; d<main.dim-1; d++){
            func.xCurrent[p][d] = 0e0;
            func.xCurrent[p][d] = main.grid.node.x[n][planeDir[d]];
        }
    }
    
    for(int i1=0; i1<2; i1++){
        for(int i2=0; i2<2; i2++){
            ShapeFunction2D::C2D4_N(func.N, g2.point[i1], g2.point[i2]);
            ShapeFunction2D::C2D4_dNdr(func.dNdr, g2.point[i1], g2.point[i2]);
            MathFEM::comp_dxdr2D(dxdr, func.dNdr, func.xCurrent, nc);
            func.detJ = MathCommon::compDeterminant_2x2(dxdr);
            func.weight = g2.weight[i1] * g2.weight[i2];
            for(int ii=0; ii<nc; ii++){
                updateRowIndexPlane(ii, ic);
                for(int jj=0; jj<nc; jj++){
                    updateColumnIndexPlane(jj, ic);
                    boundaryInGaussIntegral(Klocal, func, ii, jj);
                }
            }
            
        }
    }
}

/*******************************************************************
 * @brief Compute LHS element sfiffness matrix for boundary integral
 *        on gauss integral point. 
 */
void Adjoint::boundaryInGaussIntegral(MatrixXd &Klocal, Function &func, const int ii, const int jj)
{
    func.vol = func.detJ * func.weight;
    
    Klocal(IU, JLU) -= func.N[ii] * func.N[jj] * func.vol;
    Klocal(IV, JLV) -= func.N[ii] * func.N[jj] * func.vol;
    Klocal(IW, JLW) -= func.N[ii] * func.N[jj] * func.vol;
    Klocal(ILU, JU) -= func.N[ii] * func.N[jj] * func.vol; 
    Klocal(ILV, JV) -= func.N[ii] * func.N[jj] * func.vol; 
    Klocal(ILW, JW) -= func.N[ii] * func.N[jj] * func.vol; 
}


/* Advection term for Optimize Then Discretize
 * 
 * for(int d=0; d<main.dim; d++){
 *     Klocal(IU, JU) += 5e-1 * dNdx[ii][d] * vel[d] * N[jj] * func.vol;
 *     Klocal(IV, JV) += 5e-1 * dNdx[ii][d] * vel[d] * N[jj] * func.vol;
 *     Klocal(IW, JW) += 5e-1 * dNdx[ii][d] * vel[d] * N[jj] * func.vol;
 *  }
 *  Klocal(IU, JV) += 5e-1 * N[ii] * dvdx[1][0] * N[jj] * func.vol;
 *  Klocal(IU, JW) += 5e-1 * N[ii] * dvdx[2][0] * N[jj] * func.vol;
 *  Klocal(IV, JU) += 5e-1 * N[ii] * dvdx[0][1] * N[jj] * func.vol;
 *  Klocal(IV, JW) += 5e-1 * N[ii] * dvdx[2][1] * N[jj] * func.vol;
 *  Klocal(IW, JU) += 5e-1 * N[ii] * dvdx[0][2] * N[jj] * func.vol;
 *  Klocal(IW, JV) += 5e-1 * N[ii] * dvdx[1][2] * N[jj] * func.vol;
 */