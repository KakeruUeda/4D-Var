/**
 * @file AssemblyUSNS.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "DirectProblem.h"

/**************************
 * @brief Assemble system.
 */
void DirectProblem::matrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, Function &func, const int ic, const int t)
{
    for(int p=0; p<grid.cell.nNodesInCell; p++){
       for(int d=0; d<dim; d++){
            func.xCurrent[p][d] = 0e0;
            func.xCurrent[p][d] = grid.node.x[grid.cell(ic).node[p]][d];
       } 
    }       
             
    double he = fabs(func.xCurrent[1][0] - func.xCurrent[0][0]);
    double f = resistance * alpha * (1e0 - grid.cell(ic).phi) 
                                / (alpha + grid.cell(ic).phi);
                                
    Gauss g2(2);
    for(int i1=0; i1<2; i1++){
        for(int i2=0; i2<2; i2++){
            for(int i3=0; i3<2; i3++){
                setValuesInGaussIntegral(func, g2, he, i1, i2, i3, ic, t);
                for(int ii=0; ii<grid.cell.nNodesInCell; ii++){  
                    updateRowIndex(ii, ic);
                    for(int jj=0; jj<grid.cell.nNodesInCell; jj++){  
                        updateColumnIndex(jj, ic);
                        usnsGaussIntegralLHS(Klocal, func, f, ii, jj);
                    }
                    usnsGaussIntegralRHS(Flocal, func, f, ii);    
                }
            }
        }
    }

}

/**************************************
 * @brief Set Values in gauss integral.
 */
void DirectProblem::setValuesInGaussIntegral(Function &func, Gauss &g2, const double he,
                                            const int i1, const int i2, const int i3, 
                                            const int ic, const int t)
{
    double dxdr[3][3];
    ShapeFunction3D::C3D8_N(func.N, g2.point[i1], g2.point[i2], g2.point[i3]);
    ShapeFunction3D::C3D8_dNdr(func.dNdr, g2.point[i1], g2.point[i2], g2.point[i3]);

    MathFEM::comp_dxdr(dxdr, func.dNdr, func.xCurrent, grid.cell.nNodesInCell);
    MathFEM::comp_dNdx(func.dNdx, func.dNdr, dxdr, grid.cell.nNodesInCell);
                
    func.detJ = MathCommon::compDeterminant_3x3(dxdr);
    func.weight = g2.weight[i1] * g2.weight[i2] * g2.weight[i3];

    setVelocityValue(func, ic, t);
    tau = MathFEM::comp_tau(advgp, he, Re, dt);
}

/**********************************************************************
 * @brief Compute element sfiffness LHS matrix on gauss integral point. 
 */
void DirectProblem::usnsGaussIntegralLHS(MatrixXd &Klocal, Function &func, const double f, const int ii, const int jj)
{
    int n1, n2, n3;
    func.vol = func.detJ * func.weight;

    func.K[ii][jj] = 0e0;
    for(int d=0; d<3; d++){
       func.K[ii][jj] += func.dNdx[ii][d] * func.dNdx[jj][d];
    }

    // Mass term
    Klocal(IU, JU) += func.N[ii] * func.N[jj] / dt * func.vol;
    Klocal(IV, JV) += func.N[ii] * func.N[jj] / dt * func.vol;
    Klocal(IW, JW) += func.N[ii] * func.N[jj] / dt * func.vol;

    /*
     // Diffusion term
    for(int d=0; d<3; d++){
        if(d == 0){n1 = 2; n2 = 1; n3 = 1;}
        if(d == 1){n1 = 1; n2 = 2; n3 = 1;}
        if(d == 2){n1 = 1; n2 = 1; n3 = 2;}
        Klocal(IU, JU) += 5e-1 * n1 * func.dNdx[ii][d] * func.dNdx[jj][d] / Re * func.vol;
        Klocal(IV, JV) += 5e-1 * n2 * func.dNdx[ii][d] * func.dNdx[jj][d] / Re * func.vol;
        Klocal(IW, JW) += 5e-1 * n3 * func.dNdx[ii][d] * func.dNdx[jj][d] / Re * func.vol;
    }
    Klocal(IU, JV) += 5e-1 * func.dNdx[ii][1] * func.dNdx[jj][0] / Re * func.vol;
    Klocal(IU, JW) += 5e-1 * func.dNdx[ii][2] * func.dNdx[jj][0] / Re * func.vol;
    Klocal(IV, JU) += 5e-1 * func.dNdx[ii][0] * func.dNdx[jj][1] / Re * func.vol;
    Klocal(IV, JW) += 5e-1 * func.dNdx[ii][2] * func.dNdx[jj][1] / Re * func.vol;
    Klocal(IW, JU) += 5e-1 * func.dNdx[ii][0] * func.dNdx[jj][2] / Re * func.vol;
    Klocal(IW, JV) += 5e-1 * func.dNdx[ii][1] * func.dNdx[jj][2] / Re * func.vol;
    */
    
     // Diffusion term
    for(int d=0; d<3; d++){
        Klocal(IU, JU) += 5e-1 * func.dNdx[ii][d] * func.dNdx[jj][d] / Re * func.vol;
        Klocal(IV, JV) += 5e-1 * func.dNdx[ii][d] * func.dNdx[jj][d] / Re * func.vol;
        Klocal(IW, JW) += 5e-1 * func.dNdx[ii][d] * func.dNdx[jj][d] / Re * func.vol;
    }
    
    // Advection term
    for(int d=0; d<3; d++){
        Klocal(IU, JU) += 5e-1 * func.N[ii] * advgp[d] * func.dNdx[jj][d] * func.vol;
        Klocal(IV, JV) += 5e-1 * func.N[ii] * advgp[d] * func.dNdx[jj][d] * func.vol;
        Klocal(IW, JW) += 5e-1 * func.N[ii] * advgp[d] * func.dNdx[jj][d] * func.vol;
    }

    // Pressure term
    Klocal(IU, JP) -= func.N[jj] * func.dNdx[ii][0] * func.vol;
    Klocal(IV, JP) -= func.N[jj] * func.dNdx[ii][1] * func.vol;
    Klocal(IW, JP) -= func.N[jj] * func.dNdx[ii][2] * func.vol;

    // Continuity term
    Klocal(IP, JU) += func.N[ii] * func.dNdx[jj][0] * func.vol;
    Klocal(IP, JV) += func.N[ii] * func.dNdx[jj][1] * func.vol;
    Klocal(IP, JW) += func.N[ii] * func.dNdx[jj][2] * func.vol;

    // Darcy term
    Klocal(IU, JU) += 5e-1 * f * func.N[ii] * func.N[jj] * func.vol;
    Klocal(IV, JV) += 5e-1 * f * func.N[ii] * func.N[jj] * func.vol;
    Klocal(IW, JW) += 5e-1 * f * func.N[ii] * func.N[jj] * func.vol;

    // SUPG　mass term
    for(int d=0; d<3; d++){
        Klocal(IU, JU) += tau * func.dNdx[ii][d] * advgp[d] * func.N[jj] / dt * func.vol;
        Klocal(IV, JV) += tau * func.dNdx[ii][d] * advgp[d] * func.N[jj] / dt * func.vol;
        Klocal(IW, JW) += tau * func.dNdx[ii][d] * advgp[d] * func.N[jj] / dt * func.vol;
    }

    // SUPG advection term
    for(int d1=0; d1<3; d1++){
        for(int d2=0; d2<3; d2++){
            Klocal(IU, JU) += 5e-1 * tau * advgp[d2] * func.dNdx[ii][d2] * advgp[d1] * func.dNdx[jj][d1] * func.vol;
            Klocal(IV, JV) += 5e-1 * tau * advgp[d2] * func.dNdx[ii][d2] * advgp[d1] * func.dNdx[jj][d1] * func.vol;
            Klocal(IW, JW) += 5e-1 * tau * advgp[d2] * func.dNdx[ii][d2] * advgp[d1] * func.dNdx[jj][d1] * func.vol;
        }
    }

    // SUPG pressure term
    for(int d=0; d<3; d++){
        Klocal(IU, JP) += tau * func.dNdx[ii][d] * advgp[d] * func.dNdx[jj][0] * func.vol;
        Klocal(IV, JP) += tau * func.dNdx[ii][d] * advgp[d] * func.dNdx[jj][1] * func.vol;
        Klocal(IW, JP) += tau * func.dNdx[ii][d] * advgp[d] * func.dNdx[jj][2] * func.vol;
    }

    // PSPG mass term
    Klocal(IP, JU) += tau * func.dNdx[ii][0] * func.N[jj] / dt * func.vol;
    Klocal(IP, JV) += tau * func.dNdx[ii][1] * func.N[jj] / dt * func.vol;
    Klocal(IP, JW) += tau * func.dNdx[ii][2] * func.N[jj] / dt * func.vol;

    // PSPG asvection term
    for(int d=0; d<3; d++){
        Klocal(IP, JU) += 5e-1 * tau * func.dNdx[ii][0] * advgp[d] * func.dNdx[jj][d] * func.vol;
        Klocal(IP, JV) += 5e-1 * tau * func.dNdx[ii][1] * advgp[d] * func.dNdx[jj][d] * func.vol;
        Klocal(IP, JW) += 5e-1 * tau * func.dNdx[ii][2] * advgp[d] * func.dNdx[jj][d] * func.vol;
    }

    // PSPG pressure term
    Klocal(IP, JP) += tau * func.K[ii][jj] * func.vol;

}

/**********************************************************************
 * @brief Compute element sfiffness RHS vector on gauss integral point. 
 */
void DirectProblem::usnsGaussIntegralRHS(VectorXd &Flocal, Function &func, const double f, const int ii)
{
    int n1, n2, n3;
    func.vol = func.detJ * func.weight;

    // Mass term
    Flocal(IU) += func.N[ii] * vgp[0] / dt * func.vol;
    Flocal(IV) += func.N[ii] * vgp[1] / dt * func.vol;
    Flocal(IW) += func.N[ii] * vgp[2] / dt * func.vol;

    /*
    // Diffusion term
    for(int d=0; d<3; d++){
        if(d == 0){n1 = 2; n2 = 1; n3 = 1;}
        if(d == 1){n1 = 1; n2 = 2; n3 = 1;}
        if(d == 2){n1 = 1; n2 = 1; n3 = 2;}
        Flocal(IU) -= 5e-1 * n1 * func.dNdx[ii][d] * dvgpdx[0][d] / Re * func.vol;
        Flocal(IV) -= 5e-1 * n2 * func.dNdx[ii][d] * dvgpdx[1][d] / Re * func.vol;
        Flocal(IW) -= 5e-1 * n3 * func.dNdx[ii][d] * dvgpdx[2][d] / Re * func.vol;
    }
    Flocal(IU) -= 5e-1 * func.dNdx[ii][1] * dvgpdx[1][0] / Re * func.vol;
    Flocal(IU) -= 5e-1 * func.dNdx[ii][2] * dvgpdx[2][0] / Re * func.vol;
    Flocal(IV) -= 5e-1 * func.dNdx[ii][0] * dvgpdx[0][1] / Re * func.vol;
    Flocal(IV) -= 5e-1 * func.dNdx[ii][2] * dvgpdx[2][1] / Re * func.vol;
    Flocal(IW) -= 5e-1 * func.dNdx[ii][0] * dvgpdx[0][2] / Re * func.vol;
    Flocal(IW) -= 5e-1 * func.dNdx[ii][1] * dvgpdx[1][2] / Re * func.vol;
    */

    // Diffusion term
    for(int d=0; d<3; d++){
        Flocal(IU) -= 5e-1 * func.dNdx[ii][d] * dvgpdx[0][d] / Re * func.vol;
        Flocal(IV) -= 5e-1 * func.dNdx[ii][d] * dvgpdx[1][d] / Re * func.vol;
        Flocal(IW) -= 5e-1 * func.dNdx[ii][d] * dvgpdx[2][d] / Re * func.vol;
    }

    // Advection term
    for(int d=0;d<3;d++){
        Flocal(IU) -= 5e-1 * func.N[ii] * advgp[d] * dvgpdx[0][d] * func.vol;
        Flocal(IV) -= 5e-1 * func.N[ii] * advgp[d] * dvgpdx[1][d] * func.vol;
        Flocal(IW) -= 5e-1 * func.N[ii] * advgp[d] * dvgpdx[2][d] * func.vol;
    }

    // Darcy term
    Flocal(IU) -= 5e-1 * f * func.N[ii] * vgp[0] * func.vol;
    Flocal(IV) -= 5e-1 * f * func.N[ii] * vgp[1] * func.vol;
    Flocal(IW) -= 5e-1 * f * func.N[ii] * vgp[2] * func.vol;

    // SUPG mass term
    for(int d=0; d<3; d++){
        Flocal(IU) += tau * func.dNdx[ii][d] * advgp[d] * vgp[0] / dt * func.vol;
        Flocal(IV) += tau * func.dNdx[ii][d] * advgp[d] * vgp[1] / dt * func.vol;
        Flocal(IW) += tau * func.dNdx[ii][d] * advgp[d] * vgp[2] / dt * func.vol;
    }
    
    // SUPG advection term
    for(int d=0; d<3; d++){
        for(int nn=0; nn<3; nn++){
            Flocal(IU) -= 5e-1 * tau * vgp[nn] * func.dNdx[ii][nn] * advgp[d] * dvgpdx[0][d] * func.vol;
            Flocal(IV) -= 5e-1 * tau * vgp[nn] * func.dNdx[ii][nn] * advgp[d] * dvgpdx[1][d] * func.vol;
            Flocal(IW) -= 5e-1 * tau * vgp[nn] * func.dNdx[ii][nn] * advgp[d] * dvgpdx[2][d] * func.vol;
        }
    }
    
    // PSPG mass term
    Flocal(IP) += tau * func.dNdx[ii][0] * vgp[0] / dt * func.vol;
    Flocal(IP) += tau * func.dNdx[ii][1] * vgp[1] / dt * func.vol;
    Flocal(IP) += tau * func.dNdx[ii][2] * vgp[2] / dt * func.vol;

    // PSPG advection term
    for(int d=0; d<3; d++){
        Flocal(IP) -= 5e-1 * tau * func.dNdx[ii][0] * advgp[d] * dvgpdx[0][d] * func.vol;
        Flocal(IP) -= 5e-1 * tau * func.dNdx[ii][1] * advgp[d] * dvgpdx[1][d] * func.vol;
        Flocal(IP) -= 5e-1 * tau * func.dNdx[ii][2] * advgp[d] * dvgpdx[2][d] * func.vol;
    }
}


void DirectProblem::setVelocityValue(Function &func, const int ic, const int t)
{
    for(int d=0; d<dim; d++){
        advgp[d] = 0e0;
        vgp[d] = 0e0;
        for(int p=0; p<grid.cell.nNodesInCell; p++){
            if(t == 0){
                advgp[d] += func.N[p] * grid.node.v[grid.cell(ic).node[p]][d];
            }else{
                advgp[d] += func.N[p] * (1.5 * grid.node.v[grid.cell(ic).node[p]][d] - 0.5 * grid.node.vPrev[grid.cell(ic).node[p]][d]);
            }
            vgp[d] += func.N[p] * grid.node.v[grid.cell(ic).node[p]][d];
        }
    }
    for(int d=0; d<dim; d++){
        for(int e=0; e<dim; e++){
            dvgpdx[d][e] = 0e0;
            for(int p=0; p<grid.cell.nNodesInCell; p++){
                dvgpdx[d][e] += func.dNdx[p][e] * grid.node.v[grid.cell(ic).node[p]][d];
            }
        }
    }
}

