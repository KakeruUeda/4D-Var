#include "DirectProblem.h"

void DirectProblem::matrixAssemblyUSNS(MatrixXd &Klocal, VectorXd &Flocal, 
                                       const int ic, const int t)
{
    for(int p=0; p<grid.cell.nNodesInCell; p++){
       for(int d=0; d<dim; d++){
            xCurrent[p][d] = 0e0;
            xCurrent[p][d] = grid.node.x[grid.cell(ic).node[p]][d];
       } 
    }       
             
    double he = fabs(xCurrent[1][0] - xCurrent[0][0]);
    double f = resistance * alpha * (1e0 - grid.cell(ic).phi) / (alpha + grid.cell(ic).phi);

    Gauss g2(2);

    for(int i1=0; i1<2; i1++){
        for(int i2=0; i2<2; i2++){
            for(int i3=0; i3<2; i3++){
                ShapeFunction3D::C3D8_N(N, g2.point[i1], g2.point[i2], g2.point[i3]);
                ShapeFunction3D::C3D8_dNdr(dNdr, g2.point[i1], g2.point[i2], g2.point[i3]);

                MathFEM::comp_dxdr(dxdr, dNdr, xCurrent, grid.cell.nNodesInCell);
                MathFEM::comp_dNdx(dNdx, dNdr, dxdr, grid.cell.nNodesInCell);
                
                detJ = MathCommon::compDeterminant_3x3(dxdr);
                weight = g2.weight[i1] * g2.weight[i2] * g2.weight[i3];

                setVelocityValue(ic, t);
                tau = comp_tau(advgp, he, Re, dt);

                for(int ii=0; ii<grid.cell.nNodesInCell; ii++){  
                    updateRowIndex(ii, ic);
                    for(int jj=0; jj<grid.cell.nNodesInCell; jj++){  
                        updateColumnIndex(jj, ic);
                        mainGaussIntegralLHS(Klocal, f, ii, jj);
                    }
                    mainGaussIntegralRHS(Flocal, f, ii);    
                }
            }
        }
    }

}


void DirectProblem::mainGaussIntegralLHS(MatrixXd &Klocal, const double f, const int ii, const int jj)
{
    int n1, n2, n3;

    K[ii][jj] = 0e0;
    for(int d=0; d<3; d++){
        K[ii][jj] += dNdx[ii][d] * dNdx[jj][d];
    }

    // MASS TERM
    Klocal(IU, JU) += N[ii] * N[jj] / dt * detJ * weight;
    Klocal(IV, JV) += N[ii] * N[jj] / dt * detJ * weight;
    Klocal(IW, JW) += N[ii] * N[jj] / dt * detJ * weight;
    
     // DIFFUSION TERM
    for(int d=0; d<3; d++){
        if(d == 0){n1 = 2; n2 = 1; n3 = 1;}
        if(d == 1){n1 = 1; n2 = 2; n3 = 1;}
        if(d == 2){n1 = 1; n2 = 1; n3 = 2;}
        Klocal(IU, JU) += 5e-1 * n1 * dNdx[ii][d] * dNdx[jj][d] / Re * detJ * weight;
        Klocal(IV, JV) += 5e-1 * n2 * dNdx[ii][d] * dNdx[jj][d] / Re * detJ * weight;
        Klocal(IW, JW) += 5e-1 * n3 * dNdx[ii][d] * dNdx[jj][d] / Re * detJ * weight;
    }
    Klocal(IU, JV) += 5e-1 * dNdx[ii][1] * dNdx[jj][0] / Re * detJ * weight;
    Klocal(IU, JW) += 5e-1 * dNdx[ii][2] * dNdx[jj][0] / Re * detJ * weight;
    Klocal(IV, JU) += 5e-1 * dNdx[ii][0] * dNdx[jj][1] / Re * detJ * weight;
    Klocal(IV, JW) += 5e-1 * dNdx[ii][2] * dNdx[jj][1] / Re * detJ * weight;
    Klocal(IW, JU) += 5e-1 * dNdx[ii][0] * dNdx[jj][2] / Re * detJ * weight;
    Klocal(IW, JV) += 5e-1 * dNdx[ii][1] * dNdx[jj][2] / Re * detJ * weight;

    // ADVECTION TERM
    for(int d=0; d<3; d++){
        Klocal(IU, JU) += 5e-1 * N[ii] * advgp[d] * dNdx[jj][d] * detJ * weight;
        Klocal(IV, JV) += 5e-1 * N[ii] * advgp[d] * dNdx[jj][d] * detJ * weight;
        Klocal(IW, JW) += 5e-1 * N[ii] * advgp[d] * dNdx[jj][d] * detJ * weight;
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
    for(int d=0; d<3; d++){
        Klocal(IU, JU) += tau * dNdx[ii][d] * advgp[d] * N[jj] / dt * detJ * weight;
        Klocal(IV, JV) += tau * dNdx[ii][d] * advgp[d] * N[jj] / dt * detJ * weight;
        Klocal(IW, JW) += tau * dNdx[ii][d] * advgp[d] * N[jj] / dt * detJ * weight;
    }

    // ADVECTION TERM
    for(int d1=0; d1<3; d1++){
        for(int d2=0; d2<3; d2++){
            Klocal(IU, JU) += 5e-1 * tau * advgp[d2] * dNdx[ii][d2] * advgp[d1] * dNdx[jj][d1] * detJ * weight;
            Klocal(IV, JV) += 5e-1 * tau * advgp[d2] * dNdx[ii][d2] * advgp[d1] * dNdx[jj][d1] * detJ * weight;
            Klocal(IW, JW) += 5e-1 * tau * advgp[d2] * dNdx[ii][d2] * advgp[d1] * dNdx[jj][d1] * detJ * weight;
        }
    }

    // PRESSURE TERM
    for(int d=0; d<3; d++){
        Klocal(IU, JP) += tau * dNdx[ii][d] * advgp[d] * dNdx[jj][0] * detJ * weight;
        Klocal(IV, JP) += tau * dNdx[ii][d] * advgp[d] * dNdx[jj][1] * detJ * weight;
        Klocal(IW, JP) += tau * dNdx[ii][d] * advgp[d] * dNdx[jj][2] * detJ * weight;
    }

    // PSPG TERM
    // MASS TERM 
    Klocal(IP, JU) += tau * dNdx[ii][0] * N[jj] / dt * detJ * weight;
    Klocal(IP, JV) += tau * dNdx[ii][1] * N[jj] / dt * detJ * weight;
    Klocal(IP, JW) += tau * dNdx[ii][2] * N[jj] / dt * detJ * weight;

    // ADVECTION TERM
    for(int d=0; d<3; d++){
        Klocal(IP, JU) += 5e-1 * tau * dNdx[ii][0] * advgp[d] * dNdx[jj][d] * detJ * weight;
        Klocal(IP, JV) += 5e-1 * tau * dNdx[ii][1] * advgp[d] * dNdx[jj][d] * detJ * weight;
        Klocal(IP, JW) += 5e-1 * tau * dNdx[ii][2] * advgp[d] * dNdx[jj][d] * detJ * weight;
    }

    // PRESSURE TERM
    Klocal(IP, JP) += tau * K[ii][jj] * detJ * weight;

}


void DirectProblem::mainGaussIntegralRHS(VectorXd &Flocal, const double f, const int ii)
{
    int n1, n2, n3;

    // MASS TERM
    Flocal(IU) += N[ii] * vgp[0] / dt * detJ * weight;
    Flocal(IV) += N[ii] * vgp[1] / dt * detJ * weight;
    Flocal(IW) += N[ii] * vgp[2] / dt * detJ * weight;

    // DIFFUSION TERM
    for(int d=0; d<3; d++){
        if(d == 0){n1 = 2; n2 = 1; n3 = 1;}
        if(d == 1){n1 = 1; n2 = 2; n3 = 1;}
        if(d == 2){n1 = 1; n2 = 1; n3 = 2;}
        Flocal(IU) -= 5e-1 * n1 * dNdx[ii][d] * dvgpdx[0][d] / Re * detJ * weight;
        Flocal(IV) -= 5e-1 * n2 * dNdx[ii][d] * dvgpdx[1][d] / Re * detJ * weight;
        Flocal(IW) -= 5e-1 * n3 * dNdx[ii][d] * dvgpdx[2][d] / Re * detJ * weight;
    }
    Flocal(IU) -= 5e-1 * dNdx[ii][1] * dvgpdx[1][0] / Re * detJ * weight;
    Flocal(IU) -= 5e-1 * dNdx[ii][2] * dvgpdx[2][0] / Re * detJ * weight;
    Flocal(IV) -= 5e-1 * dNdx[ii][0] * dvgpdx[0][1] / Re * detJ * weight;
    Flocal(IV) -= 5e-1 * dNdx[ii][2] * dvgpdx[2][1] / Re * detJ * weight;
    Flocal(IW) -= 5e-1 * dNdx[ii][0] * dvgpdx[0][2] / Re * detJ * weight;
    Flocal(IW) -= 5e-1 * dNdx[ii][1] * dvgpdx[1][2] / Re * detJ * weight;

    // ADVECTION TERM
    for(int d=0;d<3;d++){
        Flocal(IU) -= 5e-1 * N[ii] * advgp[d] * dvgpdx[0][d] * detJ * weight;
        Flocal(IV) -= 5e-1 * N[ii] * advgp[d] * dvgpdx[1][d] * detJ * weight;
        Flocal(IW) -= 5e-1 * N[ii] * advgp[d] * dvgpdx[2][d] * detJ * weight;
    }

    // DARCY TERM
    Flocal(IU) -= 5e-1 * f * N[ii] * vgp[0] * detJ * weight;
    Flocal(IV) -= 5e-1 * f * N[ii] * vgp[1] * detJ * weight;
    Flocal(IW) -= 5e-1 * f * N[ii] * vgp[2] * detJ * weight;

    // SUPG TERM
    // MASS TERM
    for(int d=0; d<3; d++){
        Flocal(IU) += tau * dNdx[ii][d] * advgp[d] * vgp[0] / dt * detJ * weight;
        Flocal(IV) += tau * dNdx[ii][d] * advgp[d] * vgp[1] / dt * detJ * weight;
        Flocal(IW) += tau * dNdx[ii][d] * advgp[d] * vgp[2] / dt * detJ * weight;
    }
    
    // ADVECTION TERM
    for(int d=0; d<3; d++){
        for(int nn=0; nn<3; nn++){
            Flocal(IU) -= 5e-1 * tau * vgp[nn] * dNdx[ii][nn] * advgp[d] * dvgpdx[0][d] * detJ * weight;
            Flocal(IV) -= 5e-1 * tau * vgp[nn] * dNdx[ii][nn] * advgp[d] * dvgpdx[1][d] * detJ * weight;
            Flocal(IW) -= 5e-1 * tau * vgp[nn] * dNdx[ii][nn] * advgp[d] * dvgpdx[2][d] * detJ * weight;
        }
    }
    
    // PSPG TERM
    // MASS TERM 
    Flocal(IP) += tau * dNdx[ii][0] * vgp[0] / dt * detJ * weight;
    Flocal(IP) += tau * dNdx[ii][1] * vgp[1] / dt * detJ * weight;
    Flocal(IP) += tau * dNdx[ii][2] * vgp[2] / dt * detJ * weight;

    // ADVECTION TERM
    for(int d=0; d<3; d++){
        Flocal(IP) -= 5e-1 * tau * dNdx[ii][0] * advgp[d] * dvgpdx[0][d] * detJ * weight;
        Flocal(IP) -= 5e-1 * tau * dNdx[ii][1] * advgp[d] * dvgpdx[1][d] * detJ * weight;
        Flocal(IP) -= 5e-1 * tau * dNdx[ii][2] * advgp[d] * dvgpdx[2][d] * detJ * weight;
    }
}



void DirectProblem::setVelocityValue(const int ic, const int t)
{
    for(int d=0; d<dim; d++){
        advgp[d] = 0e0;
        vgp[d] = 0e0;
        for(int p=0; p<grid.cell.nNodesInCell; p++){
            advgp[d] += N[p] * (1.5 * grid.node.v[grid.cell(ic).node[p]][d] - 0.5 * grid.node.vPrev[grid.cell(ic).node[p]][d]);
            vgp[d] += N[p] * grid.node.v[grid.cell(ic).node[p]][d];
        }
    }
    for(int d=0; d<dim; d++){
        for(int e=0; e<dim; e++){
            dvgpdx[d][e] = 0e0;
            for(int p=0; p<grid.cell.nNodesInCell; p++){
                dvgpdx[d][e] += dNdx[p][e] * grid.node.v[grid.cell(ic).node[p]][d];
            }
        }
    }
}

