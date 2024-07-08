#include "InverseProblem.h"

void Adjoint::solveAdjoint(DirectProblem &main, std::string outputDir,
                           std::vector<std::vector<std::vector<double>>> &feedbackForceT,
                           const int nData, const int loop)
{
    PetscPrintf(MPI_COMM_WORLD, "\nADJOINT SOLVER\n");

    PetscScalar *arraySolnTmp;
    Vec  vecSEQ;
    VecScatter  ctx;
    VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

    for(int id=0; id<grid.nDofsGlobal; id++){
        grid.dirichlet.dirichletBCsValueNewInit[id] = 0e0;
        grid.dirichlet.dirichletBCsValueNew[id] = 0e0;
    }

    petsc.setMatAndVecZero(grid.cell);
    petsc.initialAssembly();
    setVariablesZero(main.dim);

    int st = 0;
    for(int t=0; t<timeMax; t++){
        petsc.setValueZero();
        grid.dirichlet.assignDirichletBCs(grid.dirichlet.vDirichletWallNew, 
                                          grid.dirichlet.pDirichletNew, 
                                          grid.node, main.dim, t);
        grid.dirichlet.applyDirichletBCsAdjoint(grid.cell, petsc);

        for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                MatrixXd  Klocal(nDofsInCell, nDofsInCell);
                VectorXd  Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();
                matrixAssemblyAdjointUSNS(main, Klocal, Flocal, ic, t);
                petsc.setValue(grid.cell(ic).dofsBCsMapWall, grid.cell(ic).dofsMap,
                               grid.cell(ic).dofsBCsMapWall, Klocal, Flocal);
            }
        }

        for(int ib=0; ib<grid.dirichlet.controlCellMap.size(); ib++){
            int ic = grid.dirichlet.controlCellMap[ib];
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                MatrixXd  Klocal(nDofsInCell, nDofsInCell);
                VectorXd  Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();
                boundaryIntegral(main, Klocal, Flocal, ic, ib);
                petsc.setMatValue(grid.cell(ic).dofsBCsMapWall, grid.cell(ic).dofsMap, Klocal);
            }
        }

        for(int in=0; in<grid.node.nNodesGlobal; in++){
            int n = grid.node.mapNew[in];
            int size = grid.node.dofsBCsMapNew[n].size();
            VectorXd  Flocal(size);
            Flocal.setZero();
            if(grid.node.subId[in] == mpi.myId){
                Flocal(0) -= feedbackForceT[t][in][0];
                Flocal(1) -= feedbackForceT[t][in][1];
                Flocal(2) -= feedbackForceT[t][in][2];
                petsc.setVecValue(grid.node.dofsBCsMapWallNew[n], Flocal);
            }
        }

        petsc.solve();

        VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecGetArray(vecSEQ, &arraySolnTmp);

        // update solution vector
        for(int id=0; id<grid.nDofsGlobal; id++)
            petsc.solution[id] = arraySolnTmp[id];

        VecRestoreArray(vecSEQ, &arraySolnTmp);
        updateVariables(outputDir, main.dim, t, loop);

        if(mpi.myId == 0){
            double timeNow = t * dt;
            printf("Adjoint Solver : Time = %f \n", timeNow);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    VecScatterDestroy(&ctx);
    VecDestroy(&vecSEQ);
}

void Adjoint::solveAdjointDO(DirectProblem &main, std::string outputDir,
                             std::vector<std::vector<std::vector<double>>> &feedbackForceT,
                             const int nData, const int loop)
{
    PetscPrintf(MPI_COMM_WORLD, "\nADJOINT SOLVER\n");

    PetscScalar *arraySolnTmp;
    Vec  vecSEQ;
    VecScatter  ctx;
    VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

    for(int id=0; id<grid.nDofsGlobal; id++){
        grid.dirichlet.dirichletBCsValueNewInit[id] = 0e0;
        grid.dirichlet.dirichletBCsValueNew[id] = 0e0;
    }
    petsc.setMatAndVecZero(grid.cell);
    petsc.initialAssembly();
    setVariablesZero(main.dim);

    int st = 0;
    for(int t=timeMax-1; t>=0; t--){
        petsc.setValueZero();
        grid.dirichlet.assignDirichletBCs(grid.dirichlet.vDirichletWallNew, 
                                          grid.dirichlet.pDirichletNew, 
                                          grid.node, main.dim, t);
        grid.dirichlet.applyDirichletBCsAdjoint(grid.cell, petsc);

        for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                MatrixXd  Klocal(nDofsInCell, nDofsInCell);
                VectorXd  Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();
                matrixAssemblyAdjoint_DO(main, Klocal, Flocal, ic, t);
                petsc.setValue(grid.cell(ic).dofsBCsMapWall, grid.cell(ic).dofsMap,
                               grid.cell(ic).dofsBCsMapWall, Klocal, Flocal);
            }
        }

        for(int ib=0; ib<grid.dirichlet.controlCellMap.size(); ib++){
            int ic = grid.dirichlet.controlCellMap[ib];
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                MatrixXd  Klocal(nDofsInCell, nDofsInCell);
                VectorXd  Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();
                boundaryIntegral(main, Klocal, Flocal, ic, ib);
                petsc.setMatValue(grid.cell(ic).dofsBCsMapWall, grid.cell(ic).dofsMap, Klocal);
            }
        }

        for(int in=0; in<grid.node.nNodesGlobal; in++){
            int n = grid.node.mapNew[in];
            int size = grid.node.dofsBCsMapNew[n].size();
            VectorXd  Flocal(size);
            Flocal.setZero();
            if(grid.node.subId[in] == mpi.myId){
                Flocal(0) -= feedbackForceT[t][in][0];
                Flocal(1) -= feedbackForceT[t][in][1];
                Flocal(2) -= feedbackForceT[t][in][2];
                petsc.setVecValue(grid.node.dofsBCsMapWallNew[n], Flocal);
            }
        }

        petsc.solve();

        VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecGetArray(vecSEQ, &arraySolnTmp);

        // update solution vector
        for(int id=0; id<grid.nDofsGlobal; id++)
            petsc.solution[id] = arraySolnTmp[id];

        VecRestoreArray(vecSEQ, &arraySolnTmp);
        updateVariables(outputDir, main.dim, t, loop);

        if(mpi.myId == 0){
            double timeNow = t * dt;
            printf("Adjoint Solver : Time = %f \n", timeNow);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    VecScatterDestroy(&ctx);
    VecDestroy(&vecSEQ);
}


void Adjoint::setVariablesZero(const int dim)
{
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            grid.node.w[in][d] = 0e0;
            grid.node.wPrev[in][d] = 0e0;
            grid.node.l[in][d] = 0e0;
        }
        grid.node.q[in] = 0e0;
    }

    for(int t=0; t<timeMax; t++){
        for(int in=0; in<grid.node.nNodesGlobal; in++){
            for(int d=0; d<dim; d++){
                grid.node.lt[t][in][d] = 0e0;
            }
        }
    }
}

void Adjoint::boundaryIntegral(DirectProblem &main, MatrixXd &Klocal, VectorXd &Flocal,
                               const int ic, const int ib)
{
    int IU, ILU, ILV, ILW;
    int JU, JLU, JLV, JLW;
    int IV, IW, IP;
    int JV, JW, JP;
    
    int nControlNodesInCell = 4;
    double dxdr[2][2]; 
    std::vector<double> N2D;
    std::vector<std::vector<double>> dNdr2D;
    std::vector<std::vector<double>> dNdx2D;
    std::vector<std::vector<double>> xCurrent2D;
    std::vector<std::vector<double>> K;
    
    VecTool::resize(N2D, nControlNodesInCell);
    VecTool::resize(dNdr2D, nControlNodesInCell, main.dim-1);
    VecTool::resize(dNdx2D,  nControlNodesInCell, main.dim-1);
    VecTool::resize(xCurrent2D, nControlNodesInCell, main.dim-1);
    VecTool::resize(K, nControlNodesInCell, nControlNodesInCell);

    Gauss gauss(2);
    double value;

    for(int p=0; p<nControlNodesInCell; p++){
        int n = grid.dirichlet.controlNodeInCell[ib][p];
        for(int d=0; d<main.dim-1; d++){
            xCurrent2D[p][d] = main.grid.node.x[n][planeDir[d]];
        }
    }
    double he = fabs(xCurrent2D[1][0] - xCurrent2D[0][0]);
    
    int s[grid.cell.nNodesInCell];
    for(int p=0; p<grid.cell.nNodesInCell; p++){
        s[p] = 0;
    }       
    for(int p=1; p<grid.cell.nNodesInCell; p++){
        s[p] = s[p-1] + grid.node.nDofsOnNode[grid.cell(ic).node[p-1]]; 
    }

    int ss[nControlNodesInCell];
    int count = 0;
    for(int p=0; p<grid.cell.nNodesInCell; p++){
        if(grid.node.nDofsOnNode[grid.cell(ic).node[p]] > main.dim+1){
            ss[count] = s[p];
            count++;
        }
    }

    double wgp[3];
    double lgp[3];

    int tmp = ss[2];
    ss[2] = ss[3];
    ss[3] = tmp;

    for(int i1=0; i1<2; i1++){
        for(int i2=0; i2<2; i2++){
            ShapeFunction2D::C2D4_N(N2D, gauss.point[i1], gauss.point[i2]);
            ShapeFunction2D::C2D4_dNdr(dNdr2D, gauss.point[i1], gauss.point[i2]);
            MathFEM::comp_dxdr2D(dxdr, dNdr2D, xCurrent2D, nControlNodesInCell);
            MathFEM::comp_dNdx2D(dNdx2D, dNdr2D, dxdr, nControlNodesInCell);

            double detJ = dxdr[0][0] * dxdr[1][1] - dxdr[0][1] * dxdr[1][0];
            double weight = gauss.weight[i1] * gauss.weight[i2];
            double tau = MathFEM::comp_tau(wgp, he, main.Re, main.dt);

            for(int ii=0; ii<count; ii++){
                IU = ss[ii];  IV = IU + 1;  IW = IU + 2;
                ILU = IU + 4; ILV = IU + 5; ILW = IU + 6;
                for(int jj=0; jj<count; jj++){
                    JU = ss[jj];  JV = JU + 1;  JW = JU + 2;
                    JLU = JU + 4; JLV = JU + 5; JLW = JU + 6;
                    Klocal(IU, JLU) -= N2D[ii] * N2D[jj] * detJ * weight;
                    Klocal(IV, JLV) -= N2D[ii] * N2D[jj] * detJ * weight;
                    Klocal(IW, JLW) -= N2D[ii] * N2D[jj] * detJ * weight;
                    Klocal(ILU, JU) -= N2D[ii] * N2D[jj] * detJ * weight; 
                    Klocal(ILV, JV) -= N2D[ii] * N2D[jj] * detJ * weight; 
                    Klocal(ILW, JW) -= N2D[ii] * N2D[jj] * detJ * weight; 
                }
            }
            
        }
    }
}

void Adjoint::updateVariables(std::string outputDir, const int dim, const int t, const int loop)
{
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        int n1 = 0;
        for(int i=0; i<grid.node.mapNew[in]; i++)
            n1 += grid.node.nDofsOnNodeNew[i];

        for(int d=0; d<dim; d++)
            grid.node.wPrev[in][d] = grid.node.w[in][d];

        for(int d=0; d<dim; d++)
            grid.node.w[in][d] = petsc.solution[n1+d];
        grid.node.q[in] = petsc.solution[n1+dim];

        if(grid.node.nDofsOnNodeNew[grid.node.mapNew[in]] > dim+1){
            for(int d=0; d<dim; d++)
                grid.node.l[in][d] = petsc.solution[n1+dim+1+d];
        }
    }
    
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            grid.node.lt[t][in][d] = grid.node.l[in][d];
        }
    }

    if(mpi.myId == 0){
        std::string vtuFile;
        vtuFile = outputDir + "/solution/w_" + to_string(loop) + "_" + to_string(t) + ".vtu";
        grid.output.exportAdjointSolutionVTU(vtuFile, grid.node, grid.cell, DataType::ADJOINT_W);
        vtuFile = outputDir + "/solution/q_" + to_string(loop) + "_" + to_string(t) + ".vtu";
        grid.output.exportAdjointSolutionVTU(vtuFile, grid.node, grid.cell, DataType::ADJOINT_Q);
        vtuFile = outputDir + "/solution/l_" + to_string(loop) + "_" + to_string(t) + ".vtu";
        grid.output.exportAdjointSolutionVTU(vtuFile, grid.node, grid.cell, DataType::ADJOINT_L);
    }
}
