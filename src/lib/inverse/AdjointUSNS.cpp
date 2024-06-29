#include "InverseProblem.h"

void Adjoint::solveAdjointEquation(DirectProblem &main, std::string outputDir,
                                   std::vector<std::vector<std::vector<double>>> &feedbackForce)
{
    PetscPrintf(MPI_COMM_WORLD, "\nADJOINT UNSTEADY NAVIER STOKES SOLVER\n");
    
    main.nu = main.mu / main.rho;
    main.Re = 1e0 / main.nu;  
    PetscScalar *arraySolnTmp;
    Vec  vecSEQ;
    VecScatter  ctx;
    VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

    petsc.solution.resize(grid.nDofsGlobal);
    grid.dirichlet.dirichletBCsValue.resize(grid.nDofsGlobal, 0e0);
    grid.dirichlet.dirichletBCsValueNew.resize(grid.nDofsGlobal, 0e0);
    grid.dirichlet.dirichletBCsValueInit.resize(grid.nDofsGlobal, 0e0);
    grid.dirichlet.dirichletBCsValueNewInit.resize(grid.nDofsGlobal, 0e0);

    grid.dirichlet.assignDirichletBCs(grid.node, main.dim);

    petsc.setMatAndVecZero(grid.cell);
    petsc.initialAssembly();

    int st = 0;
    for(int t=0; t<main.timeMax; t++){
        petsc.setValueZero();
        //grid.dirichlet.applyDirichletBCs(grid.cell, petsc);

        MPI_Barrier(MPI_COMM_WORLD);
        double timer = MPI_Wtime();
        
        for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                MatrixXd  Klocal(nDofsInCell, nDofsInCell);
                VectorXd  Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();
                matrixAssemblyAdjointUSNS(main, Klocal, Flocal, feedbackForce, st, ic, t);
                petsc.setValue(grid.cell(ic).dofsMap, grid.cell(ic).dofsMap,
                               grid.cell(ic).dofsMap, Klocal, Flocal);
            }
        }
        for(int ib=0; ib<grid.dirichlet.controlCellMap.size(); ib++){
            int ic = grid.dirichlet.controlCellMap[ib];
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                MatrixXd  Klocal(nDofsInCell, nDofsInCell);
                Klocal.setZero();
                boundaryIntegral(main, Klocal, ic, ib);
                petsc.setMatValue(grid.cell(ic).dofsMap, 
                                  grid.cell(ic).dofsMap, Klocal);
            }
        }

        timer = MPI_Wtime() - timer;
        PetscPrintf(MPI_COMM_WORLD, "\nAdjoint Matrix assembly = %f seconds\n", timer);
        petsc.currentStatus = ASSEMBLY_OK;
    
        MPI_Barrier(MPI_COMM_WORLD); 
        timer = MPI_Wtime();
        petsc.solve();
        timer = MPI_Wtime() - timer;

        PetscPrintf(MPI_COMM_WORLD, "\nAdjoint PETSc solver = %f seconds \n", timer);

        VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecGetArray(vecSEQ, &arraySolnTmp);

        // update solution vector
        for(int id=0; id<grid.nDofsGlobal; id++)
            petsc.solution[id] = arraySolnTmp[id];

        VecRestoreArray(vecSEQ, &arraySolnTmp);
        updateVariables(outputDir, main.dim, t);

        if(mpi.myId == 0){
            double timeNow = t * main.dt;
            printf("\n TIME = %f \n", timeNow);
        }
        if(t % main.snap.snapInterval == 0){
            st++;
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    VecScatterDestroy(&ctx);
    VecDestroy(&vecSEQ);
}

void Adjoint::boundaryIntegral(DirectProblem &main, MatrixXd &Klocal, 
                               const int ic, const int ib)
{
    int IU, ILU, ILV, ILW;
    int JU, JLU, JLV, JLW;
    int IV, IW, IP;
    int JV, JW, JP;
    
    int nControlNodesInCell = 4;

    std::vector<double> N2D;
    std::vector<std::vector<double>> dNdr2D;
    std::vector<std::vector<double>> xCurrent2D;

    N2D.resize(nControlNodesInCell, 0e0);
    dNdr2D.resize(nControlNodesInCell, std::vector<double>(main.dim-1, 0e0));
    xCurrent2D.resize(nControlNodesInCell, std::vector<double>(main.dim-1, 0e0));

    Gauss gauss(2);
    double value;

    for(int p=0; p<nControlNodesInCell; p++){
        int index = grid.dirichlet.controlNodeInCell[ib][p];
        for(int d=0; d<main.dim-1; d++){
            xCurrent2D[p][d] = main.grid.node.x[index][planeDir[d]];
        }
    }
    
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

    /*
    int tmp = ss[3];
    ss[3] = ss[4];
    ss[4] = tmp;
    */

    for(int i1=0; i1<2; i1++){
        for(int i2=0; i2<2; i2++){
            ShapeFunction2D::C2D4_N(N2D, gauss.point[i1], gauss.point[i2]);
            ShapeFunction2D::C2D4_dNdr(dNdr2D, gauss.point[i1], gauss.point[i2]);
            double dxdr[2][2]; 
            MathFEM::calc_dxdr2D(dxdr, dNdr2D, xCurrent2D, nControlNodesInCell);
            double detJ = dxdr[0][0] * dxdr[1][1] - dxdr[0][1] * dxdr[1][0];
            double weight = gauss.weight[i1] * gauss.weight[i2];
            for(int ii=0; ii<nControlNodesInCell; ii++){
                IU = ss[ii];
                IV = IU + 1;
                IW = IU + 2;
                ILU = IU + 4;
                ILV = IU + 5;
                ILW = IU + 6;
                for(int jj=0; jj<nControlNodesInCell; jj++){
                    JU = ss[jj];
                    JV = JU + 1;
                    JW = JU + 2;
                    JLU = JU + 4;
                    JLV = JU + 5;
                    JLW = JU + 6;
                    Klocal(IU, JLU) -= N2D[ii] * N2D[jj] * detJ * weight;
                    Klocal(IV, JLV) -= N2D[ii] * N2D[jj] * detJ * weight;
                    Klocal(IW, JLW) -= N2D[ii] * N2D[jj] * detJ * weight;
                    Klocal(ILU, JU) -= N2D[ii] * N2D[jj] * detJ * weight; 
                    Klocal(ILV, JV) -= N2D[ii] * N2D[jj] * detJ * weight; 
                    Klocal(ILW, JW) -= N2D[ii] * N2D[jj] * detJ * weight; 
                    //Klocal(ILU, ILU) -= N2D[jj] * N2D[ii] * detJ * weight;
                    //Klocal(ILV, ILV) -= N2D[jj] * N2D[ii] * detJ * weight; 
                    //Klocal(ILW, ILW) -= N2D[jj] * N2D[ii] * detJ * weight;  
                }
            }
            
        }
    }
}

void Adjoint::updateVariables(std::string outputDir, const int dim, const int t)
{
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        int n1 = 0;
        for(int i=0; i<grid.node.mapNew[in]; i++)
            n1 += grid.node.nDofsOnNodeNew[i];

        grid.node.p[in] = 0e0;

        for(int d=0; d<dim; d++)
            grid.node.v[in][d] = petsc.solution[n1+d];
        grid.node.p[in] = petsc.solution[n1+dim];

        if(grid.node.nDofsOnNodeNew[grid.node.mapNew[in]] > dim+1){
            for(int d=0; d<dim; d++)
                grid.node.lambda[in][d] = petsc.solution[n1+dim+1+d];
        }
    }

    for(int in=0; in<grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            grid.node.vt[t][in][d] =  grid.node.v[in][d];
            grid.node.lambdat[t][in][d] = grid.node.lambda[in][d];
        }
        grid.node.pt[t][in] = grid.node.p[in];
    }

    if(mpi.myId == 0){
        std::string vtuFile;
        vtuFile = outputDir + "/velocity/velocity" + to_string(t) + ".vtu";
        grid.output.exportSolutionVTU(vtuFile, grid.node, grid.cell, DataType::VELOCITY);
        vtuFile = outputDir + "/pressure/pressure" + to_string(t) + ".vtu";
        grid.output.exportSolutionVTU(vtuFile, grid.node, grid.cell, DataType::PRESSURE);
        vtuFile = outputDir + "/lambda/lambda" + to_string(t) + ".vtu";
        grid.output.exportSolutionVTU(vtuFile, grid.node, grid.cell, DataType::LAMBDA);
    }
}
