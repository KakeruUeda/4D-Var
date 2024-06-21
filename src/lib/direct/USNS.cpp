#include "DirectProblem.h"

void DirectProblem::solveUSNS()
{
    PetscPrintf(MPI_COMM_WORLD, "\nUNSTEADY NAVIER STOKES SOLVER\n");
    
    PetscScalar *arraySolnTmp;
    Vec  vecSEQ;
    VecScatter  ctx;
    VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

    petsc.solution.resize(grid.nDofsGlobal);
    grid.dirichlet.dirichletBCsValue.resize(grid.nDofsGlobal, 0e0);
    grid.dirichlet.dirichletBCsValueNew.resize(grid.nDofsGlobal, 0e0);

    grid.dirichlet.assignDirichletBCs(grid.node, dim);
    petsc.setMatAndVecZero(grid.cell);
    petsc.initialAssembly();

    // debug
    if(mpi.myId == 0){
        std::ofstream outDirichletBCsValue(outputDir
                      + "/dat/dirichletBCsValue.dat");
        for(int id=0; id<grid.nDofsGlobal; id++){
            outDirichletBCsValue << grid.dirichlet.dirichletBCsValue[id]
                                 << std::endl;
        }
        outDirichletBCsValue.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double timer = MPI_Wtime();

    for(int tItr=0; tItr<timeMax; tItr++){
        petsc.setValueZero();
        grid.dirichlet.applyDirichletBCs(grid.cell, petsc);
        for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                MatrixXd  Klocal(nDofsInCell, nDofsInCell);
                VectorXd  Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();

                if(grid.cell(ic).phi > 0.999)
                    matrixAssemblyUSNS(Klocal, Flocal, ic, tItr);
                else
                    DarcyMatrixAssemblyUSNS(Klocal, Flocal, ic, tItr);

                petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap,
                               Klocal, Flocal);
            }
        }
        timer = MPI_Wtime() - timer;
        PetscPrintf(MPI_COMM_WORLD, "\n ****** Time for matrix assembly = %f seconds ****** \n", timer);
        petsc.currentStatus = ASSEMBLY_OK;
    
        MPI_Barrier(MPI_COMM_WORLD); 
        timer = MPI_Wtime();

        petsc.solve();

        timer = MPI_Wtime() - timer;
        PetscPrintf(MPI_COMM_WORLD, "\n ****** Time for PETSc solver = %f seconds ****** \n", timer);

        VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecGetArray(vecSEQ, &arraySolnTmp);

        // update solution vector
        for(int id=0; id<grid.nDofsGlobal; id++){
            petsc.solution[id] += arraySolnTmp[id];
        }

        VecRestoreArray(vecSEQ, &arraySolnTmp);
        updataValiables(tItr);

        if(mpi.myId == 0){
            double timeNow = tItr * dt;
            printf("\n TIME = %f \n", timeNow);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    VecScatterDestroy(&ctx);
    VecDestroy(&vecSEQ);
}

void DirectProblem::updataValiables(const int tItr)
{
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        int n1 = 0;
        for(int i=0; i<grid.node.mapNew[in]; i++)
            n1 += grid.node.nDofsOnNode[i];

        for(int d=0; d<dim; d++)
            grid.node.vPrev[in][d] = 0e0;
        grid.node.p[in] = 0e0;

        for(int d=0; d<dim; d++)
            grid.node.vPrev[in][d] = grid.node.v[in][d];

        for(int d=0; d<dim; d++)
            grid.node.v[in][d] = petsc.solution[n1+d];
        grid.node.p[in] = petsc.solution[n1+dim];
    }

    if(mpi.myId == 0){
        std::string vtuFile;
        vtuFile = outputDir + "/velocity/velocity" + to_string(tItr) + ".vtu";
        grid.outputVTU.exportSolutionVTU(vtuFile, grid.node, grid.cell, DataType::VELOCITY);
        vtuFile = outputDir + "/pressure/pressure" + to_string(tItr) + ".vtu";
        grid.outputVTU.exportSolutionVTU(vtuFile, grid.node, grid.cell, DataType::PRESSURE);
    }
}