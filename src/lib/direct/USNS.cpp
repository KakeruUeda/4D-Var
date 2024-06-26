#include "DirectProblem.h"

void DirectProblem:: solveUSNS(Application &app)
{ 
    PetscPrintf(MPI_COMM_WORLD, "\nUNSTEADY NAVIER STOKES SOLVER\n");
    
    nu = mu / rho;
    Re = 1e0 / nu;  

    PetscScalar *arraySolnTmp;
    Vec  vecSEQ;
    VecScatter  ctx;
    VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

    petsc.solution.resize(grid.nDofsGlobal);
    grid.dirichlet.dirichletBCsValue.resize(grid.nDofsGlobal, 0e0);
    grid.dirichlet.dirichletBCsValueNew.resize(grid.nDofsGlobal, 0e0);
    grid.dirichlet.dirichletBCsValueInit.resize(grid.nDofsGlobal, 0e0);
    grid.dirichlet.dirichletBCsValueNewInit.resize(grid.nDofsGlobal, 0e0);

    grid.dirichlet.assignDirichletBCs(grid.node, dim);

    petsc.setMatAndVecZero(grid.cell);
    petsc.initialAssembly();

    if(app == Application::FDVAR){
        pulsatileFlow = OFF;
        snap.isSnapShot = ON;
    }

    int snapCount = 0;
    if(pulsatileFlow == ON){
        snap.v.resize(snap.nSnapShot, std::vector<std::vector<double>>
                      (grid.nNodesGlobal, std::vector<double>(dim, 0e0)));
    }

    for(int t=0; t<timeMax; t++){
        petsc.setValueZero();
        
         if(pulsatileFlow == ON)
            if(t > pulseBeginItr)
                grid.dirichlet.assignPulsatileBCs(t, dt, T, grid.nDofsGlobal);

        grid.dirichlet.applyDirichletBCs(grid.cell, petsc);
        
        MPI_Barrier(MPI_COMM_WORLD);
        double timer = MPI_Wtime();
        
        for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                MatrixXd  Klocal(nDofsInCell, nDofsInCell);
                VectorXd  Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();

                if(grid.cell(ic).phi > 0.999)
                    matrixAssemblyUSNS(Klocal, Flocal, ic, t);
                else
                    DarcyMatrixAssemblyUSNS(Klocal, Flocal, ic, t);

                petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap,
                               Klocal, Flocal);
            }
        }
        timer = MPI_Wtime() - timer;

        PetscPrintf(MPI_COMM_WORLD, "\nMatrix assembly = %f seconds\n", timer);
        petsc.currentStatus = ASSEMBLY_OK;
    
        MPI_Barrier(MPI_COMM_WORLD); 
        timer = MPI_Wtime();
        petsc.solve();
        timer = MPI_Wtime() - timer;

        PetscPrintf(MPI_COMM_WORLD, "\nPETSc solver = %f seconds \n", timer);

        VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecGetArray(vecSEQ, &arraySolnTmp);

        // update solution vector
        for(int id=0; id<grid.nDofsGlobal; id++)
            petsc.solution[id] = arraySolnTmp[id];

        VecRestoreArray(vecSEQ, &arraySolnTmp);
        updateValiables(t);

        if(snap.isSnapShot == ON){
            if(t >= snap.snapTimeBeginItr && (snapCount < snap.nSnapShot)){
                if((t - snap.snapTimeBeginItr) % snap.snapInterval == 0){
                    snap.takeSnapShot(grid.node.v, snapCount, grid.node.nNodesGlobal, dim);
                    if(mpi.myId == 0){
                        std::string vtuFile;
                        vtuFile = outputDir + "/data/reference" + to_string(snapCount) + ".vtu";
                        grid.output.exportSnapShotVTU(vtuFile, grid.node, grid.cell, snap, snapCount);
                    }
                    snapCount++;
                }
            }
            if(snapCount >= snap.nSnapShot) 
                break;
        }

        if(mpi.myId == 0){
            double timeNow = t * dt;
            printf("\n TIME = %f \n", timeNow);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    VecScatterDestroy(&ctx);
    VecDestroy(&vecSEQ);
}

void DirectProblem::updateValiables(const int t)
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
        vtuFile = outputDir + "/velocity/velocity" + to_string(t) + ".vtu";
        grid.output.exportSolutionVTU(vtuFile, grid.node, grid.cell, DataType::VELOCITY);
        vtuFile = outputDir + "/pressure/pressure" + to_string(t) + ".vtu";
        grid.output.exportSolutionVTU(vtuFile, grid.node, grid.cell, DataType::PRESSURE);
    }
}

void SnapShot::takeSnapShot(std::vector<std::vector<double>> &_v,
                            const int &snapCount, const int &nNodesGlobal, const int &dim)
{
    for(int in=0; in<nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            v[snapCount][in][d] = _v[in][d];
        }
    }
}