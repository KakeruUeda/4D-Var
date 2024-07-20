/**
 * @file USNS.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "DirectProblem.h"

/**************************************
 * @brief Solve Unsteady Navier Stokes.
 */
void DirectProblem::solveUSNS(Application &app)
{ 
    PetscPrintf(MPI_COMM_WORLD, "\nMain solver\n");

    PetscScalar *arraySolnTmp;
    Vec  vecSEQ;
    VecScatter  ctx;
    VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

    petsc.setMatAndVecZero(grid.cell);
    petsc.initialAssembly();

    for(int id=0; id<grid.nDofsGlobal; id++){
        grid.dirichlet.dirichletBCsValueNewInit[id] = 0e0;
        grid.dirichlet.dirichletBCsValueNew[id] = 0e0;
    }
    if(app == Application::FDVAR){
        pulsatileFlow = OFF;
        snap.isSnapShot = ON;
    }

    int snapCount = 0;
    for(int t=0; t<timeMax; t++){
        petsc.setValueZero();
        grid.dirichlet.assignDirichletBCs(grid.dirichlet.vDirichletNew, 
                                          grid.dirichlet.pDirichletNew, 
                                          grid.node, dim, t);
         if(pulsatileFlow == ON)
            if(t > pulseBeginItr)
                grid.dirichlet.assignPulsatileBCs(t, dt, T, grid.nDofsGlobal);

        grid.dirichlet.applyDirichletBCs(grid.cell, petsc);
        
        MPI_Barrier(MPI_COMM_WORLD);
        double timer1 = MPI_Wtime();
        
        for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                Function f3(grid.cell.nNodesInCell, dim);
                MatrixXd Klocal(nDofsInCell, nDofsInCell);
                VectorXd Flocal(nDofsInCell);
                Klocal.setZero(); 
                Flocal.setZero();
                matrixAssemblyUSNS(Klocal, Flocal, f3, ic, t);
                petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap,
                               grid.cell(ic).dofsBCsMap, Klocal, Flocal);
            }
        }
        petsc.currentStatus = ASSEMBLY_OK;
        timer1 = MPI_Wtime() - timer1;
        //PetscPrintf(MPI_COMM_WORLD, "\nMatrix assembly = %f seconds\n", timer);
        MPI_Barrier(MPI_COMM_WORLD); 

        double timer2 = MPI_Wtime();
        petsc.solve();
        timer2 = MPI_Wtime() - timer2;

        //PetscPrintf(MPI_COMM_WORLD, "PETSc solver = %f seconds \n", timer);
        VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
        VecGetArray(vecSEQ, &arraySolnTmp);

        // update solution vector
        for(int id=0; id<grid.nDofsGlobal; id++)
            petsc.solution[id] = arraySolnTmp[id];

        VecRestoreArray(vecSEQ, &arraySolnTmp);
        updateVariables(t);
        outputSolutionVTI(t);

        if(app == Application::FDVAR)
            assignTimeVariables(t);

        if(snap.isSnapShot == ON){
            if(t >= snap.snapTimeBeginItr && (snapCount < snap.nSnapShot)){
                if((t - snap.snapTimeBeginItr) % snap.snapInterval == 0){
                    snap.takeSnapShot(grid.node.v, snapCount, grid.node.nNodesGlobal, dim);
                    snapCount++;
                }
            }
            if(snapCount >= snap.nSnapShot) 
                break;
        }

        if(mpi.myId == 0){
            double timeNow = t * dt;
            printf("Assy: %fs | Solve: %fs | SimTime: %fs \n", timer1, timer2, timeNow);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    VecScatterDestroy(&ctx);
    VecDestroy(&vecSEQ);
}

/*********************************************************
 * @brief Solve Unsteady Navier Stokes. 
 *        This function takes boundary arguments and 
 *        mainly created for checking armijo criteria.
 */
void DirectProblem::solveUSNS(std::vector<std::map<int, std::vector<double>>> &vDirichletTmp,
                              std::vector<std::map<int, double>> &pDirichletTmp)
{ 
    PetscPrintf(MPI_COMM_WORLD, "\nMain Solver\n");

    PetscScalar *arraySolnTmp;
    Vec  vecSEQ;
    VecScatter  ctx;
    VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ); 
    
    petsc.setMatAndVecZero(grid.cell);
    petsc.initialAssembly();

    for(int id=0; id<grid.nDofsGlobal; id++){
        grid.dirichlet.dirichletBCsValueNewInit[id] = 0e0;
        grid.dirichlet.dirichletBCsValueNew[id] = 0e0;
    }   

    int snapCount = 0;
    for(int t=0; t<timeMax; t++){
        if(t == 0){
            assignTimeVariables(t);
            if((t - snap.snapTimeBeginItr) % snap.snapInterval == 0){;
                snap.takeSnapShot(grid.node.v, snapCount, grid.node.nNodesGlobal, dim);
                snapCount++;
            }
            if(mpi.myId == 0){
                double timeNow = t * dt;
                printf("Main Solver : Time = %f \n", timeNow);
            }
            continue;
        }
        petsc.setValueZero();
        grid.dirichlet.assignDirichletBCs(vDirichletTmp, pDirichletTmp, 
                                          grid.node, dim, t);
        if(pulsatileFlow == ON)
            if(t > pulseBeginItr)
                grid.dirichlet.assignPulsatileBCs(t, dt, T, grid.nDofsGlobal);

        grid.dirichlet.applyDirichletBCs(grid.cell, petsc); 
        
        for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                Function f3(grid.cell.nNodesInCell, dim);
                MatrixXd Klocal(nDofsInCell, nDofsInCell);
                VectorXd Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();
                matrixAssemblyUSNS(Klocal, Flocal, f3, ic, t);
                petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap,
                               grid.cell(ic).dofsBCsMap, Klocal, Flocal);
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
        updateVariables(t);
        assignTimeVariables(t);
        
        if((t - snap.snapTimeBeginItr) % snap.snapInterval == 0){;
            snap.takeSnapShot(grid.node.v, snapCount, grid.node.nNodesGlobal, dim);
            snapCount++;
        }

        if(mpi.myId == 0){
            double timeNow = t * dt;
            printf("Main Solver : Time = %f \n", timeNow);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    VecScatterDestroy(&ctx);
    VecDestroy(&vecSEQ);
}

/*******************************************************************
 * @brief Compute fully developed flow field.
 *        Usually used to get initial condition for inverse problem.
 */
void DirectProblem::compInitialCondition(std::vector<std::map<int, std::vector<double>>> &vDirichletTmp,
                                         std::vector<std::map<int, double>> &pDirichletTmp)
{
    PetscPrintf(MPI_COMM_WORLD, "\nCompute Initial Condition\n");
    
    double norm, norm0;
    PetscScalar *arraySolnTmp;
    Vec  vecSEQ;
    VecScatter  ctx;
    VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ); 

    petsc.setMatAndVecZero(grid.cell);
    petsc.initialAssembly();
    setVariablesZero();
    
    for(int id=0; id<grid.nDofsGlobal; id++){
        grid.dirichlet.dirichletBCsValueNewInit[id] = 0e0;
        grid.dirichlet.dirichletBCsValueNew[id] = 0e0;
    }   

    int snapCount = 0;
    for(int t=0; t<timeMax; t++){
        petsc.setValueZero();
        grid.dirichlet.assignConstantDirichletBCs(vDirichletTmp, pDirichletTmp, grid.node, dim, t);
        grid.dirichlet.applyDirichletBCs(grid.cell, petsc);
        
        for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                Function f3(grid.cell.nNodesInCell, dim);
                MatrixXd Klocal(nDofsInCell, nDofsInCell);
                VectorXd Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();
                matrixAssemblyUSNS(Klocal, Flocal, f3, ic, t);
                petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap,
                               grid.cell(ic).dofsBCsMap, Klocal, Flocal);
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
        updateVariables(t);

        if(t == timeMax-1)
            for(int in=0; in<grid.node.nNodesGlobal; in++)
                for(int d=0; d<dim; d++)
                    grid.node.v0[in][d] = grid.node.v[in][d];

        if(mpi.myId == 0){
            double timeNow = t * dt;
            printf("Compute initial condition : Time = %f \n", timeNow);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    VecScatterDestroy(&ctx);
    VecDestroy(&vecSEQ);
}

/********************************
 * @brief Set all solutions zero.
 */
void DirectProblem::setVariablesZero()
{
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            grid.node.v[in][d] = 0e0;
            grid.node.vPrev[in][d] = 0e0;
        }
        grid.node.p[in] = 0e0;
    }

    for(int t=0; t<timeMax; t++){
        for(int in=0; in<grid.node.nNodesGlobal; in++){
            for(int d=0; d<dim; d++){
                grid.node.vt[t][in][d] = 0e0;
            }
            grid.node.pt[t][in] = 0e0;
        }
    }
}

/******************************
 * @brief Update all solutions.
 */
void DirectProblem::updateVariables(const int t)
{
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        int n1 = 0;
        for(int i=0; i<grid.node.mapNew[in]; i++)
            n1 += grid.node.nDofsOnNode[i];

        for(int d=0; d<dim; d++)
            grid.node.vPrev[in][d] = grid.node.v[in][d];

        for(int d=0; d<dim; d++)
            grid.node.v[in][d] = petsc.solution[n1+d];
        grid.node.p[in] = petsc.solution[n1+dim];
    }
}

/*****************************
 * @brief Visualize solutions.
 */
void DirectProblem::outputSolution(const int t)
{
    if(mpi.myId > 0) return;
    std::string vtuFile;
    vtuFile = outputDir + "/solution/velocity" + to_string(t) + ".vtu";
    grid.vtk.exportSolutionVTU(vtuFile, grid.node, grid.cell, DataType::VELOCITY);
    vtuFile = outputDir + "/solution/pressure" + to_string(t) + ".vtu";
    grid.vtk.exportSolutionVTU(vtuFile, grid.node, grid.cell, DataType::PRESSURE);
}

void DirectProblem::outputSolutionVTI(const int t)
{
    if(mpi.myId > 0) return;

    std::vector<std::vector<double>> vvti;
    std::vector<double> pvti;

    int nNodesGlobalPrev = (grid.nx+1) * (grid.ny+1) * (grid.nz+1);

    VecTool::resize(vvti, nNodesGlobalPrev, dim);
    VecTool::resize(pvti, nNodesGlobalPrev);
    
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            vvti[grid.node.sortNode[in]][d] = grid.node.v[in][d];
        }
        pvti[grid.node.sortNode[in]] = grid.node.p[in];
    }

    std::string vtiFile;
    vtiFile = outputDir + "/solution/solution" + to_string(t) + ".vti";
    grid.vtk.exportSolutionVTI(vtiFile, vvti, pvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
}

/**************************************************************
 * @brief Assign time-dependent variables for adjoint equation.
 */
void DirectProblem::assignTimeVariables(const int t)
{
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            grid.node.vt[t][in][d] = grid.node.v[in][d];
        }
        grid.node.pt[t][in] = grid.node.p[in];
    }
}

/********************************************
 * @brief Take snapshots for error functions.
 */
void SnapShot::takeSnapShot(std::vector<std::vector<double>> &_v,
                            const int &snapCount, const int &nNodesGlobal, const int &dim)
{
    for(int in=0; in<nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            v[snapCount][in][d] = _v[in][d];
        }
    }
}


