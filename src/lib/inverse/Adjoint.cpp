/**
 * @file Adjoint.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "InverseProblem.h"

/*********************************
 * @brief Solve adjoint equation.
 */
void Adjoint::solveAdjoint(DirectProblem &main, std::string outputDir,
                           std::vector<std::vector<std::vector<double>>> &feedbackForceT)
{
    PetscPrintf(MPI_COMM_WORLD, "\nADJOINT SOLVER\n");
    PetscScalar *arraySolnTmp;
    Vec vecSEQ;
    VecScatter ctx;
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
                Function f3(grid.cell.nNodesInCell, dim);
                MatrixXd Klocal(nDofsInCell, nDofsInCell);
                VectorXd Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();
                matrixAssemblyAdjoint(main, Klocal, Flocal, f3, ic, t);
                petsc.setValue(grid.cell(ic).dofsBCsMapWall, grid.cell(ic).dofsMap,
                               grid.cell(ic).dofsBCsMapWall, Klocal, Flocal);
            }
        }

        for(int ib=0; ib<grid.dirichlet.controlCellMap.size(); ib++){
            int ic = grid.dirichlet.controlCellMap[ib];
            if(grid.cell(ic).subId == mpi.myId){
                int nDofsInCell = grid.cell(ic).dofsMap.size();
                Function f2(grid.dirichlet.nControlNodesInCell, dim-1);
                MatrixXd  Klocal(nDofsInCell, nDofsInCell);
                VectorXd  Flocal(nDofsInCell);
                Klocal.setZero();
                Flocal.setZero();
                boundaryIntegral(main, Klocal, Flocal, f2, ic, ib);
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
        updateSolutions();
        updateTimeSolutions(t);

        if(mpi.myId == 0){
            double timeNow = t * dt;
            printf("Adjoint Solver : Time = %f \n", timeNow);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    VecScatterDestroy(&ctx);
    VecDestroy(&vecSEQ);
}

/*********************************
 * @brief Set all solutions zero.
 */
void Adjoint::setVariablesZero(const int dim)
{
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            grid.node.w[in][d] = 0e0;
            grid.node.wPrev[in][d] = 0e0;
            grid.node.l[in][d] = 0e0;
        }
        grid.node.q[in] = 0e0;
        grid.node.qPrev[in] = 0e0;
    }

    for(int t=0; t<timeMax; t++){
        for(int in=0; in<grid.node.nNodesGlobal; in++){
            for(int d=0; d<dim; d++){
                grid.node.lt[t][in][d] = 0e0;
            }
        }
    }
}

/**********************************
 * @brief Update solutions for VTI.
 */
void Adjoint::updateSolutionsVTI()
{    
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            grid.node.wvti[grid.node.sortNode[in]][d] = grid.node.w[in][d];
            grid.node.lvti[grid.node.sortNode[in]][d] = grid.node.l[in][d];
        }
        grid.node.qvti[grid.node.sortNode[in]] = grid.node.q[in];
    }   
}

/**********************************
 * @brief Update solutions for VTI.
 */
void Adjoint::updateSolutionsVTI(const int t)
{    
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            grid.node.wvti[grid.node.sortNode[in]][d] = grid.node.wt[t][in][d];
            grid.node.lvti[grid.node.sortNode[in]][d] = grid.node.lt[t][in][d];
        }
        grid.node.qvti[grid.node.sortNode[in]] = grid.node.qt[t][in];
    }   
}

/*********************************************
 * @brief Update solutions for next time step.
 */
void Adjoint::updateSolutions()
{
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        int n1 = 0;
        for(int i=0; i<grid.node.mapNew[in]; i++)
            n1 += grid.node.nDofsOnNodeNew[i];

        for(int d=0; d<dim; d++)
            grid.node.wPrev[in][d] = grid.node.w[in][d];
        grid.node.qPrev[in] = grid.node.q[in];

        for(int d=0; d<dim; d++)
            grid.node.w[in][d] = petsc.solution[n1+d];
        grid.node.q[in] = petsc.solution[n1+dim];

        if(grid.node.nDofsOnNodeNew[grid.node.mapNew[in]] > dim+1){
            for(int d=0; d<dim; d++)
                grid.node.l[in][d] = petsc.solution[n1+dim+1+d];
        }
    }
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        if(grid.dirichlet.isBoundaryEdge[in]){
            for(int d=0; d<dim; d++){
                grid.node.l[in][d] = 0e0;
            }
        }   
    }
}

/**********************************
 * @brief Save solutions over time.
 */
void Adjoint::updateTimeSolutions(const int t)
{
    for(int in=0; in<grid.node.nNodesGlobal; in++){
        for(int d=0; d<dim; d++){
            grid.node.wt[t][in][d] = grid.node.w[in][d];
            grid.node.lt[t][in][d] = grid.node.l[in][d];
        }
        grid.node.qt[t][in] = grid.node.q[in];
    }
}

/***********************************
 * @brief Output solutions for VTU.
 */
void Adjoint::outputSolutionsVTU(const std::string &dir,const int t)
{
    if(mpi.myId > 0) return;

    std::string vtuFile;
    vtuFile = outputDir + "/" + dir + "/w_" + to_string(t) + ".vtu";
    VTK::exportVectorPointDataVTU(vtuFile, "w", grid.node, grid.cell, grid.node.w);
    vtuFile = outputDir + "/" + dir + "/q_" + to_string(t) + ".vtu";
    VTK::exportScalarPointDataVTU(vtuFile, "q", grid.node, grid.cell, grid.node.q);
    vtuFile = outputDir + "/" + dir + "/l_" + to_string(t) + ".vtu";
    VTK::exportVectorPointDataVTU(vtuFile, "l", grid.node, grid.cell, grid.node.l);
}

/***********************************
 * @brief Output solutions for VTU.
 */
void Adjoint::outputSolutionsVTU(const std::string &dir, const int t, const int loop)
{
    if(mpi.myId > 0) return;

    std::string vtuFile;
    vtuFile = outputDir + "/" + dir + "/w_" + to_string(loop) + "_" + to_string(t) + ".vtu";
    VTK::exportVectorPointDataVTU(vtuFile, "w", grid.node, grid.cell, grid.node.w);
    vtuFile = outputDir + "/" + dir + "/q_" + to_string(loop) + "_" + to_string(t) + ".vtu";
    VTK::exportScalarPointDataVTU(vtuFile, "q", grid.node, grid.cell, grid.node.q);
    vtuFile = outputDir + "/" + dir + "/l_" + to_string(loop) + "_" + to_string(t) + ".vtu";
    VTK::exportVectorPointDataVTU(vtuFile, "l", grid.node, grid.cell, grid.node.l);
}

/***********************************
 * @brief Output solutions for VTI.
 */
void Adjoint::outputSolutionsVTI(const std::string &dir, const int t)
{
    if(mpi.myId > 0) return;

    std::string vtiFile;
    vtiFile = outputDir + "/" + dir + "/w_" + to_string(t) + ".vti";
    VTK::exportVectorPointDataVTI(vtiFile, "w", grid.node.wvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
    vtiFile = outputDir + "/" + dir + "/q_" + to_string(t) + ".vti";
    VTK::exportScalarPointDataVTI(vtiFile, "q", grid.node.qvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
    vtiFile = outputDir + "/" + dir + "/l_" + to_string(t) + ".vti";
    VTK::exportVectorPointDataVTI(vtiFile, "l", grid.node.lvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
}

/***********************************
 * @brief Output solutions for VTI.
 */
void Adjoint::outputSolutionsVTI(const std::string &dir, const int t, const int loop)
{
    if(mpi.myId > 0) return;

    std::string vtiFile;
    vtiFile = outputDir + "/" + dir + "/w_" + to_string(loop) + "_" + to_string(t) + ".vti";
    VTK::exportVectorPointDataVTI(vtiFile, "w", grid.node.wvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
    vtiFile = outputDir + "/" + dir + "/q_" + to_string(loop) + "_" + to_string(t) + ".vti";
    VTK::exportScalarPointDataVTI(vtiFile, "q", grid.node.qvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
    vtiFile = outputDir + "/" + dir + "/l_" + to_string(loop) + "_" + to_string(t) + ".vti";
    VTK::exportVectorPointDataVTI(vtiFile, "l", grid.node.lvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
}

/**************************
 * @brief Update row index.
 */
void Adjoint::updateRowIndex(const int ii, const int ic)
{
    IU =grid.cell(ic).dofStart[ii]; IV = IU + 1;  IW = IU + 2;
    IP = IU + 3; ILU = IU + 4; ILV = IU + 5; ILW = IU + 6;
}

/*****************************
 * @brief Update column index.
 */
void Adjoint::updateColumnIndex(const int jj, const int ic)
{
    JU = grid.cell(ic).dofStart[jj]; JV = JU + 1;  JW = JU + 2;
    JP = JU + 3; JLU = JU + 4; JLV = JU + 5; JLW = JU + 6;
}

/*****************************
 * @brief Update 2D row index.
 */
void Adjoint::updateRowIndexPlane(const int ii, const int ic)
{
    IU = grid.cell(ic).dofStartPlane[ii]; IV = IU + 1; IW = IU + 2;
    IP = IU + 3; ILU = IU + 4; ILV = IU + 5; ILW = IU + 6;
}

/********************************
 * @brief Update 2D column index.
 */
void Adjoint::updateColumnIndexPlane(const int jj, const int ic)
{
    JU = grid.cell(ic).dofStartPlane[jj]; JV = JU + 1; JW = JU + 2;
    JP = JU + 3; JLU = JU + 4; JLV = JU + 5; JLW = JU + 6;
}
