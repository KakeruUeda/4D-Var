/**
 * @file Initialize.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "InverseProblem.h"

/************************************
 * @brief Initialize inverse problem.
 */
void InverseProblem::initialize(Config &conf)
{
    PetscPrintf(MPI_COMM_WORLD, "\n*** Main initialize ***\n\n");

    main.nu = main.mu / main.rho;
    main.Re = 1e0 / main.nu; 

    VecTool::resize(main.grid.node.v0, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(main.grid.node.v, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(main.grid.node.vPrev, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(main.grid.node.vt, main.timeMax, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(main.grid.node.p, main.grid.node.nNodesGlobal);
    VecTool::resize(main.grid.node.pt, main.timeMax, main.grid.node.nNodesGlobal);
    VecTool::resize(main.snap.v, main.snap.nSnapShot, main.grid.nNodesGlobal, dim);

    VecTool::resize(main.vgp, dim);
    VecTool::resize(main.advgp, dim);
    VecTool::resize(main.dvgpdx, dim, dim);

    main.grid.dirichlet.initialize(conf);
    main.grid.cell.initialize(conf);
    main.grid.node.initialize(conf);
    main.grid.prepareMatrix(main.petsc, main.outputDir, main.timeMax);

    VecTool::resize(main.petsc.solution, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValue, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValueNew, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValueInit, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValueNewInit, main.grid.nDofsGlobal);

    PetscPrintf(MPI_COMM_WORLD, "\n*** Adjoint initialize ***\n\n");

    adjoint.nu = adjoint.mu / adjoint.rho;
    adjoint.Re = 1e0 / adjoint.nu; 

    isConverged_X = false;
    isConverged_X0 = false;

    VecTool::resize(adjoint.grid.node.w, adjoint.grid.node.nNodesGlobal, dim);
    VecTool::resize(adjoint.grid.node.wPrev, adjoint.grid.node.nNodesGlobal, dim);
    VecTool::resize(adjoint.grid.node.q, adjoint.grid.node.nNodesGlobal);
    VecTool::resize(adjoint.grid.node.l, main.grid.nNodesGlobal, dim);
    VecTool::resize(adjoint.grid.node.lt, adjoint.timeMax, adjoint.grid.nNodesGlobal, dim);
    VecTool::resize(adjoint.grid.node.wt, adjoint.timeMax, adjoint.grid.nNodesGlobal, dim);
    VecTool::resize(adjoint.grid.node.qt, adjoint.timeMax, adjoint.grid.nNodesGlobal);

    VecTool::resize(adjoint.vk, dim);
    VecTool::resize(adjoint.vk1, dim);
    VecTool::resize(adjoint.vk2, dim);
    VecTool::resize(adjoint.advk1, dim);
    VecTool::resize(adjoint.advk2, dim);
    VecTool::resize(adjoint.dvkdx, dim, dim);
    VecTool::resize(adjoint.dvk1dx, dim, dim);
    VecTool::resize(adjoint.dvk2dx, dim, dim);
    VecTool::resize(adjoint.wk, dim);
    VecTool::resize(adjoint.wk1, dim);
    VecTool::resize(adjoint.wk2, dim);
    VecTool::resize(adjoint.dwkdx, dim, dim);
    VecTool::resize(adjoint.dwk1dx, dim, dim);
    VecTool::resize(adjoint.dwk2dx, dim, dim);
    
    adjoint.grid.dirichlet.initializeAdjoint(conf);
    adjoint.grid.cell.initializeAdjoint(conf);
    adjoint.grid.node.initializeAdjoint(conf, adjoint.grid.dirichlet.controlBoundaryMap);
    adjoint.grid.prepareMatrix(adjoint.petsc, outputDir, adjoint.timeMax);
    
    for(int ic=0; ic<adjoint.grid.cell.nCellsGlobal; ic++){
        int count = 0;
        VecTool::resize(adjoint.grid.cell(ic).dofStartPlane, adjoint.grid.dirichlet.nControlNodesInCell);
        for(int p=0; p<adjoint.grid.cell.nNodesInCell; p++){
            if(adjoint.grid.node.nDofsOnNode[adjoint.grid.cell(ic).node[p]] > dim+1){
                adjoint.grid.cell(ic).dofStartPlane[count] = adjoint.grid.cell(ic).dofStart[p];
                count++;
            }
        }
        int tmp = adjoint.grid.cell(ic).dofStartPlane[2];
        adjoint.grid.cell(ic).dofStartPlane[2] = adjoint.grid.cell(ic).dofStartPlane[3];
        adjoint.grid.cell(ic).dofStartPlane[3] = tmp;
    }   

    VecTool::resize(adjoint.petsc.solution, adjoint.grid.nDofsGlobal);
    VecTool::resize(adjoint.grid.dirichlet.dirichletBCsValue, adjoint.grid.nDofsGlobal);
    VecTool::resize(adjoint.grid.dirichlet.dirichletBCsValueNew, adjoint.grid.nDofsGlobal);
    VecTool::resize(adjoint.grid.dirichlet.dirichletBCsValueInit, adjoint.grid.nDofsGlobal);
    VecTool::resize(adjoint.grid.dirichlet.dirichletBCsValueNewInit, adjoint.grid.nDofsGlobal);
    VecTool::resize(feedbackForce, main.snap.nSnapShot, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(feedbackForceT, main.timeMax, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(gradWholeNode, main.timeMax, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(grad, main.timeMax, adjoint.grid.dirichlet.controlBoundaryMap.size(), dim);
    VecTool::resize(X, main.timeMax, adjoint.grid.dirichlet.controlBoundaryMap.size(), dim);
    VecTool::resize(X0, main.grid.node.nNodesGlobal, dim);
    VecTool::resize(gradInitVel, main.grid.node.nNodesGlobal, dim);

    VecTool::resize(adjoint.vgp, dim);
    VecTool::resize(adjoint.advgp, dim);
    VecTool::resize(adjoint.dvgpdx, dim, dim);

    data.initialize(conf, main.grid.node, main.grid.cell, dim);
}
