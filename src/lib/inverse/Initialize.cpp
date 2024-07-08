#include "InverseProblem.h"

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

    main.grid.dirichlet.initialize(conf);
    main.grid.cell.initialize(conf);
    main.grid.node.initialize(conf);
    main.grid.prepareMatrix(main.petsc, main.outputDir, main.timeMax);

    VecTool::resize(main.petsc.solution, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValue, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValueNew, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValueInit, main.grid.nDofsGlobal);
    VecTool::resize(main.grid.dirichlet.dirichletBCsValueNewInit, main.grid.nDofsGlobal);

    main.initializeFEM();

    PetscPrintf(MPI_COMM_WORLD, "\n*** Adjoint initialize ***\n\n");

    adjoint.nu = adjoint.mu / adjoint.rho;
    adjoint.Re = 1e0 / adjoint.nu; 

    VecTool::resize(adjoint.grid.node.w, adjoint.grid.node.nNodesGlobal, dim);
    VecTool::resize(adjoint.grid.node.wPrev, adjoint.grid.node.nNodesGlobal, dim);
    VecTool::resize(adjoint.grid.node.q, adjoint.grid.node.nNodesGlobal);
    VecTool::resize(adjoint.grid.node.l, main.grid.nNodesGlobal, dim);
    VecTool::resize(adjoint.grid.node.lt, adjoint.timeMax, adjoint.grid.nNodesGlobal, dim);

    VecTool::resize(adjoint.vk, dim);
    VecTool::resize(adjoint.vk1, dim);
    VecTool::resize(adjoint.vk2, dim);
    VecTool::resize(adjoint.advk1, dim);
    VecTool::resize(adjoint.advk2, dim);
    VecTool::resize(adjoint.dvkdx, dim, dim);
    VecTool::resize(adjoint.dvkdx, dim, dim);
    VecTool::resize(adjoint.dvk2dx, dim, dim);
    VecTool::resize(adjoint.dwk1dx, dim, dim);
    VecTool::resize(adjoint.dwk2dx, dim, dim);

    adjoint.grid.dirichlet.initializeAdjoint(conf);
    adjoint.grid.cell.initializeAdjoint(conf);
    adjoint.grid.node.initializeAdjoint(conf, adjoint.grid.dirichlet.controlBoundaryMap);
    adjoint.grid.prepareMatrix(adjoint.petsc, outputDir, adjoint.timeMax);

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

    adjoint.initializeFEM();

    data.initialize(conf, main.grid.node, main.grid.cell, dim);
}

void Adjoint::initializeFEM()
{
    detJ = 0e0;
    weight = 0e0;

    VecTool::resize(N, grid.cell.nNodesInCell);
    VecTool::resize(xCurrent, grid.cell.nNodesInCell, dim);
    VecTool::resize(dNdr, grid.cell.nNodesInCell, dim);
    VecTool::resize(dNdx, grid.cell.nNodesInCell, dim);
    VecTool::resize(K, grid.cell.nNodesInCell, grid.cell.nNodesInCell);

    VecTool::resize(N2D, grid.dirichlet.nControlNodesGlobal);
    VecTool::resize(xCurrent2D, grid.dirichlet.nControlNodesGlobal, dim-1);
    VecTool::resize(dNdr2D, grid.dirichlet.nControlNodesGlobal, dim-1);
    VecTool::resize(dNdx2D, grid.dirichlet.nControlNodesGlobal, dim-1);

    VecTool::resize(dofStart, grid.cell.nCellsGlobal, grid.cell.nNodesInCell);
    VecTool::resize(dofStartPlane, grid.cell.nCellsGlobal, grid.dirichlet.nControlNodesGlobal);

    VecTool::resize(vgp, dim);
    VecTool::resize(advgp, dim);
    VecTool::resize(dvgpdx, dim, dim);

    for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
        for(int p=1; p<grid.cell.nNodesInCell; p++){
            dofStart[ic][p] = dofStart[ic][p-1] + grid.node.nDofsOnNode[grid.cell(ic).node[p-1]]; 
        }
        int count = 0;
        for(int p=0; p<grid.cell.nNodesInCell; p++){
            if(grid.node.nDofsOnNode[grid.cell(ic).node[p]] > dim+1){
                dofStartPlane[ic][count] = dofStart[ic][p];
                count++;
            }
        }
        int tmp = dofStartPlane[ic][2];
        dofStartPlane[ic][2] = dofStartPlane[ic][3];
        dofStartPlane[ic][3] = tmp;
    }

}