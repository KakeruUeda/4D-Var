#include "DirectProblem.h"

void DirectProblem::initialize(Config &conf)
{
    nu = mu / rho;
    Re = 1e0 / nu; 

    VecTool::resize(grid.node.v, grid.node.nNodesGlobal, dim);
    VecTool::resize(grid.node.vPrev, grid.node.nNodesGlobal, dim);
    //VecTool::resize(grid.node.vt, timeMax, grid.node.nNodesGlobal, dim);
    VecTool::resize(grid.node.p, grid.node.nNodesGlobal);
    //VecTool::resize(grid.node.pt, timeMax, grid.node.nNodesGlobal);
    VecTool::resize(snap.v, snap.nSnapShot, grid.nNodesGlobal, dim);

    grid.dirichlet.initialize(conf);
    grid.cell.initialize(conf);
    grid.node.initialize(conf);
   
    grid.prepareMatrix(petsc, outputDir, timeMax);

    VecTool::resize(petsc.solution, grid.nDofsGlobal);
    VecTool::resize(grid.dirichlet.dirichletBCsValue, grid.nDofsGlobal);
    VecTool::resize(grid.dirichlet.dirichletBCsValueNew, grid.nDofsGlobal);
    VecTool::resize(grid.dirichlet.dirichletBCsValueInit, grid.nDofsGlobal);
    VecTool::resize(grid.dirichlet.dirichletBCsValueNewInit, grid.nDofsGlobal);
    
    initializeFEM();
}

void DirectProblem::initializeFEM()
{
    detJ = 0e0;
    weight = 0e0;

    VecTool::resize(N, grid.cell.nNodesInCell);
    VecTool::resize(xCurrent, grid.cell.nNodesInCell, dim);
    VecTool::resize(dNdr, grid.cell.nNodesInCell, dim);
    VecTool::resize(dNdx, grid.cell.nNodesInCell, dim);
    VecTool::resize(K, grid.cell.nNodesInCell, grid.cell.nNodesInCell);

    VecTool::resize(dofStart, grid.cell.nCellsGlobal, grid.cell.nNodesInCell);

    VecTool::resize(vgp, dim);
    VecTool::resize(advgp, dim);
    VecTool::resize(dvgpdx, dim, dim);

    for(int ic=0; ic<grid.cell.nCellsGlobal; ic++){
        for(int p=1; p<grid.cell.nNodesInCell; p++){
            dofStart[ic][p] = dofStart[ic][p-1] + grid.node.nDofsOnNode[grid.cell(ic).node[p-1]]; 
        }
    }
}