/**
 * @file Initialize.cpp
 * @author K.U.
 * @date July, 2024
 */

#include "DirectProblem.h"

void DirectProblem::initialize(Config &conf)
{
    nu = mu / rho;
    Re = 1e0 / nu; 
    
    VecTool::resize(grid.node.v, grid.node.nNodesGlobal, dim);
    VecTool::resize(grid.node.vPrev, grid.node.nNodesGlobal, dim);
    VecTool::resize(grid.node.p, grid.node.nNodesGlobal);
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

    VecTool::resize(vgp, dim);
    VecTool::resize(advgp, dim);
    VecTool::resize(dvgpdx, dim, dim);
}
