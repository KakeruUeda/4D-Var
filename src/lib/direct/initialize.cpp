/**
 * @file Initialize.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "DirectProblem.h"

/***********************************
 * @brief Initialize direct problem.
 */
void DirectProblem::initialize(Config &conf)
{
  nu = mu / rho;
  Re = 1e0 / nu;

  grid.cell.initialize(conf);
  grid.node.initialize(conf);
  dirichlet.initialize(conf);

  grid.prepareMatrix(dirichlet, petsc, outputDir, timeMax);
  dirichlet.getNewArray(grid.node.mapNew);
}

void DirectProblem::resize()
{
  VecTool::resize(grid.node.v, grid.node.nNodesGlobal, dim);
  VecTool::resize(grid.node.vPrev, grid.node.nNodesGlobal, dim);
  VecTool::resize(grid.node.p, grid.node.nNodesGlobal);

  v.resize(grid.node.nNodesGlobal, dim);
  vPrev.resize(grid.node.nNodesGlobal, dim);
  p.resize(grid.node.nNodesGlobal);

	vrt.resize(grid.node.nStrNodesGlobal, dim);

  dirichlet.values.resize(grid.nDofsGlobal);
  dirichlet.initialValues.resize(grid.nDofsGlobal);

  v.fillZero();
  vPrev.fillZero();
  p.fillZero();
	vrt.fillZero();

  dirichlet.values.fillZero();
  dirichlet.initialValues.fillZero();

  if(grid.gridType == GridType::STRUCTURED) {
    vvti.resize(grid.node.nNodesStrGlobal, dim);
    pvti.resize(grid.node.nNodesStrGlobal);
  }

  VecTool::resize(petsc.solution, grid.nDofsGlobal);
  VecTool::resize(grid.dirichlet.dirichletBCsValue, grid.nDofsGlobal);
  VecTool::resize(grid.dirichlet.dirichletBCsValueNew, grid.nDofsGlobal);
  VecTool::resize(grid.dirichlet.dirichletBCsValueInit, grid.nDofsGlobal);
  VecTool::resize(grid.dirichlet.dirichletBCsValueNewInit, grid.nDofsGlobal);

  VecTool::resize(vgp, dim);
  VecTool::resize(advgp, dim);
  VecTool::resize(dvgpdx, dim, dim);
}
