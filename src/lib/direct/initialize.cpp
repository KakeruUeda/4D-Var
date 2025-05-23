/**
 * @file Initialize.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "DirectProblem.h"

/**
 * @brief Initialize direct problem.
 */
void DirectProblem::initialize(Config &conf)
{
  grid.cell.initialize(conf);
  grid.node.initialize(conf);
  dirichlet.initialize(conf);

  comp_Re(dirichlet.velocitySet);

  if(mpi.myId == 0){
    std::cout << "Re = " << Re << std::endl;
  }

  grid.prepareMatrix(dirichlet, petsc, outputDir, timeMax);
  dirichlet.getNewArray(grid.node.mapNew);

	resizeVar();
	initializeVarZero();
}

void DirectProblem::resizeVar()
{
  VecTool::resize(grid.node.v, grid.node.nNodesGlobal, dim);
  VecTool::resize(grid.node.vPrev, grid.node.nNodesGlobal, dim);
  VecTool::resize(grid.node.p, grid.node.nNodesGlobal);

  v.allocate(grid.node.nNodesGlobal, 3);
  vPrev.allocate(grid.node.nNodesGlobal, 3);
  p.allocate(grid.node.nNodesGlobal);

	vrt.allocate(grid.node.nStrNodesGlobal, dim);

  dirichlet.values.allocate(grid.nDofsGlobal);
  dirichlet.initialValues.allocate(grid.nDofsGlobal);

  if(grid.gridType == GridType::STRUCTURED) {
    vvti.allocate(grid.node.nNodesStrGlobal, dim);
    pvti.allocate(grid.node.nNodesStrGlobal);
  }

  velCurrent.allocate(grid.cell.nNodesInCell, 3);

  VecTool::resize(petsc.solution, grid.nDofsGlobal);
}

void DirectProblem::initializeVarZero()
{
  v.fillZero();
  vPrev.fillZero();
  p.fillZero();
	vrt.fillZero();

	dirichlet.values.fillZero();
  dirichlet.initialValues.fillZero();

	if(grid.gridType == GridType::STRUCTURED) {
    vvti.fillZero();
    pvti.fillZero();
  }
}