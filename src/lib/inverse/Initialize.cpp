/**
 * @file Initialize.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "InverseProblem.h"

/**
 * @brief Initialize inverse problem.
 */
void InverseProblem::initialize(Config &conf)
{
  PetscPrintf(MPI_COMM_WORLD, "\n*** Main initialize ***\n\n");
  // initialize main
  main.dirichlet.initialize(conf);
  main.grid.cell.initialize(conf);
  main.grid.node.initialize(conf);
  main.grid.prepareMatrix(main.dirichlet, main.petsc, main.outputDir, main.timeMax);
  main.dirichlet.getNewArray(main.grid.node.mapNew);

  main.comp_Re(main.dirichlet.velocitySet);

  if(mpi.myId == 0)
    std::cout << "Re: " << main.Re << std::endl;

  PetscPrintf(MPI_COMM_WORLD, "\n*** Adjoint initialize ***\n\n");
  // initialize adjoint
  adjoint.grid.cell.initializeAdjoint(conf);
  adjoint.dirichlet.initialize(conf);
  inletCB.initialize(conf);
  adjoint.dirichlet.eraseControlNodes(adjoint.grid.cell, inletCB);
  adjoint.grid.node.initializeAdjoint(conf, inletCB.CBNodeMap);
  adjoint.grid.prepareMatrix(adjoint.dirichlet, adjoint.petsc, outputDir, adjoint.timeMax);
  adjoint.dirichlet.getNewArray(adjoint.grid.node.mapNew);

  for(int ic = 0; ic < adjoint.grid.cell.nCellsGlobal; ic++) {
    int count = 0;
    VecTool::resize(adjoint.grid.cell(ic).dofStartPlane, 4);
    for(int p = 0; p < adjoint.grid.cell.nNodesInCell; p++) {
      if(adjoint.grid.node.nDofsOnNode[adjoint.grid.cell(ic).node[p]] > dim + 1) {
        adjoint.grid.cell(ic).dofStartPlane[count] = adjoint.grid.cell(ic).dofStart[p];
        count++;
      }
    }
    int tmp = adjoint.grid.cell(ic).dofStartPlane[2];
    adjoint.grid.cell(ic).dofStartPlane[2] = adjoint.grid.cell(ic).dofStartPlane[3];
    adjoint.grid.cell(ic).dofStartPlane[3] = tmp;
  }

  // initialize data
  data.initialize(conf);

  // initialize variables
  resize();
  initializeVarZero();
}

void InverseProblem::resize()
{
  // main params
  main.v0.allocate(main.grid.node.nNodesGlobal, 3);
  main.v.allocate(main.grid.node.nNodesGlobal, 3);
  main.vPrev.allocate(main.grid.node.nNodesGlobal, 3);
  main.p.allocate(main.grid.node.nNodesGlobal);

  main.vt.allocate(main.timeMax, main.grid.node.nNodesGlobal, 3);
  main.pt.allocate(main.timeMax, main.grid.node.nNodesGlobal);

  main.dirichlet.values.allocate(main.grid.nDofsGlobal);
  main.dirichlet.initialValues.allocate(main.grid.nDofsGlobal);

  if(main.grid.gridType == GridType::STRUCTURED) {
    main.vvti.allocate(main.grid.node.nNodesStrGlobal, 3);
    main.pvti.allocate(main.grid.node.nNodesStrGlobal);
  }

  main.velCurrent.allocate(main.grid.cell.nNodesInCell, 3);

  VecTool::resize(main.petsc.solution, main.grid.nDofsGlobal);

  // adjoint params
  adjoint.w.allocate(adjoint.grid.node.nNodesGlobal, dim);
  adjoint.wPrev.allocate(adjoint.grid.node.nNodesGlobal, dim);
  adjoint.q.allocate(adjoint.grid.node.nNodesGlobal);
  adjoint.qPrev.allocate(adjoint.grid.node.nNodesGlobal);
  adjoint.l.allocate(main.grid.nNodesGlobal, dim);

  adjoint.lt.allocate(adjoint.timeMax, adjoint.grid.nNodesGlobal, dim);
  adjoint.wt.allocate(adjoint.timeMax, adjoint.grid.nNodesGlobal, dim);
  adjoint.qt.allocate(adjoint.timeMax, adjoint.grid.nNodesGlobal);
  
  adjoint.dirichlet.values.allocate(adjoint.grid.nDofsGlobal);
  adjoint.dirichlet.initialValues.allocate(adjoint.grid.nDofsGlobal);

  if(adjoint.grid.gridType == GridType::STRUCTURED) {
    adjoint.wvti.allocate(main.grid.node.nNodesStrGlobal, 3);
    adjoint.qvti.allocate(main.grid.node.nNodesStrGlobal);
    adjoint.lvti.allocate(main.grid.node.nNodesStrGlobal, 3);
  }

  adjoint.feedbackForce.allocate(main.snap.nSnapShot, main.grid.node.nNodesGlobal, dim);
  adjoint.feedbackForceT.allocate(main.timeMax, main.grid.node.nNodesGlobal, dim);
  
  // inverse params
  gradX0.allocate(main.grid.node.nNodesGlobal, dim);
  gradX.allocate(main.timeMax, main.grid.node.nNodesGlobal, dim);
  X.allocate(main.timeMax, main.grid.node.nNodesGlobal, dim);
  X0.allocate(main.grid.node.nNodesGlobal, dim);

  if(main.grid.gridType == GridType::STRUCTURED) {
    X0vti.allocate(main.grid.node.nNodesStrGlobal, 3);
    Xvti.allocate(main.timeMax, main.grid.node.nNodesStrGlobal, 3);
  }

  VecTool::resize(adjoint.petsc.solution, adjoint.grid.nDofsGlobal);

  // snap
  main.snap.vSnap.allocate(main.snap.nSnapShot, main.grid.nNodesGlobal, 3);

}

void InverseProblem::initializeVarZero()
{
  // main params
  main.v0.fillZero();
  main.v.fillZero();
  main.vPrev.fillZero();
  main.p.fillZero();

  main.vt.fillZero();
  main.pt.fillZero();

  main.dirichlet.values.fillZero();
  main.dirichlet.initialValues.fillZero();

  if(main.grid.gridType == GridType::STRUCTURED) {
    main.vvti.fillZero();
    main.pvti.fillZero();
  }

  // adjoint params
  adjoint.w.fillZero();
  adjoint.wPrev.fillZero();
  adjoint.q.fillZero();
  adjoint.qPrev.fillZero();
  adjoint.l.fillZero();

  adjoint.lt.fillZero();
  adjoint.wt.fillZero();
  adjoint.qt.fillZero();

  adjoint.dirichlet.values.fillZero();
  adjoint.dirichlet.initialValues.fillZero();

  if(adjoint.grid.gridType == GridType::STRUCTURED) {
    adjoint.vvti.fillZero();
    adjoint.pvti.fillZero();
  }

  adjoint.feedbackForce.fillZero();
  adjoint.feedbackForceT.fillZero();

  // inverse prams
  gradX0.fillZero();
  gradX.fillZero();
  X.fillZero();
  X0.fillZero();
}

