/**
 * @file Adjoint.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "InverseProblem.h"

/**
 * @brief Solve adjoint equation.
 * @param main Direct problem.
 * @param outputDir Output directory.
 * @param cb Control boundary.
 */
void Adjoint::solveAdjoint(DirectProblem &main, ControlBoundary &cb)
{
  PetscPrintf(MPI_COMM_WORLD, "\nAdjoit Solver\n");
  PetscScalar *arraySolnTmp;
  Vec vecSEQ;
  VecScatter ctx;

  VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

  petsc.setMatAndVecZero(grid.cell);
  petsc.initialAssembly();

  setVariablesZero(main.dim);
  dirichlet.setValuesZero(grid.nDofsGlobal);

  for(int t = timeMax - 1; t >= 0; t--) {
    petsc.setValueZero();

    dirichlet.assignBCs(grid.node);
    dirichlet.applyBCs(grid.cell, petsc);

    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      if(grid.cell(ic).subId == mpi.myId) {
        int nDofsInCell = grid.cell(ic).dofsMap.size();
        MatrixXd Klocal(nDofsInCell, nDofsInCell);
        VectorXd Flocal(nDofsInCell);
        matrixAssemblyAdjoint(main, Klocal, Flocal, ic, t);
        petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap, grid.cell(ic).dofsBCsMap, Klocal, Flocal);
      }
    }

    for(int ib = 0; ib < cb.CBCellMap.size(); ib++) {
      int ic = cb.CBCellMap[ib];
      if(grid.cell(ic).subId == mpi.myId) {
        int nDofsInCell = grid.cell(ic).dofsMap.size();
        MatrixXd Klocal(nDofsInCell, nDofsInCell);
        VectorXd Flocal(nDofsInCell);
        boundaryIntegral(main, Klocal, Flocal, cb, ic, ib);
        petsc.setMatValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap, Klocal);
      }
    }

    for(int in = 0; in < grid.node.nNodesGlobal; in++) {
      int n = grid.node.mapNew[in];
      int size = grid.node.dofsBCsMapNew[n].size();
      VectorXd Flocal(size);
      Flocal.setZero();
      if(grid.node.subId[in] == mpi.myId) {
        Flocal(0) -= feedbackForceT(t, in, 0);
        Flocal(1) -= feedbackForceT(t, in, 1);
        Flocal(2) -= feedbackForceT(t, in, 2);
        petsc.setVecValue(grid.node.dofsBCsMapNew[n], Flocal);
      }
    }
    petsc.solve();

    VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(vecSEQ, &arraySolnTmp);

    // update solution vector
    for(int id = 0; id < grid.nDofsGlobal; id++)
      petsc.solution[id] = arraySolnTmp[id];

    VecRestoreArray(vecSEQ, &arraySolnTmp);
    updateSolutions();
    updateTimeSolutions(t);

    if(mpi.myId == 0) {
      double timeNow = t * dt;
      printf("Adjoint: Step: %d/%d | Time: %.2fs \n", t, timeMax, timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);
}

/**
 * @brief Set all solutions zero.
 */
void Adjoint::setVariablesZero(const int dim)
{
  w.fillZero();
  q.fillZero();
  wPrev.fillZero();
  qPrev.fillZero();
  l.fillZero();
  wt.fillZero();
  qt.fillZero();
  lt.fillZero();
}

/**
 * @brief Update solutions for VTI.
 */
void Adjoint::updateSolutionsVTI()
{
  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < dim; d++) {
      wvti(grid.vecFluidUniqueNodes[in], d) = w(in, d);
      lvti(grid.vecFluidUniqueNodes[in], d) = l(in, d);
    }
    qvti(grid.vecFluidUniqueNodes[in]) = q(in);
  }
}

/**
 * @brief Update solutions for VTI.
 */
void Adjoint::updateSolutionsVTI(const int t)
{
  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < dim; d++) {
      wvti(grid.vecFluidUniqueNodes[in], d) = wt(t, in, d);
      lvti(grid.vecFluidUniqueNodes[in], d) = lt(t, in, d);
    }
    qvti(grid.vecFluidUniqueNodes[in]) = qt(t, in);
  }
}

/**
 * @brief Update solutions for next time step.
 */
void Adjoint::updateSolutions()
{
  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    int n1 = 0;
    for(int i = 0; i < grid.node.mapNew[in]; i++) {
      n1 += grid.node.nDofsOnNodeNew[i];
    }

    for(int d = 0; d < dim; d++) {
      wPrev(in, d) = w(in, d);
    }
    qPrev(in) = q(in);

    for(int d = 0; d < dim; d++) {
      w(in, d) = petsc.solution[n1 + d];
    }
    q(in) = petsc.solution[n1 + dim];

    if(grid.node.nDofsOnNodeNew[grid.node.mapNew[in]] > dim + 1) {
      for(int d = 0; d < dim; d++)
        l(in, d) = petsc.solution[n1 + dim + 1 + d];
    }
  }
}

/**
 * @brief Save solutions over time.
 */
void Adjoint::updateTimeSolutions(const int t)
{
  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < dim; d++) {
      wt(t, in, d) = w(in, d);
      lt(t, in, d) = l(in, d);
    }
    qt(t, in) = q(in);
  }
}

/**
 * @brief Output solutions for VTU.
 */
void Adjoint::outputSolutionsVTU(const std::string &dir, const int t)
{
  if(mpi.myId > 0)
    return;

  std::string vtuFile;
  vtuFile = outputDir + "/" + dir + "/w_" + to_string(t) + ".vtu";
  EXPORT::exportVectorPointDataVTU(vtuFile, "w", grid.node, grid.cell, w);
  vtuFile = outputDir + "/" + dir + "/q_" + to_string(t) + ".vtu";
  EXPORT::exportScalarPointDataVTU(vtuFile, "q", grid.node, grid.cell, q);
  vtuFile = outputDir + "/" + dir + "/l_" + to_string(t) + ".vtu";
  EXPORT::exportVectorPointDataVTU(vtuFile, "l", grid.node, grid.cell, l);
}

/**
 * @brief Output solutions for VTU.
 */
void Adjoint::outputSolutionsVTU(const std::string &dir, const int t, const int loop)
{
  if(mpi.myId > 0)
    return;

  std::string vtuFile;
  vtuFile = outputDir + "/" + dir + "/w_" + to_string(loop) + "_" + to_string(t) + ".vtu";
  EXPORT::exportVectorPointDataVTU(vtuFile, "w", grid.node, grid.cell, wt, t);
  vtuFile = outputDir + "/" + dir + "/q_" + to_string(loop) + "_" + to_string(t) + ".vtu";
  EXPORT::exportScalarPointDataVTU(vtuFile, "q", grid.node, grid.cell, qt, t);
  vtuFile = outputDir + "/" + dir + "/l_" + to_string(loop) + "_" + to_string(t) + ".vtu";
  EXPORT::exportVectorPointDataVTU(vtuFile, "l", grid.node, grid.cell, lt, t);
}

/**
 * @brief Output solutions for VTI.
 */
void Adjoint::outputSolutionsVTI(const std::string &dir, const int t)
{
  if(mpi.myId > 0)
    return;

  std::string vtiFile;
  vtiFile = outputDir + "/" + dir + "/w_" + to_string(t) + ".vti";
  EXPORT::exportVectorPointDataVTI(vtiFile, "w", wvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
  vtiFile = outputDir + "/" + dir + "/q_" + to_string(t) + ".vti";
  EXPORT::exportScalarPointDataVTI(vtiFile, "q", qvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
  vtiFile = outputDir + "/" + dir + "/l_" + to_string(t) + ".vti";
  EXPORT::exportVectorPointDataVTI(vtiFile, "l", lvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
}

/**
 * @brief Output solutions for VTI.
 */
void Adjoint::outputSolutionsVTI(const std::string &dir, const int t, const int loop)
{
  if(mpi.myId > 0)
    return;

  std::string vtiFile;
  vtiFile = outputDir + "/" + dir + "/w_" + to_string(loop) + "_" + to_string(t) + ".vti";
  EXPORT::exportVectorPointDataVTI(vtiFile, "w", wvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
  vtiFile = outputDir + "/" + dir + "/q_" + to_string(loop) + "_" + to_string(t) + ".vti";
  EXPORT::exportScalarPointDataVTI(vtiFile, "q", qvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
  vtiFile = outputDir + "/" + dir + "/l_" + to_string(loop) + "_" + to_string(t) + ".vti";
  EXPORT::exportVectorPointDataVTI(vtiFile, "l", lvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
}
