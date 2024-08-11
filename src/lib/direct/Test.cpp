/**
 * @file USNS.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "DirectProblem.h"

/**************************************
 * @brief Solve Unsteady Navier Stokes.
 */
void DirectProblem::solveNavierStokes()
{
  PetscPrintf(MPI_COMM_WORLD, "\nMain solver\n");

  PetscScalar *arraySolnTmp;
  Vec vecSEQ;
  VecScatter ctx;
  VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

  petsc.setMatAndVecZero(grid.cell);
  petsc.initialAssembly();

  for (int id = 0; id < grid.nDofsGlobal; id++)
  {
    grid.dirichlet.dirichletBCsValueNewInit[id] = 0e0;
    grid.dirichlet.dirichletBCsValueNew[id] = 0e0;
  }

  int snapCount = 0;
  for (int t = 0; t < timeMax; t++)
  {
    petsc.setValueZero();
    grid.dirichletBCs.assignBCs(grid.node, t);

    if (pulsatileFlow == ON)
    {
      if (t >= pulseBeginItr)
      {
        grid.dirichlet.assignPulsatileBCs(t, dt, T, pulseBeginItr, grid.nDofsGlobal);
      }
    }
    grid.dirichlet.applyDirichletBCs(grid.cell, petsc);

    MPI_Barrier(MPI_COMM_WORLD);
    double timer1 = MPI_Wtime();

    for (int ic = 0; ic < grid.cell.nCellsGlobal; ic++)
    {
      if (grid.cell(ic).subId == mpi.myId)
      {
        int nDofsInCell = grid.cell(ic).dofsMap.size();
        MathTools3D tools(grid.cell.nNodesInCell);
        MatrixXd Klocal(nDofsInCell, nDofsInCell);
        VectorXd Flocal(nDofsInCell);
        Klocal.setZero();
        Flocal.setZero();
        matrixAssemblyUSNS(Klocal, Flocal, tools, ic, t);
        petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap,
                       grid.cell(ic).dofsBCsMap, Klocal, Flocal);
      }
    }
    petsc.currentStatus = ASSEMBLY_OK;
    timer1 = MPI_Wtime() - timer1;

    MPI_Barrier(MPI_COMM_WORLD);

    double timer2 = MPI_Wtime();
    petsc.solve();
    timer2 = MPI_Wtime() - timer2;

    VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(vecSEQ, &arraySolnTmp);

    // update solution vector
    for (int id = 0; id < grid.nDofsGlobal; id++)
      petsc.solution[id] = arraySolnTmp[id];

    VecRestoreArray(vecSEQ, &arraySolnTmp);
    updateSolutions();

    updateSolutionsVTI();
    std::string binFile = outputDir + "/input/velocity_" + to_string(t) + ".bin";
    BIN::exportVectorDataBIN(binFile, grid.node.vvti);

    // visualize
    switch (grid.gridType)
    {
    case GridType::STRUCTURED:
      updateSolutionsVTI();
      outputSolutionsVTI("solution", t);
      compVorticity(t);
      break;
    case GridType::UNSTRUCTURED:
      outputSolutionsVTU("solution", t);
      break;
    default:
      PetscPrintf(MPI_COMM_WORLD, "\nUndifined gridType\n");
      exit(1);
      break;
    }

    if (mpi.myId == 0)
    {
      double timeNow = t * dt;
      printf("Assy: %fs | Solve: %fs | SimTime: %fs \n", timer1, timer2, timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);
}

/******************************************************
 * @brief Solve Unsteady Navier Stokes.
 *        This function takes boundary arguments
 *        for the purpose of armijo criteria.
 */
void DirectProblem::solveUSNS(std::vector<std::map<int, std::vector<double>>> &vDirichletTmp,
                              std::vector<std::map<int, double>> &pDirichletTmp, std::vector<std::vector<double>> &v0Tmp)
{
  PetscPrintf(MPI_COMM_WORLD, "\nMain Solver\n");

  PetscScalar *arraySolnTmp;
  Vec vecSEQ;
  VecScatter ctx;
  VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

  petsc.setMatAndVecZero(grid.cell);
  petsc.initialAssembly();
  setVariablesZero();

  for (int id = 0; id < grid.nDofsGlobal; id++)
  {
    grid.dirichlet.dirichletBCsValueNewInit[id] = 0e0;
    grid.dirichlet.dirichletBCsValueNew[id] = 0e0;
  }

  // Give initial velocity value
  for (int in = 0; in < grid.node.nNodesGlobal; in++)
  {
    for (int d = 0; d < dim; d++)
    {
      grid.node.v[in][d] = v0Tmp[in][d];
    }
  }

  int snapCount = 0;
  for (int t = 0; t < timeMax; t++)
  {
    petsc.setValueZero();
    grid.dirichlet.assignDirichletBCs(vDirichletTmp, pDirichletTmp,
                                      grid.node, dim, t);
    grid.dirichlet.applyDirichletBCs(grid.cell, petsc);

    for (int ic = 0; ic < grid.cell.nCellsGlobal; ic++)
    {
      if (grid.cell(ic).subId == mpi.myId)
      {
        int nDofsInCell = grid.cell(ic).dofsMap.size();
        MathTools3D tools(grid.cell.nNodesInCell);
        MatrixXd Klocal(nDofsInCell, nDofsInCell);
        VectorXd Flocal(nDofsInCell);
        Klocal.setZero();
        Flocal.setZero();
        matrixAssemblyUSNS(Klocal, Flocal, tools, ic, t);
        petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap,
                       grid.cell(ic).dofsBCsMap, Klocal, Flocal);
      }
    }
    petsc.solve();

    VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(vecSEQ, &arraySolnTmp);

    // update solution vector
    for (int id = 0; id < grid.nDofsGlobal; id++)
      petsc.solution[id] = arraySolnTmp[id];

    VecRestoreArray(vecSEQ, &arraySolnTmp);
    updateSolutions();
    updateTimeSolutions(t);

    if ((t - snap.snapTimeBeginItr) % snap.snapInterval == 0)
    {
      snap.takeSnapShot(grid.node.v, snapCount, grid.node.nNodesGlobal, dim);
      snapCount++;
    }

    if (mpi.myId == 0)
    {
      double timeNow = t * dt;
      printf("Main Solver : Time = %f \n", timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);
}
