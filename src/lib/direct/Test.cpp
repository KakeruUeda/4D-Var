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

  PetscScalar *arraySolution;
  Vec vecSEQ;
  VecScatter ctx;
  VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

  petsc.setMatAndVecZero(grid.cell);
  petsc.initialAssembly();

  dirichlet.setValuesZero(grid.nDofsGlobal);

  for(int t = 0; t < timeMax; t++) {
    petsc.setValueZero();
    dirichlet.assignBCs(grid.node, t);

    if(pulsatileFlow == ON) {
      if(t >= pulseBeginItr) {
        double pulse = comp_pulse(t);
        dirichlet.assignPulsatileBCs(pulse, grid.nDofsGlobal);
      }
    }
    dirichlet.applyBCs(grid.cell, petsc);

    MPI_Barrier(MPI_COMM_WORLD);
    double timer1 = MPI_Wtime();

    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      int nDofs = grid.cell(ic).dofsMap.size();

      if(grid.cell(ic).subId == mpi.myId) {
        MatrixXd Klocal(nDofs, nDofs);
        VectorXd Flocal(nDofs);
        matrixAssemblyUSNS(Klocal, Flocal, ic, t);
        petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap, 
                       grid.cell(ic).dofsBCsMap, Klocal, Flocal);
      }
    }
    timer1 = MPI_Wtime() - timer1;

    MPI_Barrier(MPI_COMM_WORLD);
    double timer2 = MPI_Wtime();

    petsc.solve();

    timer2 = MPI_Wtime() - timer2;

    VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(vecSEQ, &arraySolution);

    // update solution vector
    for(int id = 0; id < grid.nDofsGlobal; id++) {
      petsc.solution[id] = arraySolution[id];
    }

    VecRestoreArray(vecSEQ, &arraySolution);

    updateSolutions();
    outputSolutions(t);

    if(mpi.myId == 0) {
      double timeNow = t * dt;
      printf("Assy: %fs | Solve: %fs | SimTime: %fs \n", timer1, timer2, timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);
}

/******************************************************
 * @brief Solve Unsteady Navier Stokes
 *        for the purpose of checking armijo criteria.
 */
void DirectProblem::solveFowardNavierStokes(Array2D<double> &X0, Array3D<double> &X)
{
  PetscPrintf(MPI_COMM_WORLD, "\nMain Solver\n");

  PetscScalar *arraySolnTmp;
  Vec vecSEQ;
  VecScatter ctx;
  VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

  petsc.setMatAndVecZero(grid.cell);
  petsc.initialAssembly();
  setVariablesZero();

  updateInitialVelocity(X0);

  int snapCount = 0;
  for(int t = 0; t < timeMax; t++) {
    petsc.setValueZero();

    dirichlet.updateValues(X, t);
    dirichlet.getNewArray(grid.node.mapNew);

    dirichlet.assignBCs(grid.node, t);
    dirichlet.applyBCs(grid.cell, petsc);

    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      int nDofs = grid.cell(ic).dofsMap.size();

      if(grid.cell(ic).subId == mpi.myId) {
        MatrixXd Klocal(nDofs, nDofs);
        VectorXd Flocal(nDofs);
        matrixAssemblyUSNS(Klocal, Flocal, ic, t);
        petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap, grid.cell(ic).dofsBCsMap, Klocal, Flocal);
      }
    }
    petsc.solve();

    VecScatterBegin(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, petsc.solnVec, vecSEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(vecSEQ, &arraySolnTmp);

    // update solution vector
    for(int id = 0; id < grid.nDofsGlobal; id++) {
      petsc.solution[id] = arraySolnTmp[id];
    }

    VecRestoreArray(vecSEQ, &arraySolnTmp);
    updateSolutions();
    updateTimeSolutions(t);

    if((t - snap.snapTimeBeginItr) % snap.snapInterval == 0) {
      snap.takeSnapShot(grid.node.v, snapCount, grid.node.nNodesGlobal, dim);
      snapCount++;
    }

    if(mpi.myId == 0) {
      double timeNow = t * dt;
      printf("Main Solver : Time = %f \n", timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);
}

void DirectProblem::updateInitialVelocity(Array2D<double> &X0)
{
  for(int i = 0; i < grid.node.nNodesGlobal; i++) {
    for(int d = 0; d < dim; d++) {
      v(i, d) = X0(i, d);
      vPrev(i, d) = X0(i, d);
    }
  }
}

void DirectProblem::outputSolutions(const int t)
{
  switch(grid.gridType) {
  case GridType::STRUCTURED:
    updateSolutionsVTI();
    outputSolutionsVTI("solution", t);
    outputSolutionsBIN("input", t);
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
}

void DirectProblem::compVorticity(const int t)
{
  if(mpi.myId > 0) {
    return;
  }

  std::function<double(int, int, int)> u_func = [this](int i, int j, int k) {
    int in = i + j * (grid.nx + 1) + k * (grid.nx + 1) * (grid.ny + 1);
    return vvti(in, 0);
  };
  std::function<double(int, int, int)> v_func = [this](int i, int j, int k) {
    int in = i + j * (grid.nx + 1) + k * (grid.nx + 1) * (grid.ny + 1);
    return vvti(in, 1);
  };
  std::function<double(int, int, int)> w_func = [this](int i, int j, int k) {
    int in = i + j * (grid.nx + 1) + k * (grid.nx + 1) * (grid.ny + 1);
    return vvti(in, 2);
  };

  for(int k = 1; k < grid.nz; k++) {
    for(int j = 1; j < grid.ny; j++) {
      for(int i = 1; i < grid.nx; i++) {
        double dvdz = CDM::zDerivative(v_func, i, j, k, grid.dz);
        double dwdy = CDM::yDerivative(w_func, i, j, k, grid.dy);
        double dudz = CDM::zDerivative(u_func, i, j, k, grid.dz);
        double dwdx = CDM::xDerivative(w_func, i, j, k, grid.dx);
        double dvdx = CDM::xDerivative(v_func, i, j, k, grid.dx);
        double dudy = CDM::yDerivative(u_func, i, j, k, grid.dy);

        int in = i + j * (grid.nx + 1) + k * (grid.nx + 1) * (grid.ny + 1);

        vrt(in, 0) = dwdy - dvdz;
        vrt(in, 1) = dudz - dwdx;
        vrt(in, 2) = dvdx - dudy;
      }
    }
  }

  std::string vtiFile = outputDir + "/solution/vorticity" + std::to_string(t) + ".vti";
  EXPORT::exportVectorPointDataVTI<double>(vtiFile, "vorticity", vrt, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy,
                                           grid.dz);
}