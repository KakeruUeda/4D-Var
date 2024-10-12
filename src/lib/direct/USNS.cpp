/**
 * @file USNS.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "DirectProblem.h"

/**
 * @brief Solve Unsteady Navier Stokes.
 */
void DirectProblem::solveNavierStokes()
{
  PetscPrintf(MPI_COMM_WORLD, "\nMain Solver\n");

  PetscScalar *arraySolution;
  Vec vecSEQ;
  VecScatter ctx;
  VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

  petsc.setMatAndVecZero(grid.cell);
  petsc.initialAssembly();

  dirichlet.setValuesZero(grid.nDofsGlobal);
  dirichlet.velocitySetInit = dirichlet.velocitySet;

  for(int t = 0; t < timeMax; t++) {
    petsc.setValueZero();
    dirichlet.assignBCs(grid.node);

    if(pulsatileFlow == ON) {
      if(t >= pulseBeginItr) {
        double pulse = comp_pulse(t);
        //double pulse = dirichlet.comp_pulse(t * dt);
        dirichlet.assignPulsatileBCs(pulse);
        dirichlet.getNewArray(grid.node.mapNew);
        comp_Re(dirichlet.velocitySet);
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
        petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap, grid.cell(ic).dofsBCsMap, Klocal, Flocal);
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

    /* upadate solution vector */
    for(int id = 0; id < grid.nDofsGlobal; id++) {
      petsc.solution[id] = 0e0;
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

/**
 * @brief Compute fully developed flow field.
 */
void DirectProblem::solveNaveirStokes(const int stepMax, std::vector<std::array<double, 2>> &velArr)
{
  PetscPrintf(MPI_COMM_WORLD, "\nMain Solver Opt Initial\n");

  PetscScalar *arraySolution;
  Vec vecSEQ;
  VecScatter ctx;
  VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

  petsc.setMatAndVecZero(grid.cell);
  petsc.initialAssembly();

  dirichlet.setValuesZero(grid.nDofsGlobal);
  dirichlet.velocitySetInit = dirichlet.velocitySet;

  for(int t = 0; t < stepMax; t++) {
    petsc.setValueZero();

    dirichlet.assignBCs(grid.node);
    dirichlet.applyBCs(grid.cell, petsc);

    MPI_Barrier(MPI_COMM_WORLD);
    double timer1 = MPI_Wtime();

    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      int nDofs = grid.cell(ic).dofsMap.size();
      if(grid.cell(ic).subId == mpi.myId) {
        MatrixXd Klocal(nDofs, nDofs);
        VectorXd Flocal(nDofs);
        matrixAssemblyUSNS(Klocal, Flocal, ic, t);
        petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap, grid.cell(ic).dofsBCsMap, Klocal, Flocal);
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

    /* update solution vector */
    for(int id = 0; id < grid.nDofsGlobal; id++) {
      petsc.solution[id] = arraySolution[id];
    }
    VecRestoreArray(vecSEQ, &arraySolution);
    updateSolutions();

    if(mpi.myId == 0) {
      double timeNow = t * dt;
      printf("Initial : Step: %d/%d | Time: %f \n", t, stepMax, timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);

  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < 3; d++) {
      v0(in, d) = v(in, d);
    }
  }
}

/**
 * @brief tmp
 */
void DirectProblem::solveNaveirStokes(std::vector<std::array<double, 2>> &velArr, 
                                      std::vector<std::map<int, std::vector<double>>> &vectorVelocitySet)
{
  PetscPrintf(MPI_COMM_WORLD, "\nMain Solver\n");

  PetscScalar *arraySolution;
  Vec vecSEQ;
  VecScatter ctx;
  VecScatterCreateToAll(petsc.solnVec, &ctx, &vecSEQ);

  petsc.setMatAndVecZero(grid.cell);
  petsc.initialAssembly();

  dirichlet.setValuesZero(grid.nDofsGlobal);
  //dirichlet.velocitySetInit = dirichlet.velocitySet;

  int snapCount = 0;

  for(int t = 0; t < timeMax; t++) {
    petsc.setValueZero();
    
    double pulse = dirichlet.comp_pulse2(t * dt, velArr);

    dirichlet.assignPulsatileBCs(pulse);
    dirichlet.getNewArray(grid.node.mapNew);

    dirichlet.assignBCs(grid.node);
    dirichlet.applyBCs(grid.cell, petsc);

    vectorVelocitySet.push_back(dirichlet.velocitySet);

    MPI_Barrier(MPI_COMM_WORLD);
    double timer1 = MPI_Wtime();

    for(int ic = 0; ic < grid.cell.nCellsGlobal; ic++) {
      int nDofs = grid.cell(ic).dofsMap.size();

      if(grid.cell(ic).subId == mpi.myId) {
        MatrixXd Klocal(nDofs, nDofs);
        VectorXd Flocal(nDofs);
        matrixAssemblyUSNS(Klocal, Flocal, ic, t);
        petsc.setValue(grid.cell(ic).dofsBCsMap, grid.cell(ic).dofsMap, grid.cell(ic).dofsBCsMap, Klocal, Flocal);
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

    /* upadate solution vector */
    for(int id = 0; id < grid.nDofsGlobal; id++) {
      petsc.solution[id] = 0e0;
      petsc.solution[id] = arraySolution[id];
    }

    VecRestoreArray(vecSEQ, &arraySolution);

    updateSolutions();
    updateTimeSolutions(t);
    outputSolutionsVTU("other", t);

    if((t - snap.snapTimeBeginItr) % snap.snapInterval == 0) {
      snap.takeSnapShot(v, grid.node.nNodesGlobal, snapCount);
      snapCount++;
    }

    if(mpi.myId == 0) {
      double timeNow = t * dt;
      printf("Assy: %fs | Solve: %fs | SimTime: %fs \n", timer1, timer2, timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);
}

/**
 * @brief Solve Unsteady Navier Stokes
 *        for the purpose of checking armijo criteria.
 */
void DirectProblem::solveNavierStokes(Array2D<double> &X0, Array3D<double> &X)
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

    dirichlet.assignBCs(grid.node);
    dirichlet.applyBCs(grid.cell, petsc);

    comp_Re(dirichlet.velocitySet);

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

    /* update solution vector */
    for(int id = 0; id < grid.nDofsGlobal; id++) {
      petsc.solution[id] = arraySolnTmp[id];
    }

    VecRestoreArray(vecSEQ, &arraySolnTmp);
    updateSolutions();
    updateTimeSolutions(t);

    if((t - snap.snapTimeBeginItr) % snap.snapInterval == 0) {
      snap.takeSnapShot(v, grid.node.nNodesGlobal, snapCount);
      snapCount++;
    }

    if(mpi.myId == 0) {
      double timeNow = t * dt;
      printf("Main: Step: %d/%d | Time: %.2fs \n", t, timeMax, timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);
}

void DirectProblem::updateInitialVelocity(Array2D<double> &X0)
{
  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < dim; d++) {
      v(in, d) = X0(in, d);
      vPrev(in, d) = X0(in, d);
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
    //compVorticity(t);
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

/**
 * @brief Set all solutions zero.
 */
void DirectProblem::setVariablesZero()
{
  v.fillZero();
  vPrev.fillZero();
  p.fillZero();
  vt.fillZero();
  pt.fillZero();
}

/**
 * @brief Update all solutions.
 */
void DirectProblem::updateSolutions()
{
  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    int n1 = 0;
    for(int i = 0; i < grid.node.mapNew[in]; i++) {
      n1 += grid.node.nDofsOnNode[i];
    }

    for(int d = 0; d < dim; d++) {
      vPrev(in, d) = 0e0;
      vPrev(in, d) = v(in, d);
    }

    for(int d = 0; d < dim; d++) {
      v(in, d) = 0e0;
      v(in, d) = petsc.solution[n1 + d];
    }
    p(in) = petsc.solution[n1 + dim];
  }
}

/**
 * @brief Assign time-dependent variables for adjoint equation.
 */
void DirectProblem::updateTimeSolutions(const int t)
{
  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < dim; d++) {
      vt(t, in, d) = v(in, d);
    }
    pt(t, in) = p(in);
  }
}

/**
 * @brief Update solutions for VTI.
 */
void DirectProblem::updateSolutionsVTI()
{
  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < dim; d++) {
      vvti(grid.vecFluidUniqueNodes[in], d) = v(in, d);
    }
    pvti(grid.vecFluidUniqueNodes[in]) = p(in);
  }
}

/**
 * @brief Update solutions for VTI.
 */
void DirectProblem::updateSolutionsVTI(const int t)
{
  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < dim; d++) {
      vvti(grid.vecFluidUniqueNodes[in], d) = vt(t, in, d);
    }
    pvti(grid.vecFluidUniqueNodes[in]) = pt(t, in);
  }
}

/**
 * @brief Export solutions for VTI.
 */
void DirectProblem::outputSolutionsVTI(const std::string &dir, const int t)
{
  if(mpi.myId > 0) {
    return;
  }
  std::string vtiFile;
  vtiFile = outputDir + "/" + dir + "/velocity_" + to_string(t) + ".vti";
  EXPORT::exportVectorPointDataVTI(vtiFile, "velocity", vvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
  vtiFile = outputDir + "/" + dir + "/pressure_" + to_string(t) + ".vti";
  EXPORT::exportScalarPointDataVTI(vtiFile, "pressure", pvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
}

/**
 * @brief Export solutions for VTI.
 */
void DirectProblem::outputSolutionsVTI(const std::string &dir, const int t, const int loop)
{
  if(mpi.myId > 0) {
    return;
  }
  std::string vtiFile;
  vtiFile = outputDir + "/" + dir + "/velocity_" + to_string(loop) + "_" + to_string(t) + ".vti";
  EXPORT::exportVectorPointDataVTI(vtiFile, "velocity", vvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
  vtiFile = outputDir + "/" + dir + "/pressure_" + to_string(loop) + "_" + to_string(t) + ".vti";
  EXPORT::exportScalarPointDataVTI(vtiFile, "pressure", pvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
}

/**
 * @brief Export solutions for VTU.
 */
void DirectProblem::outputSolutionsVTU(const std::string &dir, const int t)
{
  if(mpi.myId > 0) {
    return;
  }
  std::string vtuFile;
  vtuFile = outputDir + "/" + dir + "/velocity_" + to_string(t) + ".vtu";
  EXPORT::exportVectorPointDataVTU<double>(vtuFile, "velocity", grid.node, grid.cell, v);
  vtuFile = outputDir + "/" + dir + "/pressure_" + to_string(t) + ".vtu";
  EXPORT::exportScalarPointDataVTU<double>(vtuFile, "pressure", grid.node, grid.cell, p);
}

/**
 * @brief Export solutions for VTU.
 */
void DirectProblem::outputSolutionsVTU(const std::string &dir, const int t, const int loop)
{
  if(mpi.myId > 0) {
    return;
  }
  std::string vtuFile;
  vtuFile = outputDir + "/" + dir + "/velocity_" + to_string(loop) + "_" + to_string(t) + ".vtu";
  EXPORT::exportVectorPointDataVTU<double>(vtuFile, "velocity", grid.node, grid.cell, grid.node.v);
  vtuFile = outputDir + "/" + dir + "/pressure_" + to_string(loop) + "_" + to_string(t) + ".vtu";
  EXPORT::exportScalarPointDataVTU<double>(vtuFile, "pressure", grid.node, grid.cell, grid.node.p);
}

/**
 * @brief Export solutions for BIN.
 */
void DirectProblem::outputSolutionsBIN(const std::string &dir, const int t)
{
  if(mpi.myId > 0) {
    return;
  }
  std::string binFile;
  binFile = outputDir + "/" + dir + "/velocity_" + to_string(t) + ".bin";
  vvti.exportBIN(binFile);
  binFile = outputDir + "/" + dir + "/pressure_" + to_string(t) + ".bin";
  pvti.exportBIN(binFile);
}

/**
 * @brief Take snapshots for error functions.
 */
void SnapShot::takeSnapShot(Array2D<double> &v, const int nNodesGlobal, const int snapCount)
{
  for(int in = 0; in < nNodesGlobal; in++) {
    for(int d = 0; d < 3; d++) {
      vSnap(snapCount, in, d) = v(in, d);
    }
  }
}

/**
 * @brief Take snapshots for error functions.
 */
void SnapShot::takeSnapShot(Array3D<double> &v, const int nNodesGlobal, const int snapCount, const int step)
{
  for(int in = 0; in < nNodesGlobal; in++) {
    for(int d = 0; d < 3; d++) {
      vSnap(snapCount, in, d) = v(step, in, d);
    }
  }
}