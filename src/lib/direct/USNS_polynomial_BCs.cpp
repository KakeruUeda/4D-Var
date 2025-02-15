#include "DirectProblem.h"

/**
 * @brief Compute fully developed flow field.
 */
void DirectProblem::solve_NavierStokes_polynmial_BCs(DataGrid &data, ControlBoundary &cb, const int t_max)
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

  for(int t = 0; t < t_max; t++) {
    petsc.setValueZero();

    VectorXd coeff_x = Polynomial3D::fitPolynomial3D(data.x_slice, data.z_slice, data.u_interpolated_slice[0]);
    VectorXd coeff_y = Polynomial3D::fitPolynomial3D(data.x_slice, data.z_slice, data.v_interpolated_slice[0]);
    VectorXd coeff_z = Polynomial3D::fitPolynomial3D(data.x_slice, data.z_slice, data.w_interpolated_slice[0]);

    std::vector<double> x_new, y_new;
    std::vector<double> u_new, v_new, w_new;

    for(auto &entry : dirichlet.velocitySet) {
      std::vector<double> &velocities = entry.second;

      if(grid.node.x[entry.first][1] < 1e-12) {
        x_new.push_back(grid.node.x[entry.first][0]);
        y_new.push_back(grid.node.x[entry.first][2]);
      }
    }

    interpolateGrid(data, coeff_x, x_new, y_new, u_new);
    interpolateGrid(data, coeff_y, x_new, y_new, v_new);
    interpolateGrid(data, coeff_z, x_new, y_new, w_new);

    int tmp = 0;
    for(auto &entry : dirichlet.velocitySet) {
      std::vector<double> &velocities = entry.second;

      if(grid.node.x[entry.first][1] < 1e-12) {
        velocities[0] = u_new[tmp];
        velocities[1] = v_new[tmp];
        velocities[2] = w_new[tmp];
        tmp++;
      }
    }

    for(int ie = 0; ie < cb.CBEdgeNodeMap.size(); ie++) {
      int in = cb.CBEdgeNodeMap[ie];

      if(dirichlet.velocitySet.find(in) != dirichlet.velocitySet.end()) {
        dirichlet.velocitySet[in] = {0.0, 0.0, 0.0};
      }
    }

    dirichlet.getNewArray(grid.node.mapNew);

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
    outputTemporaryVelocityVTU("other", t);

    if(mpi.myId == 0) {
      double timeNow = t * dt;
      printf("Assy: %fs | Solve: %fs | SimTime: %fs \n", timer1, timer2, timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for(int in = 0; in < grid.node.nNodesGlobal; in++) {
    for(int d = 0; d < 3; d++) {
      v0(in, d) = v(in, d);
    }
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);
}


/**
 * @brief Compute fully developed flow field.
 */
void DirectProblem::solve_NavierStokes_polynmial_BCs(DataGrid &data, ControlBoundary &cb,
                                                     std::vector<std::map<int, std::vector<double>>> &vectorVelocitySet)
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

  for(int t = 0; t < timeMax; t++) {
    petsc.setValueZero();

    VectorXd coeff_x = Polynomial3D::fitPolynomial3D(data.x_slice, data.z_slice, data.u_interpolated_slice[t]);
    VectorXd coeff_y = Polynomial3D::fitPolynomial3D(data.x_slice, data.z_slice, data.v_interpolated_slice[t]);
    VectorXd coeff_z = Polynomial3D::fitPolynomial3D(data.x_slice, data.z_slice, data.w_interpolated_slice[t]);

    std::vector<double> x_new, y_new;
    std::vector<double> u_new, v_new, w_new;

    for(auto &entry : dirichlet.velocitySet) {
      std::vector<double> &velocities = entry.second;

      if(grid.node.x[entry.first][1] < 1e-12) {
        x_new.push_back(grid.node.x[entry.first][0]);
        y_new.push_back(grid.node.x[entry.first][2]);
      }
    }

    interpolateGrid(data, coeff_x, x_new, y_new, u_new);
    interpolateGrid(data, coeff_y, x_new, y_new, v_new);
    interpolateGrid(data, coeff_z, x_new, y_new, w_new);

    int tmp = 0;
    for(auto &entry : dirichlet.velocitySet) {
      std::vector<double> &velocities = entry.second;

      if(grid.node.x[entry.first][1] < 1e-12) {
        velocities[0] = u_new[tmp];
        velocities[1] = v_new[tmp];
        velocities[2] = w_new[tmp];
        tmp++;
      }
    }

    for(int ie = 0; ie < cb.CBEdgeNodeMap.size(); ie++) {
      int in = cb.CBEdgeNodeMap[ie];

      if(dirichlet.velocitySet.find(in) != dirichlet.velocitySet.end()) {
        dirichlet.velocitySet[in] = {0.0, 0.0, 0.0};
      }
    }

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

    /* update solution vector */
    for(int id = 0; id < grid.nDofsGlobal; id++) {
      petsc.solution[id] = arraySolution[id];
    }
    VecRestoreArray(vecSEQ, &arraySolution);
    updateSolutions();
    updateTimeSolutions(t);
    outputTemporaryVelocityVTU("other", t);

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
void DirectProblem::solve_NavierStokes_polynmial_BCs(DataGrid &data, ControlBoundary &cb, const int c_cycle,
                                                     std::vector<std::map<int, std::vector<double>>> &vectorVelocitySet)
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

  for(int t = 0; t < timeMax; t++) {
    petsc.setValueZero();

    VectorXd coeff_x = Polynomial3D::fitPolynomial3D(data.x_slice, data.z_slice, data.u_interpolated_slice[t]);
    VectorXd coeff_y = Polynomial3D::fitPolynomial3D(data.x_slice, data.z_slice, data.v_interpolated_slice[t]);
    VectorXd coeff_z = Polynomial3D::fitPolynomial3D(data.x_slice, data.z_slice, data.w_interpolated_slice[t]);

    std::vector<double> x_new, y_new;
    std::vector<double> u_new, v_new, w_new;

    for(auto &entry : dirichlet.velocitySet) {
      std::vector<double> &velocities = entry.second;

      if(grid.node.x[entry.first][1] < 1e-12) {
        x_new.push_back(grid.node.x[entry.first][0]);
        y_new.push_back(grid.node.x[entry.first][2]);
      }
    }

    interpolateGrid(data, coeff_x, x_new, y_new, u_new);
    interpolateGrid(data, coeff_y, x_new, y_new, v_new);
    interpolateGrid(data, coeff_z, x_new, y_new, w_new);

    int tmp = 0;
    for(auto &entry : dirichlet.velocitySet) {
      std::vector<double> &velocities = entry.second;

      if(grid.node.x[entry.first][1] < 1e-12) {
        velocities[0] = u_new[tmp];
        velocities[1] = v_new[tmp];
        velocities[2] = w_new[tmp];
        tmp++;
      }
    }

    for(int ie = 0; ie < cb.CBEdgeNodeMap.size(); ie++) {
      int in = cb.CBEdgeNodeMap[ie];

      if(dirichlet.velocitySet.find(in) != dirichlet.velocitySet.end()) {
        dirichlet.velocitySet[in] = {0.0, 0.0, 0.0};
      }
    }

    dirichlet.getNewArray(grid.node.mapNew);

    dirichlet.assignBCs(grid.node);
    dirichlet.applyBCs(grid.cell, petsc);

    if(c_cycle == 2) vectorVelocitySet.push_back(dirichlet.velocitySet);

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
    updateTimeSolutions(t);
    outputTemporaryVelocityVTU("other", t);

    if(mpi.myId == 0) {
      double timeNow = t * dt;
      printf("Assy: %fs | Solve: %fs | SimTime: %fs \n", timer1, timer2, timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);
}

void DirectProblem::interpolateGrid(DataGrid &data, const VectorXd &coeff, std::vector<double> &x_new,
                                    std::vector<double> &y_new, std::vector<double> &z_new)
{
  for(int i = 0; i < x_new.size(); ++i) {
    int mask_i = static_cast<int>(x_new[i] / data.dxData);
    int mask_j = static_cast<int>(y_new[i] / data.dzData);

    int iv = (mask_j * data.nxData * data.nyData) + mask_i;

    if(data.voxel(iv).mask < 1e-12) {
      z_new.push_back(0e0);
    } else {
      double xi = x_new[i];
      double yi = y_new[i];

      double tmp = coeff[0] + coeff[1] * xi + coeff[2] * yi + coeff[3] * xi * xi + coeff[4] * xi * yi +
                   coeff[5] * yi * yi + coeff[6] * xi * xi * xi + coeff[7] * xi * xi * yi + coeff[8] * xi * yi * yi +
                   coeff[9] * yi * yi * yi;
      z_new.push_back(tmp);
    }
  }
}