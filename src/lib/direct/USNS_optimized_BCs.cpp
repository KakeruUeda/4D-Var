#include "DirectProblem.h"

double DirectProblem::interpolated_velocity(Array2D<double> &vel_opt, double px, double pz, int dir)
{
  int nx_opt = 64;
  int ny_opt = 64;
  int nz_opt = 64;

  double lx_opt = 0.050625;
  double lz_opt = 0.0513;

  double dx_opt = lx_opt / nx_opt;
  double dz_opt = lz_opt / nz_opt;

  int ix = px / dx_opt + EPS;
  int iz = pz / dz_opt + EPS;

  if(px < 0 || px >= lx_opt || pz < 0 || pz >= lz_opt) {
    return 0e0;
  }

  double s = (px - ((ix * dx_opt) + dx_opt / 2e0)) / (dx_opt / 2e0);
  double u = (pz - ((iz * dz_opt) + dz_opt / 2e0)) / (dz_opt / 2e0);

  if(s < -1 - EPS || s > 1 + EPS || u < -1 - EPS || u > 1 + EPS) {
    if(mpi.myId == 0) {
      std::cout << "px = " << px << " pz = " << pz << " ix = " << ix << " iz = " << iz << " s = " << s << " u = " << u
                << std::endl;
    }
    throw std::runtime_error("Interpolation error: s or u out of bounds.");
  }

  Array1D<double> N(4);
  ShapeFunction2D::C2D4_N(N, s, u);

  double vel_int[3];

  int i1 = ix + iz * (nx_opt+1) * (ny_opt+1);
  int i2 = (ix + 1) + iz * (nx_opt+1) * (ny_opt+1);
  int i3 = (ix + 1) + (iz + 1) * (nx_opt+1) * (ny_opt+1);
  int i4 = ix + (iz + 1) * (nx_opt+1) * (ny_opt+1);

  for(int d = 0; d < 3; d++) {
    vel_int[d] = N(0) * vel_opt(i1, d) + N(1) * vel_opt(i2, d) +
                 N(2) * vel_opt(i3, d) + N(3) * vel_opt(i4, d);
  }

  if(dir == 0) return vel_int[0];
  if(dir == 1) return vel_int[1];
  if(dir == 2) return vel_int[2];

  return 0e0;
}

/**
 * @brief Compute fully developed flow field.
 */
void DirectProblem::solve_NavierStokes_optimized_BCs()
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

  std::vector<Array2D<double>> velocity_opt;
  velocity_opt.resize(156);

  for(int t = 0; t < velocity_opt.size(); t++) {
    std::string file_dir = "../../inverse/4dvar/output/Ubend_inlet_space_wave_time_wave_reg1e-1/optimized/";
    std::string file_name = file_dir + "velocity_" + to_string(t) + ".bin";
    velocity_opt[t].importBIN(file_name);
  }

  for(int t = 0; t < timeMax; t++) {
    petsc.setValueZero();

    int tc = 0;

    if(t >= 0 && t < 156) tc = t;
    if(t >= 156 && t < 312) tc = t-156;
    if(t >= 312 && t < 468) tc = t-312;
    if(t >= 468 && t < timeMax) tc = t-468;
    for(auto &entry : dirichlet.velocitySet) {
      std::vector<double> &velocities = entry.second;

      double x_coord = grid.node.x[entry.first][0];
      double y_coord = grid.node.x[entry.first][1];
      double z_coord = grid.node.x[entry.first][2];

      if(y_coord < 1e-12) {
        velocities[0] = interpolated_velocity(velocity_opt[tc], x_coord, z_coord, 0);
        velocities[1] = interpolated_velocity(velocity_opt[tc], x_coord, z_coord, 1);
        velocities[2] = interpolated_velocity(velocity_opt[tc], x_coord, z_coord, 2);
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
    outputSolutions(t);

    if(mpi.myId == 0) {
      double timeNow = t * dt;
      printf("Assy: %fs | Solve: %fs | SimItr: %d | SimTime: %fs \n", timer1, timer2, t, timeNow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vecSEQ);
}
