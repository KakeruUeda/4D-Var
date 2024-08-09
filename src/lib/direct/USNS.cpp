/**
 * @file USNS.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "DirectProblem.h"

/**************************************
 * @brief Solve Unsteady Navier Stokes.
 */
void DirectProblem::solveUSNS(Application &app)
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
		grid.dirichlet.assignDirichletBCs(grid.dirichlet.vDirichletNew,
																			grid.dirichlet.pDirichletNew,
																			grid.node, dim, t);
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
		// PetscPrintf(MPI_COMM_WORLD, "\nMatrix assembly = %f seconds\n", timer);
		MPI_Barrier(MPI_COMM_WORLD);

		double timer2 = MPI_Wtime();
		petsc.solve();
		timer2 = MPI_Wtime() - timer2;

		// PetscPrintf(MPI_COMM_WORLD, "PETSc solver = %f seconds \n", timer);
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

/*******************************************************************
 * @brief Compute fully developed flow field.
 *        Usually used to get initial condition for inverse problem.
 */
void DirectProblem::compInitialCondition(std::vector<std::map<int, std::vector<double>>> &vDirichletTmp,
																				 std::vector<std::map<int, double>> &pDirichletTmp)
{
	PetscPrintf(MPI_COMM_WORLD, "\nCompute Initial Condition\n");

	double norm, norm0;
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

	int snapCount = 0;
	for (int t = 0; t < timeMax * 3; t++)
	{
		petsc.setValueZero();
		grid.dirichlet.assignConstantDirichletBCs(vDirichletTmp, pDirichletTmp, grid.node, dim, t);
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

		if (t == timeMax * 3 - 1)
		{
			for (int in = 0; in < grid.node.nNodesGlobal; in++)
			{
				for (int d = 0; d < dim; d++)
				{
					grid.node.v0[in][d] = grid.node.v[in][d];
				}
			}
		}

		if (mpi.myId == 0)
		{
			double timeNow = t * dt;
			printf("Compute initial condition : Time = %f \n", timeNow);
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}
	VecScatterDestroy(&ctx);
	VecDestroy(&vecSEQ);
}

/********************************
 * @brief Set all solutions zero.
 */
void DirectProblem::setVariablesZero()
{
	for (int in = 0; in < grid.node.nNodesGlobal; in++)
	{
		for (int d = 0; d < dim; d++)
		{
			grid.node.v[in][d] = 0e0;
			grid.node.vPrev[in][d] = 0e0;
		}
		grid.node.p[in] = 0e0;
	}

	for (int t = 0; t < timeMax; t++)
	{
		for (int in = 0; in < grid.node.nNodesGlobal; in++)
		{
			for (int d = 0; d < dim; d++)
			{
				grid.node.vt[t][in][d] = 0e0;
			}
			grid.node.pt[t][in] = 0e0;
		}
	}
}

/******************************
 * @brief Update all solutions.
 */
void DirectProblem::updateSolutions()
{
	for (int in = 0; in < grid.node.nNodesGlobal; in++)
	{
		int n1 = 0;
		for (int i = 0; i < grid.node.mapNew[in]; i++)
			n1 += grid.node.nDofsOnNode[i];

		for (int d = 0; d < dim; d++)
			grid.node.vPrev[in][d] = grid.node.v[in][d];

		for (int d = 0; d < dim; d++)
			grid.node.v[in][d] = petsc.solution[n1 + d];
		grid.node.p[in] = petsc.solution[n1 + dim];
	}
}

/**************************************************************
 * @brief Assign time-dependent variables for adjoint equation.
 */
void DirectProblem::updateTimeSolutions(const int t)
{
	for (int in = 0; in < grid.node.nNodesGlobal; in++)
	{
		for (int d = 0; d < dim; d++)
		{
			grid.node.vt[t][in][d] = grid.node.v[in][d];
		}
		grid.node.pt[t][in] = grid.node.p[in];
	}
}

/**********************************
 * @brief Update solutions for VTI.
 */
void DirectProblem::updateSolutionsVTI()
{
	for (int in = 0; in < grid.node.nNodesGlobal; in++)
	{
		for (int d = 0; d < dim; d++)
		{
			grid.node.vvti[grid.node.sortNode[in]][d] = grid.node.v[in][d];
		}
		grid.node.pvti[grid.node.sortNode[in]] = grid.node.p[in];
	}
}

/**********************************
 * @brief Update solutions for VTI.
 */
void DirectProblem::updateSolutionsVTI(const int t)
{
	for (int in = 0; in < grid.node.nNodesGlobal; in++)
	{
		for (int d = 0; d < dim; d++)
		{
			grid.node.vvti[grid.node.sortNode[in]][d] = grid.node.vt[t][in][d];
		}
		grid.node.pvti[grid.node.sortNode[in]] = grid.node.pt[t][in];
	}
}

/**********************************
 * @brief Export solutions for VTI.
 */
void DirectProblem::outputSolutionsVTI(const std::string &dir, const int t)
{
	if (mpi.myId > 0)
		return;

	std::string vtiFile;
	vtiFile = outputDir + "/" + dir + "/velocity_" + to_string(t) + ".vti";
	VTK::exportVectorPointDataVTI(vtiFile, "velocity", grid.node.vvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
	vtiFile = outputDir + "/" + dir + "/pressure_" + to_string(t) + ".vti";
	VTK::exportScalarPointDataVTI(vtiFile, "pressure", grid.node.pvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
}

/**********************************
 * @brief Export solutions for VTU.
 */
void DirectProblem::outputSolutionsVTU(const std::string &dir, const int t)
{
	if (mpi.myId > 0)
		return;

	std::string vtuFile;
	vtuFile = outputDir + "/" + dir + "/velocity_" + to_string(t) + ".vtu";
	VTK::exportVectorPointDataVTU(vtuFile, "velocity", grid.node, grid.cell, grid.node.v);
	vtuFile = outputDir + "/" + dir + "/pressure_" + to_string(t) + ".vtu";
	VTK::exportScalarPointDataVTU(vtuFile, "pressure", grid.node, grid.cell, grid.node.p);
}

/**********************************
 * @brief Export solutions for VTI.
 */
void DirectProblem::outputSolutionsVTI(const std::string &dir, const int t, const int loop)
{
	if (mpi.myId > 0)
		return;

	std::string vtiFile;
	vtiFile = outputDir + "/" + dir + "/velocity_" + to_string(loop) + "_" + to_string(t) + ".vti";
	VTK::exportVectorPointDataVTI(vtiFile, "velocity", grid.node.vvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
	vtiFile = outputDir + "/" + dir + "/pressure_" + to_string(loop) + "_" + to_string(t) + ".vti";
	VTK::exportScalarPointDataVTI(vtiFile, "pressure", grid.node.pvti, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
}

/**********************************
 * @brief Export solutions for VTU.
 */
void DirectProblem::outputSolutionsVTU(const std::string &dir, const int t, const int loop)
{
	if (mpi.myId > 0)
		return;

	std::string vtuFile;
	vtuFile = outputDir + "/" + dir + "/velocity_" + to_string(loop) + "_" + to_string(t) + ".vtu";
	VTK::exportVectorPointDataVTU(vtuFile, "velocity", grid.node, grid.cell, grid.node.v);
	vtuFile = outputDir + "/" + dir + "/velocity_" + to_string(loop) + "_" + to_string(t) + ".vtu";
	VTK::exportScalarPointDataVTU(vtuFile, "pressure", grid.node, grid.cell, grid.node.p);
}

/********************************************
 * @brief Take snapshots for error functions.
 */
void SnapShot::takeSnapShot(std::vector<std::vector<double>> &_v,
														const int &snapCount, const int &nNodesGlobal, const int &dim)
{
	for (int in = 0; in < nNodesGlobal; in++)
	{
		for (int d = 0; d < dim; d++)
		{
			v[snapCount][in][d] = _v[in][d];
		}
	}
}

void DirectProblem::compVorticity(const int t)
{
	if (mpi.myId != 0)
		return;

	std::vector<std::vector<double>> omega((grid.nx + 1) * (grid.ny + 1) * (grid.nz + 1), std::vector<double>(3, 0e0));
	for (int k = 1; k < grid.nz; k++)
	{
		for (int j = 1; j < grid.ny; j++)
		{
			for (int i = 1; i < grid.nx; i++)
			{
				int n = i + j * (grid.nx + 1) + k * (grid.nx + 1) * (grid.ny + 1);
				int n_iminus1 = i - 1 + j * (grid.nx + 1) + k * (grid.nx + 1) * (grid.ny + 1);
				int n_iplus1 = i + 1 + j * (grid.nx + 1) + k * (grid.nx + 1) * (grid.ny + 1);
				int n_jminus1 = i + (j - 1) * (grid.nx + 1) + k * (grid.nx + 1) * (grid.ny + 1);
				int n_jplus1 = i + (j + 1) * (grid.nx + 1) + k * (grid.nx + 1) * (grid.ny + 1);
				int n_kminus1 = i + j * (grid.nx + 1) + (k - 1) * (grid.nx + 1) * (grid.ny + 1);
				int n_kplus1 = i + j * (grid.nx + 1) + (k + 1) * (grid.nx + 1) * (grid.ny + 1);
				double dw_dy = (grid.node.vvti[n_jplus1][2] - grid.node.vvti[n_jminus1][2]) / (2.0 * grid.dy);
				double dv_dz = (grid.node.vvti[n_kplus1][1] - grid.node.vvti[n_kminus1][1]) / (2.0 * grid.dz);
				double du_dz = (grid.node.vvti[n_kplus1][0] - grid.node.vvti[n_kminus1][0]) / (2.0 * grid.dz);
				double dw_dx = (grid.node.vvti[n_iplus1][2] - grid.node.vvti[n_iminus1][2]) / (2.0 * grid.dx);
				double dv_dx = (grid.node.vvti[n_iplus1][1] - grid.node.vvti[n_iminus1][1]) / (2.0 * grid.dx);
				double du_dy = (grid.node.vvti[n_jplus1][0] - grid.node.vvti[n_jminus1][0]) / (2.0 * grid.dy);

				omega[n][0] = dw_dy - dv_dz;
				omega[n][1] = du_dz - dw_dx;
				omega[n][2] = dv_dx - du_dy;
			}
		}
	}

	std::string vtiFile;
	vtiFile = outputDir + "/solution/vorticity" + to_string(t) + ".vti";
	VTK::exportVectorPointDataVTI(vtiFile, "vorticity", omega, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
}
