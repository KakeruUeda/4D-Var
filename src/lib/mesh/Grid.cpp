/**
 * @file Grid.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "Grid.h"

Grid::Grid(Config &conf) : gridType(conf.gridType), cell(conf), node(conf), dirichlet(conf),
													 nNodesGlobal(conf.nNodesGlobal), nCellsGlobal(conf.nCellsGlobal),
													 nDofsGlobal(0), nDofsLocal(0), dim(conf.dim)
{
	if (gridType == GridType::STRUCTURED)
	{
		nx = conf.nx;
		ny = conf.ny;
		nz = conf.nz;
		lx = conf.lx;
		ly = conf.ly;
		lz = conf.lz;
		dx = conf.dx;
		dy = conf.dy;
		dz = conf.dz;
	}
}

void Grid::setStructuredGrid(const int nxCells, const int nyCells, const int nzCells,
														 const int nxNodes, const int nyNodes, const int nzNodes,
														 const double dx, const double dy, const double dz,
														 const int nNodesInCell, const int dim,
														 Cell &cell, Node &node)
{
	for (int k = 0; k < nzCells; k++)
	{
		for (int j = 0; j < nyCells; j++)
		{
			for (int i = 0; i < nxCells; i++)
			{
				for (int p = 0; p < nNodesInCell; p++)
				{
					cell(k * nxCells * nyCells + j * nxCells + i).node[p] 
					= structuredGridNodeSet(nxNodes, nyNodes, nzNodes, i, j, k, p);
				}
			}
		}
	}

	for (int k = 0; k < nzNodes; k++)
	{
		for (int j = 0; j < nyNodes; j++)
		{
			for (int i = 0; i < nxNodes; i++)
			{
				for (int d = 0; d < dim; d++)
				{
					node.x[k * nxNodes * nyNodes + j * nxNodes + i][d] 
					= structuredGridCoordinateSet(dx, dy, dz, i, j, k, d);
				}
			}
		}
	}
}

int Grid::structuredGridNodeSet(const int nxNodes, const int nyNodes, const int nzNodes,
																const int i, const int j, const int k, const int p)
{
	std::vector<int> nodeSet(8);
	nodeSet[0] = i + j * nxNodes + k * nxNodes * nyNodes;
	nodeSet[1] = i + 1 + j * nxNodes + k * nxNodes * nyNodes;
	nodeSet[2] = i + 1 + (j + 1) * nxNodes + k * nxNodes * nyNodes;
	nodeSet[3] = i + (j + 1) * nxNodes + k * nxNodes * nyNodes;
	nodeSet[4] = i + j * nxNodes + (k + 1) * nxNodes * nyNodes;
	nodeSet[5] = i + 1 + j * nxNodes + (k + 1) * nxNodes * nyNodes;
	nodeSet[6] = i + 1 + (j + 1) * nxNodes + (k + 1) * nxNodes * nyNodes;
	nodeSet[7] = i + (j + 1) * nxNodes + (k + 1) * nxNodes * nyNodes;

	return nodeSet[p];
}

double Grid::structuredGridCoordinateSet(const double dx, const double dy, const double dz,
																				 const int i, const int j, const int k, const int d)
{
	std::vector<double> coordinateSet(3);
	coordinateSet[0] = i * dx;
	coordinateSet[1] = j * dy;
	coordinateSet[2] = k * dz;

	return coordinateSet[d];
}

void Grid::prepareMatrix(PetscSolver &petsc, std::string outputDir, const int timeMax)
{
	for (auto &pair : dirichlet.vDirichlet[0])
	{
		int count = 0;
		for (auto &value : pair.second)
		{
			node.isDirichlet[pair.first][count] = true;
			count++;
		}
	}

	for (auto &pair : dirichlet.pDirichlet[0])
		node.isDirichlet[pair.first][dim] = true;

	for (int in = 0; in < node.nNodesGlobal; in++)
		for (int id = 0; id < node.nDofsOnNode[in]; id++)
			if (node.isDirichlet[in][id])
				node.dofsBCsMap[in][id] = -1;

	for (int ic = 0; ic < cell.nCellsGlobal; ic++)
	{
		VecTool::resize(cell(ic).dofStart, cell.nNodesInCell);
		for (int p = 1; p < cell.nNodesInCell; p++)
		{
			cell(ic).dofStart[p] = cell(ic).dofStart[p - 1] + node.nDofsOnNode[cell(ic).node[p - 1]];
		}
	}

	//// add //////////////////////////////////////
	for (auto &pair : dirichlet.vDirichletWall[0])
	{
		int count = 0;
		for (auto &value : pair.second)
		{
			node.isDirichletWall[pair.first][count] = true;
			count++;
		}
	}
	for (int in = 0; in < node.nNodesGlobal; in++)
		for (int id = 0; id < node.nDofsOnNode[in]; id++)
			if (node.isDirichletWall[in][id])
				node.dofsBCsMapWall[in][id] = -1;
	////////////////////////////////////////////////

	nDofsGlobal = 0;
	for (int in = 0; in < node.nNodesGlobal; in++)
	{
		for (int id = 0; id < node.nDofsOnNode[in]; id++)
		{
			nDofsGlobal++;
		}
	}

	if (mpi.nId == 1)
	{
		setForSerial();
	}
	else if (mpi.nId > 1)
	{
		//divideWholeGrid();
		distributeToLocal(timeMax);
	}

	for (int in = 0; in < node.nNodesGlobal; in++)
	{
		for (int id = 0; id < node.nDofsOnNodeNew[in]; id++)
		{
			node.dofsMapNew1D.push_back(node.dofsMapNew[in][id]);
			node.dofsBCsMapNew1D.push_back(node.dofsBCsMapNew[in][id]);
		}
	}

	int size;
	for (int ic = 0; ic < cell.nCellsGlobal; ic++)
	{
		if (cell(ic).subId == mpi.myId)
		{
			size = 0;
			for (int p = 0; p < cell.nNodesInCell; p++)
			{
				size += node.nDofsOnNodeNew[cell(ic).nodeNew[p]];
			}
			cell(ic).dofsMap.resize(size);
			cell(ic).dofsBCsMap.resize(size);
			int i = 0;
			int j = 0;
			for (int p = 0; p < cell.nNodesInCell; p++)
			{
				j = cell(ic).nodeNew[p];
				for (int q = 0; q < node.nDofsOnNodeNew[cell(ic).nodeNew[p]]; q++)
				{
					cell(ic).dofsMap[i + q] = node.dofsMapNew[j][q];
					cell(ic).dofsBCsMap[i + q] = node.dofsBCsMapNew[j][q];
				}
				i += node.nDofsOnNodeNew[cell(ic).nodeNew[p]];
			}
		}
	}

	/// add ///
	for (int ic = 0; ic < cell.nCellsGlobal; ic++)
	{
		if (cell(ic).subId == mpi.myId)
		{
			size = 0;
			for (int p = 0; p < cell.nNodesInCell; p++)
			{
				size += node.nDofsOnNodeNew[cell(ic).nodeNew[p]];
			}
			cell(ic).dofsMapWall.resize(size);
			cell(ic).dofsBCsMapWall.resize(size);
			int i = 0;
			int j = 0;
			for (int p = 0; p < cell.nNodesInCell; p++)
			{
				j = cell(ic).nodeNew[p];
				for (int q = 0; q < node.nDofsOnNodeNew[cell(ic).nodeNew[p]]; q++)
				{
					cell(ic).dofsMapWall[i + q] = node.dofsMapWallNew[j][q];
					cell(ic).dofsBCsMapWall[i + q] = node.dofsBCsMapWallNew[j][q];
				}
				i += node.nDofsOnNodeNew[cell(ic).nodeNew[p]];
			}
		}
	}
	//////////

	int *tt, tmpInt;
	int r, kk, nSize;
	int countDiag, countOffDiag;

	std::vector<std::set<int>> forAssyMatFluid;
	std::set<int>::iterator it;

	forAssyMatFluid.resize(nDofsGlobal);

	for (int ic = 0; ic < cell.nCellsGlobal; ic++)
	{
		if (cell(ic).subId == mpi.myId)
		{
			tt = &(cell(ic).dofsMap[0]);
			nSize = cell(ic).dofsMap.size();

			for (int i = 0; i < nSize; i++)
			{
				r = tt[i];
				if (r != -1)
				{
					if (r >= rowStart && r <= rowEnd)
					{
						for (int j = 0; j < nSize; j++)
						{
							if (tt[j] != -1)
							{
								forAssyMatFluid[r].insert(tt[j]);
							}
						}
					}
				}
			}
		}
	}
	PetscMalloc1(nDofsLocal, &petsc.diag_nnz);
	PetscMalloc1(nDofsLocal, &petsc.offdiag_nnz);

	kk = 0;
	petsc.nnz_max_row = 0;
	for (int i = rowStart; i <= rowEnd; i++)
	{
		nSize = forAssyMatFluid[i].size();
		petsc.nnz_max_row = std::max(petsc.nnz_max_row, nSize);
		countDiag = 0, countOffDiag = 0;
		for (it = forAssyMatFluid[i].begin(); it != forAssyMatFluid[i].end(); it++)
		{
			tmpInt = *it;
			if (tmpInt >= rowStart && tmpInt <= rowEnd)
				countDiag++;
			else
				countOffDiag++;
		}
		petsc.diag_nnz[kk] = countDiag;
		petsc.offdiag_nnz[kk] = countOffDiag;
		kk++;
	}

	petsc.initialize(nDofsLocal, nDofsGlobal);
}

void Grid::setForSerial()
{
	nCellsLocal = cell.nCellsGlobal;
	nNodesLocal = node.nNodesGlobal;
	nDofsLocal = nDofsGlobal;
	rowStart = 0;
	rowEnd = nDofsGlobal - 1;

	node.isDirichletNew = node.isDirichlet;
	node.dofsBCsMapNew = node.dofsBCsMap;
	node.dofsMapNew = node.dofsMap;
	node.nDofsOnNodeNew = node.nDofsOnNode;
	dirichlet.vDirichletNew = dirichlet.vDirichletNew;

	for (int ic = 0; ic < cell.nCellsGlobal; ic++)
		cell(ic).nodeNew = cell(ic).node;
}

void Grid::divideWholeGrid()
{
	std::vector<int> cellId(nCellsGlobal, 0);
	std::vector<int> nodeId(nNodesGlobal, 0);

	if (mpi.myId == 0)
	{
		int kk = 0;
		int nparts = mpi.nId;

		std::unique_ptr<int[]> eptr;
		std::unique_ptr<int[]> eind;

		eptr = std::make_unique<int[]>(nCellsGlobal + 1);
		eind = std::make_unique<int[]>(nCellsGlobal * cell.nNodesInCell);

		for (int in = 0; in < nCellsGlobal + 1; in++)
			eptr[in] = 0;

		for (int ic = 0; ic < nCellsGlobal; ic++)
			for (int p = 0; p < cell.nNodesInCell; p++)
				eind[ic + p] = 0;

		for (int ic = 0; ic < nCellsGlobal; ic++)
		{
			eptr[ic + 1] = (ic + 1) * cell.nNodesInCell;
			for (int p = 0; p < cell.nNodesInCell; p++)
				eind[kk + p] = cell(ic).node[p];
			kk += cell.nNodesInCell;
		}
		// 8-noded hexa element
		int ncommon_nodes = 4;

		idx_t objval;
		idx_t options[METIS_NOPTIONS];

		METIS_SetDefaultOptions(options);

		// Specifies the partitioning method.
		options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;

		// Total communication volume minimization
		options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

		// C-style numbering is assumed that starts from 0
		options[METIS_OPTION_NUMBERING] = 0;

		// METIS partition routine
		int ret = METIS_PartMeshDual(&nCellsGlobal, &nNodesGlobal, eptr.get(), eind.get(), NULL,
																 NULL, &ncommon_nodes, &nparts, NULL, options,
																 &objval, &cellId[0], &nodeId[0]);

		if (ret == METIS_OK)
			std::cout << "METIS partition routine success " << std::endl;
		else
			std::cout << "METIS partition routine failed " << std::endl;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(cellId.data(), nCellsGlobal, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(nodeId.data(), nNodesGlobal, MPI_INT, 0, MPI_COMM_WORLD);

	for (int ic = 0; ic < nCellsGlobal; ic++)
		cell(ic).subId = cellId[ic];
	for (int in = 0; in < nNodesGlobal; in++)
		node.subId[in] = nodeId[in];

	nCellsLocal = count(cellId.begin(), cellId.end(), mpi.myId);
	nNodesLocal = count(nodeId.begin(), nodeId.end(), mpi.myId);

	printf("nCellsLocal =  %5d \t numOfId = %5d \t myId = %5d \n", nCellsLocal, mpi.nId, mpi.myId);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("nNodesLocal =  %5d \t numOfId = %5d \t myId = %5d \n", nNodesLocal, mpi.nId, mpi.myId);
}

void Grid::distributeToLocal(const int timeMax)
{
	int kk = 0;
	std::vector<int> nodeListLocal(node.nNodesLocal);

	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "e " << mpi.myId << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);

	for (int in = 0; in < nNodesGlobal; in++)
	{
		if (node.subId[in] == mpi.myId)
		{
			nodeListLocal[kk++] = in;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << mpi.myId  << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);

	std::vector<int> nNodesLocalVector(mpi.nId);
	std::vector<int> nNodesLocalSum(mpi.nId);

	MPI_Allgather(&node.nNodesLocal, 1, MPI_INT, &nNodesLocalVector[0], 1, MPI_INT, MPI_COMM_WORLD);
std::cout << "2_2" << std::endl;
	nNodesLocalSum = nNodesLocalVector;
	for (int i = 1; i < mpi.nId; i++)
		nNodesLocalSum[i] += nNodesLocalSum[i - 1];

	int nodeStart = 0;
	int nodeEnd = 0;

	if (mpi.myId > 0)
		nodeStart = nNodesLocalSum[mpi.myId - 1];
	nodeEnd = nNodesLocalSum[mpi.myId] - 1;
std::cout << "2_3" << std::endl;
	printf("nodeStart = %5d \t nodeEnd = %5d \t myId = %5d \n", nodeStart, nodeEnd, mpi.myId);

	std::vector<int> displs(mpi.nId);
std::cout << "2_4" << std::endl;
	displs[0] = 0;
	for (int i = 0; i < mpi.nId - 1; i++)
		displs[i + 1] = displs[i] + nNodesLocalVector[i];

	std::vector<int> tmp(nNodesGlobal);

	MPI_Allgatherv(&nodeListLocal[0], node.nNodesLocal, MPI_INT, &node.map[0], &nNodesLocalVector[0], &displs[0], MPI_INT, MPI_COMM_WORLD);
std::cout << "2_5" << std::endl;
	node.initializeNew();

	for (int ic = 0; ic < nCellsGlobal; ic++)
		for (int p = 0; p < cell.nNodesInCell; p++)
			cell(ic).nodeNew[p] = node.mapNew[cell(ic).node[p]];
std::cout << "2_6" << std::endl;
	int n1;
	int count;

	dirichlet.vDirichletNew.resize(timeMax);
	dirichlet.pDirichletNew.resize(timeMax);

	/// add ///
	dirichlet.vDirichletWallNew.resize(timeMax);
	//////////
std::cout << "2_7" << std::endl;
	for (int t = 0; t < timeMax; t++)
	{
		for (auto &pair : dirichlet.vDirichlet[t])
		{
			std::vector<double> vecTmp;
			n1 = node.mapNew[pair.first];
			count = 0;
			for (auto &value : pair.second)
			{
				vecTmp.push_back(value);
				if (t == 0)
					node.isDirichletNew[n1][count] = true;
				count++;
			}
			dirichlet.vDirichletNew[t][n1] = vecTmp;
		}
		/// add ///
		for (auto &pair : dirichlet.vDirichletWall[t])
		{
			std::vector<double> vecTmp;
			n1 = node.mapNew[pair.first];
			count = 0;
			for (auto &value : pair.second)
			{
				vecTmp.push_back(value);
				if (t == 0)
					node.isDirichletWallNew[n1][count] = true;
				count++;
			}
			dirichlet.vDirichletWallNew[t][n1] = vecTmp;
		}
		//////////

		for (auto &pair : dirichlet.pDirichlet[t])
		{
			n1 = node.mapNew[pair.first];
			dirichlet.pDirichletNew[t][n1] = pair.second;
			if (t == 0)
				node.isDirichletNew[n1][dim] = true;
		}
	}
std::cout << "2_8" << std::endl;
	for (int in = 0; in < node.nNodesGlobal; in++)
		for (int id = 0; id < node.nDofsOnNodeNew[in]; id++)
			if (node.isDirichletNew[in][id])
				node.dofsBCsMapNew[in][id] = -1;

	/// add ///
	for (int in = 0; in < node.nNodesGlobal; in++)
		for (int id = 0; id < node.nDofsOnNodeNew[in]; id++)
			if (node.isDirichletWallNew[in][id])
				node.dofsBCsMapWallNew[in][id] = -1;
	///////////

	rowStart = 0;
	for (int in = 0; in < nodeStart; in++)
		rowStart += node.nDofsOnNodeNew[in];

	rowEnd = 0;
	for (int in = 0; in <= nodeEnd; in++)
		rowEnd += node.nDofsOnNodeNew[in];

	rowEnd = rowEnd - 1;

	for (int in = nodeStart; in <= nodeEnd; in++)
		for (int id = 0; id < node.nDofsOnNodeNew[in]; id++)
			nDofsLocal++;

	printf("nDofsLocal = %5d/%5d \t rowStart  = %5d \t rowEnd  = %5d \t myId  = %5d \n",
				 nDofsLocal, nDofsGlobal, rowStart, rowEnd, mpi.myId);

	int nDofsGlobalCheck;
	MPI_Allreduce(&nDofsLocal, &nDofsGlobalCheck, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if (nDofsGlobalCheck != nDofsGlobal)
		std::cout << "Sum of local problem sizes is not equal to global size" << std::endl;
}
