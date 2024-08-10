/**
 * @file SetStrGrid.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

void Config::setStrGrid()
{
	cell.resize(nCellsGlobal);
	for(int ic = 0; ic < nCellsGlobal; ic++)
	{
		cell[ic].resize(nNodesInCell);
	}

	node.resize(nNodesGlobal);
	for(int in = 0; in < nNodesGlobal; in++)
	{
		node[in].resize(dim);
	}

	for (int k = 0; k < nz; k++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int i = 0; i < nx; i++)
			{
        int ic = k * ny * nx + j * nx + i;
				for (int p = 0; p < nNodesInCell; p++)
				{
          cell[ic][p] = setStrNode(i, j, k, p);
				}
			}
		}
	}

	for (int k = 0; k < nz+1; k++)
	{
		for (int j = 0; j < ny+1; j++)
		{
			for (int i = 0; i < nx+1; i++)
			{
        int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
				for (int d = 0; d < dim; d++)
				{
          node[in][d] = setStrCoordinate(i, j, k, d);
				}
			}
		}
	}
}

int Config::setStrNode(const int i, const int j, const int k, const int p)
{
	std::vector<int> nodeSet(8);

	nodeSet[0] = i + j * (nx+1) + k * (nx+1) * (ny+1);
	nodeSet[1] = i + 1 + j * (nx+1) + k * (nx+1) * (ny+1);
	nodeSet[2] = i + 1 + (j + 1) * (nx+1) + k * (nx+1) * (ny+1);
	nodeSet[3] = i + (j + 1) * (nx+1) + k * (nx+1) * (ny+1);
	nodeSet[4] = i + j * (nx+1) + (k + 1) * (nx+1) * (ny+1);
	nodeSet[5] = i + 1 + j * (nx+1) + (k + 1) * (nx+1) * (ny+1);
	nodeSet[6] = i + 1 + (j + 1) * (nx+1) + (k + 1) * (nx+1) * (ny+1);
	nodeSet[7] = i + (j + 1) * (nx+1) + (k + 1) * (nx+1) * (ny+1);

	return nodeSet[p];
}

double Config::setStrCoordinate(const int i, const int j, const int k, const int d)
{
	std::vector<double> coordinateSet(3);

	coordinateSet[0] = i * dx;
	coordinateSet[1] = j * dy;
	coordinateSet[2] = k * dz;

	return coordinateSet[d];
}
