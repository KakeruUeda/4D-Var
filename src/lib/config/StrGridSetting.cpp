/**
 * @file SetStrGrid.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

void Config::setStrGrid()
{
	cell.resize(nCellsGlobal);
	for(int ic=0; ic<nCellsGlobal; ic++)
	{
		cell[ic].resize(nNodesInCell);
	}

	node.resize(nNodesGlobal);
	for(int in=0; in<nNodesGlobal; in++)
	{
		node[in].resize(dim);
	}

	for(int k=0; k<nz; k++)
	{
		for(int j=0; j<ny; j++)
		{
			for(int i=0; i<nx; i++)
			{
				int ic=k*ny*nx+j*nx+i;
				for(int p=0; p<nNodesInCell; p++)
				{
					cell[ic][p]=setStrNode(i, j, k, p);
				}
			}
		}
	}

	for(int k=0; k<nz+1; k++)
	{
		for(int j=0; j<ny+1; j++)
		{
			for(int i=0; i<nx+1; i++)
			{
				int in=k*(ny+1)*(nx+1)+j*(nx+1)+i;
				for(int d=0; d<dim; d++)
				{
					node[in][d]=setStrCoordinate(i, j, k, d);
				}
			}
		}
	}
}

int Config::setStrNode(const int i, const int j, const int k, const int p)
{
	switch(p)
	{
	case 0:
		return i+j*(nx+1)+k*(nx+1)*(ny+1);
	case 1:
		return i+1+j*(nx+1)+k*(nx+1)*(ny+1);
	case 2:
		return i+1+(j+1)*(nx+1)+k*(nx+1)*(ny+1);
	case 3:
		return i+(j+1)*(nx+1)+k*(nx+1)*(ny+1);
	case 4:
		return i+j*(nx+1)+(k+1)*(nx+1)*(ny+1);
	case 5:
		return i+1+j*(nx+1)+(k+1)*(nx+1)*(ny+1);
	case 6:
		return i+1+(j+1)*(nx+1)+(k+1)*(nx+1)*(ny+1);
	case 7:
		return i+(j+1)*(nx+1)+(k+1)*(nx+1)*(ny+1);
	default:
		throw std::out_of_range("Invalid node index");
	}
}

double Config::setStrCoordinate(const int i, const int j, const int k, const int d)
{
	switch(d)
	{
	case 0:
		return i*dx;
	case 1:
		return j*dy;
	case 2:
		return k*dz;
	default:
		throw std::out_of_range("Invalid dimension index");
	}
}