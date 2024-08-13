/**
 * @file DataGrid.h
 * @author K.Ueda
 * @date July, 2024
 */

#ifndef DATAGRID_H
#define DATAGRID_H

#include <iostream>
#include "metis.h"
#include "Config.h"
#include "Node.h"
#include "Cell.h"
#include "Gauss.h"
#include "ShapeFunction.h"
#include "PetscSolver.h"
#include "Tool.h"
#include "Function.h"
#include "MathCommon.h"

struct VoxelInfo
{
	int centerCell;
	bool isIncluded;
	int nNodesInCell;
	double dx,dy,dz;

	std::vector<std::vector<double>> vCFD;
	std::vector<std::vector<double>> vMRI;
	std::vector<std::vector<double>> ve;
	std::vector<double> center;
	std::vector<int> cellChildren;
	std::vector<double> vEX;

	void setNearCell(Node &node,Cell &cell,const double range,const int dim);
	void setCellOnCenterPoint(Node &node,Cell &cell,const int dim);
	void average(Cell &cell,std::vector<std::vector<double>> &_v,
	             const int t,const int dim);
	double compSmoothing(Function &func,const double a,const int dim,const int p);
	void interpolate(Node &node,Cell &cell,std::vector<std::vector<double>> &_v,
	                 const int &t,const int &dim);
	void gaussIntegral(Function &func,std::vector<std::vector<double>> &velCurrent,std::vector<double> &smoothing,
	                   double &weightIntegral,const int nNodesInCell,const int t,const int dim);
};

class DataGrid
{
public:
  DataGrid() {
  }
  DataGrid(Config &conf);
  
  int dim;
  int nx,ny,nz;
  double dx,dy,dz;
  double lx,ly,lz;
  double xOrigin,yOrigin,zOrigin;
  double range;
  
  int nData;
  
  int nCellsGlobal;
  int nNodesGlobal;
  int nNodesInCell;
  
  int nSnapShot;
  int snapInterval;
  
  std::vector<std::vector<std::vector<std::vector<double>> > > vEX;
  
  inline VoxelInfo &operator()(int x)
  {
  	return data[x];
  }
  inline VoxelInfo &operator()(int y,int x)
  {
  	return data[y*nx+x];
  }
  inline VoxelInfo &operator()(int z,int y,int x)
  {
  	return data[z*nx*ny+y*nx+x];
  }
  
  inline int size()
  {
  	return data.size();
  }
  inline void resize(int n)
  {
  	data.resize(n);
  }

  void initialize(Config &conf,Node &node,Cell &cell,const int &dim);
  void compEdgeValue(const int t);

private:
  std::vector<VoxelInfo> data;
};

#endif