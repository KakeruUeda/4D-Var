#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <set>
#include "metis.h"
#include "Config.h"
#include "Cell.h"
#include "Node.h"
#include "Array.h"
#include "PetscSolver.h"
#include "Boundary.h"
#include "Output.h"
#include "Config.h"
#include "Gauss.h"
#include "ShapeFunction.h"
#include "MathFEM.h"
 
class Grid
{
    public:
        Grid(){}
        Grid(Config &conf);
        virtual ~Grid(){}

        GridType gridType;
        Cell cell; 
        Node node;
        DirichletBoundary dirichlet;
        Output output;
        
        int nx, ny, nz;
        double lx, ly, lz, dx, dy, dz;

        int dim;
        int extractFluid;
        int rowStart, rowEnd;

        int nNodesGlobal, nCellsGlobal, nDofsGlobal;
        int nNodesLocal, nCellsLocal, nDofsLocal;

        void setStructuredGrid(
            const int nxCells, const int nyCells, const int nzCells, 
            const int nxNodes, const int nyNodes, const int nzNodes, 
            const double dx, const double dy, const double dz,
            const int nNodesInCell, const int dim, 
            Cell &cell, Node &node);

        void prepareMatrix(PetscSolver &petsc, std::string outputDir, const int timeMax);
        void setForSerial();
        void divideWholeGrid();
        void distributeToLocal(const int timeMax);

    private:
        int structuredGridNodeSet(
            const int nxNodes, const int nyNodes, const int nzNodes,
            const int i, const int j, const int k, const int p);
        double structuredGridCoordinateSet(
            const double dx, const double dy, const double dz,
            const int i, const int j, const int k, const int d);
        
};

#endif