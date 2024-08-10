/**
 * @file main.cpp
 * @author K.Ueda
 * @date Mar, 2024
 */

#include <iostream>
#include <mpi.h>
#include <map>
#include <memory>
#include <sys/stat.h>
#include "Config.h"
#include "MyMPI.h"
#include "Boundary.h"
#include "FileIO.h"
#include "Grid.h"
#include "PetscSolver.h"
MyMPI mpi;

int main()
{
  
}

/*
void setControlBoundary(Cell &cell, ControlBoundary &controlBoundary,
                        std::vector<int> &controlBoundaryMap,
                        std::vector<std::vector<int>> &controlNodeInCell,
                        std::vector<int> &controlCellMap,
                        int nxNodes, int nyNodes, int nzNodes,
                        int nxCells, int nyCells, int nzCells)
{
  int nNodesInCell = 8;
  if (controlBoundary == ControlBoundary::left)
  {
    for (int k = 0; k < nzNodes; k++)
    {
      for (int j = 0; j < nyNodes; j++)
      {
        for (int i = 0; i < nxNodes; i++)
        {
          if (i == 0)
          {
            controlBoundaryMap.push_back(i + j * nxNodes + k * nxNodes * nyNodes);
          }
        }
      }
    }
    for (int k = 0; k < nzCells; k++)
    {
      for (int j = 0; j < nyCells; j++)
      {
        for (int i = 0; i < nxCells; i++)
        {
          if (i == 0)
          {
            std::vector<int> vecTmp(4, 0);
            vecTmp[0] = cell(i + j * nxCells + k * nxCells * nyCells).node[0];
            vecTmp[1] = cell(i + j * nxCells + k * nxCells * nyCells).node[3];
            vecTmp[2] = cell(i + j * nxCells + k * nxCells * nyCells).node[7];
            vecTmp[3] = cell(i + j * nxCells + k * nxCells * nyCells).node[4];
            controlNodeInCell.push_back(vecTmp);
            controlCellMap.push_back(i + j * nxCells + k * nxCells * nyCells);
          }
        }
      }
    }
  }
}
*/
