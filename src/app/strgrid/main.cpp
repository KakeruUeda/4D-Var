/**
 * @file main.cpp
 * @author k.ueda
 * @date Mar, 2024
*/

#include <iostream>
#include <mpi.h>
#include <map>
#include <memory>
#include <sys/stat.h>
#include "TextParser.h"
#include "Config.h"
#include "MyMPI.h"
#include "Boundary.h"
#include "Output.h"
#include "Grid.h"
#include "PetscSolver.h"
MyMPI mpi;

void setControlBoundary(Cell &cell, ControlBoundary &controlBoundary, 
                        std::vector<int> &controlBoundaryMap,
                        std::vector<std::vector<int>> &controlNodeInCell,
                        std::vector<int> &controlCellMap,
                        int nxNodes, int nyNodes, int nzNodes,
                        int nxCells, int nyCells, int nzCells)
{
   int nNodesInCell = 8;
    if(controlBoundary == ControlBoundary::left){
        for(int k=0; k<nzNodes; k++){
            for(int j=0; j<nyNodes; j++){
                for(int i=0; i<nxNodes; i++){
                    if(i == 0){
                        controlBoundaryMap.push_back(i+j*nxNodes+k*nxNodes*nyNodes);
                    }
                }
            }
        }
        for(int k=0; k<nzCells; k++){
            for(int j=0; j<nyCells; j++){
                for(int i=0; i<nxCells; i++){
                    if(i == 0){
                        std::vector<int> vecTmp(4, 0);
                        vecTmp[0] = cell(i+j*nxCells+k*nxCells*nyCells).node[0];
                        vecTmp[1] = cell(i+j*nxCells+k*nxCells*nyCells).node[3];
                        vecTmp[2] = cell(i+j*nxCells+k*nxCells*nyCells).node[7];
                        vecTmp[3] = cell(i+j*nxCells+k*nxCells*nyCells).node[4];
                        controlNodeInCell.push_back(vecTmp);
                        controlCellMap.push_back(i+j*nxCells+k*nxCells*nyCells);
                    }
                }
            }
        }
    }
}

int main(int argc, char* argv[])
{ 
    // Based on Message Passing Interface
    MPI_Init(NULL, NULL); 
    mpi.setSizeAndRank();

    std::string inputFile = argv[1]; 
    std::string appName = "STRGRID";
    Config conf(inputFile, appName);
    if(conf.isReadingError) return EXIT_FAILURE;

    std::string output = "output";
    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    conf.outputDir = "output/" + conf.outputDir;
    if(mpi.myId == 0)
        mkdir(conf.outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    
    Grid grid; 

    grid.cell.resize(conf.nCellsGlobal);
    for(int ic=0; ic<conf.nCellsGlobal; ic++)
        grid.cell(ic).node.resize(conf.nNodesInCell);

    grid.node.x.resize(conf.nNodesGlobal);
    for(int in=0; in<conf.nNodesGlobal; in++)
        grid.node.x[in].resize(conf.dim);

    // Set structured grid components (cell, coordinate informateion)
    grid.setStructuredGrid(conf.nxCells, conf.nyCells, conf.nzCells, conf.nxNodes, 
                           conf.nyNodes, conf.nzNodes, conf.dx, conf.dy, conf.dz, 
                           conf.nNodesInCell, conf.dim, grid.cell, grid.node);

     // Define structured boundary faces set
    std::vector<StructuredBoundaryFace> bdFaces;
    bdFaces.reserve(conf.bdStr.size());

    for(auto str : conf.bdStr)
        bdFaces.emplace_back(str);
    for(auto &vec : bdFaces)
        vec.setNodesOnBoundaryFace(conf.nxNodes, conf.nyNodes, conf.nzNodes);

    int iteration = 0;

    // Insert dirichlet info from config parameter
    for(auto &vec : bdFaces){
        vec.setDirichletInfo(conf.bdType, conf.bdValue, conf.dim, iteration);
        iteration++;
    }

    ControlBoundary controlBoundary = conf.controlBoundary;
    std::vector<int> controlBoundaryMap;
    std::vector<std::vector<int>> controlNodeInCell;
    std::vector<int> controlCellMap;
    setControlBoundary(grid.cell, controlBoundary, controlBoundaryMap, controlNodeInCell, controlCellMap,
                       conf.nxNodes, conf.nyNodes, conf.nzNodes, conf.nxCells, conf.nyCells, conf.nzCells);

    // Expoort all results to dat file
    if(mpi.myId == 0){
        std::ofstream ofsCell(conf.outputDir + "/cell.dat");
        std::ofstream ofsNode(conf.outputDir + "/node.dat");
        std::ofstream ofsVelDirichlet(conf.outputDir + "/velocityDirichlet.dat");
        std::ofstream ofsPreDirichlet(conf.outputDir + "/pressureDirichlet.dat");
        std::ofstream ofsControlBoundary(conf.outputDir + "/controlBoundary.dat");
        std::ofstream ofsControlNodeInCell(conf.outputDir + "/controlNodeInCell.dat");
        std::ofstream ofsControlCellMap(conf.outputDir + "/controlCellMap.dat");
    
        // Cell dat
        for(int ic=0; ic<conf.nCellsGlobal; ic++){
            for(int p=0; p<conf.nNodesInCell; p++){
                ofsCell << grid.cell(ic).node[p] << " ";
            }
            ofsCell << std::endl;
        }
        ofsCell.close();

        for(int ic=0; ic<controlNodeInCell.size(); ic++){
            for(int p=0; p<4; p++){
                ofsControlNodeInCell << controlNodeInCell[ic][p] << " ";
            }
            ofsControlNodeInCell << std::endl;
        }
        ofsControlNodeInCell.close();

        for(int ic=0; ic<controlCellMap.size(); ic++){
            ofsControlCellMap << controlCellMap[ic] << std::endl;
        }
        ofsControlCellMap.close();
    
        // Node coordinates dat
        for(int in=0; in<conf.nNodesGlobal; in++){
            for(int d=0; d<conf.dim; d++){
                ofsNode << grid.node.x[in][d] << " ";
            }
            ofsNode << std::endl;
        }
        ofsNode.close();

        // Dirichlet dat
        for(auto vec : bdFaces){
            for(int in=0; in<vec.node.size(); in++){
                if(vec.dirichletType[in] == "v"){
                    ofsVelDirichlet << vec.node.at(in) << " ";
                    for(int d=0; d<vec.dirichletValue[in].size(); d++)
                        ofsVelDirichlet << vec.dirichletValue[in][d] << " ";
                    ofsVelDirichlet << std::endl;
                }
                else if(vec.dirichletType.at(in) == "p"){
                    ofsPreDirichlet << vec.node[in] << " ";
                    ofsPreDirichlet << vec.dirichletValue[in][0] << std::endl;
                }
            }
        }
        ofsVelDirichlet.close();
        ofsPreDirichlet.close();

        for(int ib=0; ib<controlBoundaryMap.size(); ib++){
            ofsControlBoundary << controlBoundaryMap[ib] << std::endl;
        }
        ofsControlBoundary.close();
    }

    return EXIT_SUCCESS;
}
