/**
 * @file main.cpp
 * @brief create Structured Grid
 * @author K.Ueda
 * @date Mar 24, 2024
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

void setControlBoundary(ControlBoundary &controlBoundary, std::vector<int> &controlBoundaryMap,
                        int nx, int ny, int nz)
{
    if(controlBoundary == ControlBoundary::left){
        for(int k=0; k<nz; k++){
            for(int j=0; j<ny; j++){
                for(int i=0; i<nx; i++){
                    if(i == 0){
                        controlBoundaryMap.push_back(i + j * nx + k * nx * ny);
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
    setControlBoundary(controlBoundary, controlBoundaryMap, conf.nxNodes, conf.nyNodes, conf.nzNodes);

    // Expoort all results to dat file
    if(mpi.myId == 0){
        std::ofstream ofsCell(conf.outputDir + "/cell.dat");
        std::ofstream ofsNode(conf.outputDir + "/node.dat");
        std::ofstream ofsVelDirichlet(conf.outputDir + "/velocityDirichlet.dat");
        std::ofstream ofsPreDirichlet(conf.outputDir + "/pressureDirichlet.dat");
        std::ofstream ofsControlBoundary(conf.outputDir + "/controlBoundary.dat");
    
        // Cell dat
        for(int ic=0; ic<conf.nCellsGlobal; ic++){
            for(int p=0; p<conf.nNodesInCell; p++){
                ofsCell << grid.cell(ic).node[p] << " ";
            }
            ofsCell << std::endl;
        }
        ofsCell.close();
    
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
