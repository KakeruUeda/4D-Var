/**
 * @file main.cpp
 * @brief create Structured Grid
 * @author K.U
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
#include "OutputTools.h"
#include "Grid.h"
MyMPI mpi;


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

    grid.node.resize(conf.nNodesGlobal);
    for(int in=0; in<conf.nNodesGlobal; in++)
        grid.node(in).x.resize(conf.dim);

    // Set structured grid components (cell, coordinate informateion)
    grid.setStructuredGrid(conf.nxCells, conf.nyCells, conf.nzCells, conf.nxNodes, 
                           conf.nyNodes, conf.nzNodes, conf.dx, conf.dy, conf.dz, 
                           conf.nNodesInCell, conf.dim, grid.cell, grid.node);

     // Define structured boundary faces
     std::vector<StructuredBoundaryFace> bdFaces;
     bdFaces.reserve(conf.bdStr.size());

    for(auto str : conf.bdStr)
        bdFaces.emplace_back(str);
    for(auto &vec : bdFaces)
        vec.setNodesOnBoundaryFace(conf.nxNodes, conf.nyNodes, conf.nzNodes);

    size_t iteration = 0;

    // Insert dirichlet info from config parameter
    for(auto &vec : bdFaces){
        vec.setDirichletInfo(conf.bdType, conf.bdValue, conf.dim, iteration);
        iteration++;
    }
    
    // Expoort all results to dat file
    if(mpi.myId == 0){
        std::ofstream ofsCell(conf.outputDir + "/cell.dat");
        std::ofstream ofsNode(conf.outputDir + "/node.dat");
        std::ofstream ofsVelDirichlet(conf.outputDir + "/veocityDirichlet.dat");
        std::ofstream ofsPreDirichlet(conf.outputDir + "/pressureDirichlet.dat");
    
        // Cell dat
        for(size_t ic=0; ic<conf.nCellsGlobal; ic++){
            for(size_t p=0; p<conf.nNodesInCell; p++){
                ofsCell << grid.cell(ic).node(p) << " ";
            }
            ofsCell << std::endl;
        }
        ofsCell.close();
    
        // Node coordinates dat
        for(size_t in=0; in<conf.nCellsGlobal; in++){
            for(size_t p=0; p<conf.dim; p++){
                ofsNode << grid.node(in).x(p) << " ";
            }
            ofsNode << std::endl;
        }
        ofsNode.close();

        // Dirichlet dat
        for(auto vec : bdFaces){
            for(size_t in=0; in<vec.node.size(); in++){
                if(vec.dirichletType.at(in) == "v"){
                    ofsVelDirichlet << vec.node.at(in) << " ";
                    for(int d=0; d<vec.dirichletValue.at(in).size(); d++)
                        ofsVelDirichlet << vec.dirichletValue.at(in).at(d) << " ";
                    ofsVelDirichlet << std::endl;
                }
                else if(vec.dirichletType.at(in)== "p"){
                    ofsPreDirichlet << vec.node.at(in) << " ";
                    ofsPreDirichlet << vec.dirichletValue.at(in).at(0) << std::endl;
                }
            }
        }
        ofsVelDirichlet.close();
        ofsPreDirichlet.close();
    }

    return EXIT_SUCCESS;
}
