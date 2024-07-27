/**
 * @file DirectProblem.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "DirectProblem.h"

/***************************************************
 * @brief Construct direct problem from config file.
 */
DirectProblem::DirectProblem(Config &conf) : 
app(conf.app), dim(conf.dim), outputDir(conf.outputDir),
nOMP(conf.nOMP), grid(conf), snap(conf),
rho(conf.rho), mu(conf.mu), dt(conf.dt), NRtolerance(conf.NRtolerance),
timeMax(conf.timeMax), pulsatileFlow(conf.pulsatileFlow), 
pulseBeginItr(conf.pulseBeginItr), T(conf.T),
alpha(conf.alpha), resistance(conf.resistance)  
{
    if(app != Application::FDVAR){
        std::string dir;
        std::string output = "output";
        mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        outputDir = "output/" + outputDir;
        mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        dir = outputDir + "/domain";
        mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        dir = outputDir + "/data";
        mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        dir = outputDir + "/solution";
        mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        dir = outputDir + "/dat";
        mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    }
}

/**********************************************
 * @brief Visualize partitioned domain and phi.
 */
void DirectProblem::outputDomain()
{
    if(mpi.myId > 0) return;
    std::string vtuFile;

    std::string vtiFile;

    vtuFile = outputDir + "/domain/meshPartition.vtu";
    VTK::exportMeshPartitionVTU(vtuFile, grid.node, grid.cell);
    vtuFile = outputDir + "/domain/phi.vtu";
    VTK::exportPhiVTU(vtuFile, grid.node, grid.cell);

    /*
    std::vector<double> vtiNodeScalar;
    vtiNodeScalar.resize(grid.node.nNodesStructuredGlobal, 1e0);
    std::vector<std::vector<double>> vtiNodeVector;
    vtiNodeVector.resize(grid.node.nNodesStructuredGlobal, std::vector<double>(3, 2e0));
    std::vector<double> vtiCellScalar;
    vtiCellScalar.resize(grid.cell.nCellsStructuredGlobal, 3e0);
    std::vector<std::vector<double>> vtiCellVector;
    vtiCellVector.resize(grid.cell.nCellsStructuredGlobal, std::vector<double>(3, 4e0));
    std::vector<double> vtuNodeScalar;
    vtuNodeScalar.resize(grid.node.nNodesGlobal, 5e0);
    std::vector<std::vector<double>> vtuNodeVector;
    vtuNodeVector.resize(grid.node.nNodesGlobal, std::vector<double>(3, 6e0));
    std::vector<double> vtuCellScalar;
    vtuCellScalar.resize(grid.cell.nCellsGlobal, 7e0);
    std::vector<std::vector<double>>  vtuCellVector;
    vtuCellVector.resize(grid.cell.nCellsGlobal, std::vector<double>(3, 8e0));

    vtiFile = outputDir + "/domain/vtiNodeScalar.vti";
    VTK::exportScalarPointDataVTI(vtiFile, "vtiNodeScalar", vtiNodeScalar, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
    vtiFile = outputDir + "/domain/vtiNodeVector.vti";
    VTK::exportVectorPointDataVTI(vtiFile, "vtiNodeVector", vtiNodeVector, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
    vtiFile = outputDir + "/domain/vtiCellScalar.vti";
    VTK::exportScalarCellDataVTI(vtiFile, "vtiCellScalar", vtiCellScalar, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
    vtiFile = outputDir + "/domain/vtiCellVector.vti";
    VTK::exportVectorCellDataVTI(vtiFile, "vtiCellVector", vtiCellVector, grid.nx, grid.ny, grid.nz, grid.dx, grid.dy, grid.dz);
    vtuFile = outputDir + "/domain/vtuNodeScalar.vtu";
    VTK::exportScalarPointDataVTU(vtuFile, "vtuNodeScalar", grid.node, grid.cell, vtuNodeScalar);
    vtuFile = outputDir + "/domain/vtuNodeVector.vtu";
    VTK::exportVectorPointDataVTU(vtuFile, "vtuNodeVector", grid.node, grid.cell, vtuNodeVector);
    vtuFile = outputDir + "/domain/vtuCellScalar.vtu";
    VTK::exportScalarCellDataVTU(vtuFile, "vtuCellScalar", grid.node, grid.cell, vtuCellScalar);
    vtuFile = outputDir + "/domain/vtuCellVector.vtu";
    VTK::exportVectorCellDataVTU(vtuFile, "vtuCellVector", grid.node, grid.cell, vtuCellVector);
    */
}

/**************************************************
 * @brief Simulate Unsteady Navier Stokes Equation.
 */
void DirectProblem::runSimulation()
{
    outputDomain();
    solveUSNS(app);
}

/****************************************************************
 * @brief Update row index when creating element stifness matrix.
 */
void DirectProblem::updateRowIndex(const int ii, const int ic)
{
    IU = grid.cell(ic).dofStart[ii]; IV = IU + 1;  IW = IU + 2;
    IP = IU + 3; ILU = IU + 4; ILV = IU + 5; ILW = IU + 6;
}

/******************************************************************
 * @brief Update column index when creating element stifness matrix.
 */
void DirectProblem::updateColumnIndex(const int jj, const int ic)
{
    JU = grid.cell(ic).dofStart[jj]; JV = JU + 1;  JW = JU + 2;
    JP = JU + 3; JLU = JU + 4; JLV = JU + 5; JLW = JU + 6;
}
