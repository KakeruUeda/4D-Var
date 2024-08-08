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
    if(app == Application::USNS){
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
        dir = outputDir + "/input";
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
