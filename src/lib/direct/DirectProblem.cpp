/**
 * @file DirectProblem.cpp
 * @author k.ueda
 * @date July, 2024
 */

#include "DirectProblem.h"

DirectProblem::DirectProblem(Config &conf) : 
app(conf.app), dim(conf.dim), outputDir(conf.outputDir), 
nOMP(conf.nOMP), grid(conf), snap(conf),
rho(conf.rho), mu(conf.mu), dt(conf.dt), NRtolerance(conf.NRtolerance),
timeMax(conf.timeMax), pulsatileFlow(conf.pulsatileFlow), 
pulseBeginItr(conf.pulseBeginItr), T(conf.T),
alpha(conf.alpha), resistance(conf.resistance)  
{
    std::string dir;
    std::string output = "output";
    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    outputDir = "output/" + outputDir;
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    outputDir = outputDir + "/main";
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/solution";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/domain";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/data";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/dat";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}


void DirectProblem::outputDomain()
{
    if(mpi.myId > 0) return;
    std::string vtuFile;

    vtuFile = outputDir + "/domain/meshPartition.vtu";
    grid.output.exportMeshPartitionVTU(vtuFile, grid.node, grid.cell);
    vtuFile = outputDir + "/domain/phi.vtu";
    grid.output.exportPhiVTU(vtuFile, grid.node, grid.cell);
}


void DirectProblem::runSimulation()
{
    outputDomain();
    solveUSNS(app);
}

void DirectProblem::updateRowIndex(const int ii, const int ic)
{
    IU = grid.cell(ic).dofStart[ii]; IV = IU + 1;  IW = IU + 2;
    IP = IU + 3; ILU = IU + 4; ILV = IU + 5; ILW = IU + 6;
}

void DirectProblem::updateColumnIndex(const int jj, const int ic)
{
    JU = grid.cell(ic).dofStart[jj]; JV = JU + 1;  JW = JU + 2;
    JP = JU + 3; JLU = JU + 4; JLV = JU + 5; JLW = JU + 6;
}
