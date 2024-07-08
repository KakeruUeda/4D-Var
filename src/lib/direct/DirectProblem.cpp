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