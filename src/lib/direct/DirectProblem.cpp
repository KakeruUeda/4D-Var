#include "DirectProblem.h"

DirectProblem::DirectProblem(Config &conf) : 
dim(conf.dim), outputDir(conf.outputDir), 
nOMP(conf.nOMP), grid(conf)
{
    std::string dir;
    std::string output = "output";
    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    outputDir = "output/" + outputDir;
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/velocity";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/pressure";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/domain";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/dat";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    grid.dirichlet.initialize(conf);
    grid.cell.initialize(conf);
    grid.node.initialize(conf);
}


void DirectProblem::visualizeDomain()
{
    if(mpi.myId > 0) return;

    std::string vtuFile;
    vtuFile = outputDir + "/domain/meshPartition.vtu";
    grid.outputVTU.exportMeshPartitionVTU(vtuFile, grid.node, grid.cell);

    vtuFile = outputDir + "/domain/phi.vtu";
    grid.outputVTU.exportPhiVTU(vtuFile, grid.node, grid.cell);
}

void DirectProblem::runSimulation()
{   
    preprocess();
    visualizeDomain();

    return;
}