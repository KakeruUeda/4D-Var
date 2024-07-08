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

void DirectProblem::initialize(Config &conf)
{
    nu = mu / rho;
    Re = 1e0 / nu; 

    VecTool::resize(grid.node.v, grid.node.nNodesGlobal, dim);
    VecTool::resize(grid.node.vPrev, grid.node.nNodesGlobal, dim);
    //VecTool::resize(grid.node.vt, timeMax, grid.node.nNodesGlobal, dim);
    VecTool::resize(grid.node.p, grid.node.nNodesGlobal);
    //VecTool::resize(grid.node.pt, timeMax, grid.node.nNodesGlobal);
    VecTool::resize(snap.v, snap.nSnapShot, grid.nNodesGlobal, dim);

    grid.dirichlet.initialize(conf);
    grid.cell.initialize(conf);
    grid.node.initialize(conf);
   
    grid.prepareMatrix(petsc, outputDir, timeMax);

    VecTool::resize(petsc.solution, grid.nDofsGlobal);
    VecTool::resize(grid.dirichlet.dirichletBCsValue, grid.nDofsGlobal);
    VecTool::resize(grid.dirichlet.dirichletBCsValueNew, grid.nDofsGlobal);
    VecTool::resize(grid.dirichlet.dirichletBCsValueInit, grid.nDofsGlobal);
    VecTool::resize(grid.dirichlet.dirichletBCsValueNewInit, grid.nDofsGlobal);
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