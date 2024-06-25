#include "InverseProblem.h"

InverseProblem::InverseProblem(Config &conf):
main(conf), adjoint(conf), data(conf),
dim(conf.dim), outputDir(conf.outputDir), 
nOMP(conf.nOMP), timeMax(conf.timeMax),
rho(conf.rho), mu(conf.mu), dt(conf.dt),
alpha(conf.alpha), resistance(conf.resistance),
aCF(conf.aCF), bCF1(conf.bCF1), bCF2(conf.bCF2), gCF(conf.gCF)
{
    std::string dir;
    std::string output = "output";
    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    outputDir = "output/" + outputDir;
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    outputDir = outputDir + "/adjoint";
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/velocity";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/pressure";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/domain";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    dir = outputDir + "/dat";
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}

void InverseProblem::initialize(Config &conf)
{
    main.grid.node.v.resize(main.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    main.grid.node.vPrev.resize(main.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    main.grid.node.vt.resize(timeMax, std::vector<std::vector<double>>(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0)));
    main.grid.node.pt.resize(timeMax, std::vector<double>(main.grid.node.nNodesGlobal, 0e0));
    
    main.grid.dirichlet.initialize(conf);
    main.grid.cell.initialize(conf);
    main.grid.node.initialize(conf);
    
    adjoint.grid.node.v.resize(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    adjoint.grid.node.vPrev.resize(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    adjoint.grid.node.vt.resize(timeMax, std::vector<std::vector<double>>(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0)));
    adjoint.grid.node.pt.resize(timeMax, std::vector<double>(adjoint.grid.node.nNodesGlobal, 0e0));
    adjoint.grid.node.lambda.resize(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    adjoint.grid.node.lambdat.resize(timeMax, std::vector<std::vector<double>>(adjoint.grid.node.nNodesGlobal, std::vector<double>(dim, 0e0)));
    
    adjoint.grid.dirichlet.initializeAdjoint(conf);
    adjoint.grid.cell.initializeAdjoint(conf);
    adjoint.grid.node.initializeAdjoint(conf, adjoint.grid.dirichlet.controlBoundaryMap);
}

void InverseProblem::runSimulation()
{
    PetscPrintf(MPI_COMM_WORLD, "\n*** Main initialize ***\n\n");
    main.grid.prepareMatrix(main.petsc, main.outputDir);

    PetscPrintf(MPI_COMM_WORLD, "\n*** Adjoint initialize ***\n\n");
    adjoint.grid.prepareMatrix(adjoint.petsc, outputDir);

    if(mpi.myId == 0){
        std::string vtuFile;
        vtuFile = outputDir + "/domain/meshPartition.vtu";
        adjoint.grid.vtu.exportMeshPartitionVTU(vtuFile, adjoint.grid.node, adjoint.grid.cell);
        vtuFile = outputDir + "/domain/phi.vtu";
        adjoint.grid.vtu.exportPhiVTU(vtuFile, adjoint.grid.node, adjoint.grid.cell);
    }

}
