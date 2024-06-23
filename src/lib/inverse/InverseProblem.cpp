#include "InverseProblem.h"

InverseProblem::InverseProblem(Config &conf)
{
    /*
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

    grid.node.v.resize(grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    grid.node.vPrev.resize(grid.node.nNodesGlobal, std::vector<double>(dim, 0e0));
    grid.node.p.resize(grid.node.nNodesGlobal, 0e0);

    grid.dirichlet.initialize(conf);
    grid.cell.initialize(conf);
    grid.node.initialize(conf);
    */
}

void InverseProblem::runSimulation()
{
    prepareMatrix();
}

void InverseProblem::prepareMatrix()
{

}