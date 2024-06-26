#include <unistd.h>
#include "InverseProblem.h"
#include "MyMPI.h"
MyMPI mpi;

int main(int argc, char *argv[])
{
    std::string petscfile = argv[2];
    PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);

    mpi.setSizeAndRank();
    mpi.printSizeAndRank();

    std::string inputFile = argv[1];
    std::string appName = "FDVAR";

    // Configurate using Text Parser
    Config* conf = new Config(inputFile, appName);
    if(conf->isReadingError) return EXIT_FAILURE;

    if(conf->gridType == GridType::STRUCTURED)
        conf->setSolidBoundary();

    if(conf->extractFluid == ON)
        conf->setFluidDomain();


    InverseProblem inverse(*conf);
    inverse.initialize(*conf);
    delete conf;
    
    inverse.runSimulation();

    PetscFinalize(); 

    return EXIT_SUCCESS;
}
