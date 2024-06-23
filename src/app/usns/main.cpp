#include <unistd.h>
#include "DirectProblem.h"
#include "Postprocess.h"
#include "MyMPI.h"
MyMPI mpi;

int main(int argc, char *argv[])
{
    std::string petscfile = argv[2];
    MPI_Init(NULL, NULL);
    PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);

    mpi.setSizeAndRank();
    mpi.printSizeAndRank();

    std::string inputFile = argv[1];
    std::string appName = "USNS";

    // Configurate using Text Parser
    Config* conf = new Config(inputFile, appName);
    if(conf->isReadingError) 
        return EXIT_FAILURE;

    if(conf->gridType == GridType::STRUCTURED)
        conf->setSolidBoundary();

    if(conf->extractFluid == ON)
        conf->setFluidDomain();

    DirectProblem direct(*conf);

    // Solve Unstready Navier Stokes
    direct.runSimulation();

    Postprocess post;
    //post.extractOutletVelocity(direct);
    post.makeObservedData(direct);

    PetscFinalize(); 
    MPI_Finalize();

    return EXIT_SUCCESS;
}
