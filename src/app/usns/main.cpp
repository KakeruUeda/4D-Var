#include <unistd.h>
#include "DirectProblem.h"
#include "MyMPI.h"
MyMPI mpi;

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);
    mpi.setSizeAndRank();
    mpi.printSizeAndRank();

    std::string inputFile = argv[1];
    std::string appName = "USNS";

    // Configurate using Text Parser
    Config* conf = new Config(inputFile, appName);
    if(conf->isReadingError) return EXIT_FAILURE;

    DirectProblem direct(*conf);
    delete conf;

    // Solve Unstready Navier Stokes
    direct.runSimulation();

    return EXIT_SUCCESS;
}
