/**
 * @file main.cpp
 * @author K.Ueda
 * @date Jun, 2024
*/

#include <unistd.h>
#include "DirectProblem.h"
#include "MyMPI.h"
MyMPI mpi;

int main(int argc, char *argv[])
{
    std::string petscfile = argv[2];
    PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);

    mpi.setSizeAndRank();
    mpi.printSizeAndRank();

    std::string input = argv[1];
    std::string appName = "USNS";

    auto conf = std::make_unique<Config>(input, appName);

    if(conf->isReadingError) 
        return EXIT_FAILURE;
    
    if(conf->gridType == GridType::STRUCTURED){
        conf->setSolidBoundary();
        if(conf->extractFluid == ON){
            conf->setFluidDomain();
        }
    } 
    
    DirectProblem direct(*conf);
    direct.initialize(*conf);

    conf.reset();
    direct.runSimulation();
    
    PetscPrintf(MPI_COMM_WORLD, "\nTerminated.\n");
    PetscFinalize(); 

    return EXIT_SUCCESS;
}
