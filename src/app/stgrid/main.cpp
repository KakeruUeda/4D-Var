#include <iostream>
#include <mpi.h>
#include "TextParser.h"
#include "Config.h"
#include "MyMPI.h"
MyMPI mpi;

int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);
    mpi.setSizeAndRank();
    mpi.printSizeAndRank();

    std::string inputFile = argv[1];
    std::string appName = "STGRID";
    Config conf(inputFile, appName);

    return EXIT_SUCCESS;
}