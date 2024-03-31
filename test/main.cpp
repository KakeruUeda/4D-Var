#include <unistd.h>
#include "DirectProblem.h"
#include "PreDirect.h"
#include "MyMPI.h"
MyMPI mpi;

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);
    mpi.setSizeAndRank();
    mpi.printSizeAndRank();

    std::string inputFile = argv[1];
    std::string appName = "USNS";
    Direct::Preprocess pre(inputFile, appName);

    if(pre.conf.isReadingError)
        return EXIT_FAILURE;

    DirectProblem direct(pre.conf);

    return EXIT_SUCCESS;
}
