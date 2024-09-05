/**
 * @file main.cpp
 * @author K.Ueda
 * @date Jun, 2024
 */

#include "DirectProblem.h"
#include "MyMPI.h"
MyMPI mpi;

int main(int argc, char *argv[])
{
  std::string petscfile = argv[2];
  PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);

  mpi.setSizeAndRank();
  mpi.printSizeAndRank();

  if(argc < 2) {
    if(mpi.myId == 0) {
      std::cerr << "argc error" << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  std::string input = argv[1];
  std::string appName = "USNS";

  std::unique_ptr<Config> conf(new Config(input, appName));

  if(conf->isReadingError) {
    if(mpi.myId == 0) {
      std::cerr << "Reading error." << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  DirectProblem direct(*conf);

  try {
    direct.initialize(*conf);
    direct.runSimulation();
  } catch(const std::runtime_error &e) {
    if(mpi.myId == 0) {
      std::cerr << "Exception : " << e.what() << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  PetscPrintf(MPI_COMM_WORLD, "\nTerminated.\n");
  PetscFinalize();

  return EXIT_SUCCESS;
}
