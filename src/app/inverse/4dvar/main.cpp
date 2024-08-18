/**
 * @file main.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "InverseProblem.h"
#include "MyMPI.h"
#include <unistd.h>
MyMPI mpi;


/// @brief Main function.
/// @param argc 
/// @param argv Text input file name
/// @return 

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
  std::string appName = "FDVAR";

  std::unique_ptr<Config> conf;
  conf.reset(new Config(input, appName));

  if(conf->isReadingError) {
    std::cerr << "Reading error." << std::endl;
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  InverseProblem inverse(*conf);
  inverse.initialize(*conf);

  conf.reset();

  try {
    inverse.runSimulation();
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
