/*=================================
  Execute with single MPI process.
===================================*/
/**
 * @file main.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "DirectProblem.h"
#include "MyMPI.h"
#include "VoxelDataCreation.h"
#include <unistd.h>
MyMPI mpi;

int main(int argc, char *argv[])
{
  if(argc < 2) {
    if(mpi.myId == 0) {
      std::cerr << "argc error" << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  std::string inputFile = argv[1];
  std::string appName = "VOXELDATACREATION";

  Config conf(inputFile, appName);
  if(conf.isReadingError) {
    std::cerr << "Error reading configuration file." << std::endl;
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  VoxelDataCreation vdc(conf);

  try {
    vdc.Initialize(conf);
    vdc.createRefAndData();
  } catch(const std::runtime_error &e) {
    std::cerr << "Exception : " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  vdc.outputVTK();
  vdc.outputBIN();

  return EXIT_SUCCESS;
}
