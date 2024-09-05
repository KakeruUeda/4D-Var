/**
 * @file GridCreation.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include "GridCreation.h"
MyMPI mpi;

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  mpi.setSizeAndRank();
  mpi.printSizeAndRank();

  if(argc < 2) {
    if(mpi.myId == 0) {
      std::cerr << "argc error" << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  std::string inputFile = argv[1];
  std::string appName = "GRIDCREATION";

  Config conf(inputFile, appName);
  if(conf.isReadingError) {
    if(mpi.myId == 0) {
      std::cerr << "Error while reading." << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  conf.setSolidDirichletValue();
  conf.filterFluidGrid();

  GridCreation gc(conf);

  try {
    gc.initialize(conf);
    gc.divideWholeGrid();
    gc.collectLocalGrid();
  } catch(const std::runtime_error &e) {
    if(mpi.myId == 0) {
      std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  gc.outputDat();
  gc.outputVTU();  /* debug */

  MPI_Finalize();
  return EXIT_SUCCESS;
}