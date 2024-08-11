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

  if (argc < 2)
  {
    if (mpi.myId == 0)
    {
      std::cerr << "argc error" << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  std::string input = argv[1];
  std::string appName = "USNS";

  std::unique_ptr<Config> conf;
  conf.reset(new Config(input, appName));

  if (conf->isReadingError)
  {
    std::cerr << "Error reading configuration file." << std::endl;
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  DirectProblem direct(*conf);
  direct.initialize(*conf);
  conf.reset();

  try
  {
    //direct.resize();
    direct.runSimulation();
  }
  catch (const std::runtime_error &e)
  {
    if (mpi.myId == 0)
    {
      std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  PetscPrintf(MPI_COMM_WORLD, "\nTerminated.\n");
  PetscFinalize();

  return EXIT_SUCCESS;
}

/*
int main(int argc, char *argv[])
{
  std::string petscfile = argv[2];
  PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);

  mpi.setSizeAndRank();
  mpi.printSizeAndRank();

  std::string input = argv[1];
  std::string appName = "USNS";
  auto conf = std::make_unique<Config>(input, appName);

  if (conf->isReadingError)
    return EXIT_FAILURE;

  if (conf->gridType == GridType::STRUCTURED)
  {
    conf->setSolidBoundary();
    if (conf->extractFluid == ON)
    {
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
*/
