/**
 * @file main.cpp
 * @author K.Ueda
 * @date July, 2024
 */

#include <unistd.h>
#include "InverseProblem.h"
#include "MyMPI.h"
MyMPI mpi;

/**********************************
 * @brief Main function.
 * @param argc Number of arguments.
 * @param argv Arguments.
 * @return EXIT_SUCCESS.
 */
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
  std::string appName = "FDVAR";

  std::unique_ptr<Config> conf;
  conf.reset(new Config(input, appName));

  if (conf->isReadingError)
  {
    std::cerr << "Reading error." << std::endl;
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  
  InverseProblem inverse(*conf);
  inverse.initialize(*conf);

  conf.reset();

  try
  {
    inverse.runSimulation();
  }
  catch (const std::runtime_error &e)
  {
    if (mpi.myId == 0)
    {
      std::cerr << "Exception : " << e.what() << std::endl;
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
  std::string appName = "FDVAR";

  auto conf = std::make_unique<Config>(input, appName);
  if (conf->isReadingError)
    return EXIT_FAILURE;

  conf->velocityData.resize(conf->nSnapShot);

  // Read data from bin
  for (int step = 0; step < conf->nSnapShot; step++)
  {
    std::string velFile = conf->inputDir + "/data_" + to_string(step) + ".bin";
    try
    {
      BIN::importVectorDataBIN(velFile, conf->velocityData[step]);
    }
    catch (const std::runtime_error &e)
    {
      return EXIT_FAILURE;
    }
  }

  if (conf->gridType == GridType::STRUCTURED)
  {
    conf->setSolidBoundary();
    if (conf->extractFluid == ON)
    {
      conf->setFluidDomain();
    }
  }

  InverseProblem inverse(*conf);
  inverse.initialize(*conf);

  conf.reset();
  inverse.runSimulation();

  PetscPrintf(MPI_COMM_WORLD, "\nTerminated.\n");
  PetscFinalize();

  return EXIT_SUCCESS;
}
*/