/**
 * @file Config.cpp
 * @author K.Ueda
 * @date May, 2024
 */

#include "Config.h"

/**
 * @brief Construct config.
 */
Config::Config(std::string inputFile, std::string appName)
{
  setApplication(appName);
  tryOpenConfigFile(inputFile);
  tryReadConfigFile();
}

/**
 * @brief Set application.
 */
void Config::setApplication(std::string appName)
{
  if(appName == "SNS") {
    app = Application::SNS;
  } else if(appName == "USNS") {
    app = Application::USNS;
  } else if(appName == "VOXELDATACREATION") {
    app = Application::VOXELDATACREATION;
  } else if(appName == "FDVAR") {
    app = Application::FDVAR;
  } else if(appName == "POSTFDVAR") {
    app = Application::POSTFDVAR;
  } else if(appName == "GRIDCREATION") {
    app = Application::GRIDCREATION;
  } else if(mpi.myId == 0) {
    std::cout << "Unknown appName" << std::endl;
  }
}

/**
 * @brief Open config file.
 */
void Config::tryOpenConfigFile(std::string inputFile)
{
  try {
    int error;
    if((error = tp.read(inputFile)) != TP_NO_ERROR) {
      throw std::runtime_error("Open error");
    }
  } catch(const std::runtime_error &e) {
    if(mpi.myId == 0) {
      std::cout << e.what() << std::endl;
    }
    isReadingError = true;
  }
}

/**
 * @brief Read config file.
 */
void Config::tryReadConfigFile()
{
  try {
    readConfigFile();
  } catch(const std::runtime_error &e) {
    if(mpi.myId == 0) {
      std::cout << e.what() << std::endl;
    }
    isReadingError = true;
  }
}

/**
 * @brief Read config file.
 */
void Config::readConfigFile()
{
  switch(app) {
  case Application::GRIDCREATION: {
    TextReaderGridCreation reader;
    reader.readBasicInfo(*this);
    reader.readGridInfo(*this);
    reader.readStructuredBoundaryInfo(*this);
    break;
  }
  case Application::USNS: {
    TextReaderUSNS reader;
    reader.readBasicInfo(*this);
    reader.readGridInfo(*this);
    reader.readBoundaryInfo(*this);
    reader.readPhysicalInfo(*this);
    reader.readDarcyInfo(*this);
    reader.readTimeInfo(*this);
    break;
  }
  case Application::VOXELDATACREATION: {
    TextReaderVoxelDataCreation reader;
    reader.readBasicInfo(*this);
    reader.readGridInfo(*this);
    reader.readSnapInfo(*this);
    reader.readOriginalInfo(*this);
    break;
  }
  case Application::FDVAR: {
    TextReader4DVar reader;
    reader.readBasicInfo(*this);
    reader.readGridInfo(*this);
    reader.readBoundaryInfo(*this);
    reader.readPhysicalInfo(*this);
    reader.readDarcyInfo(*this);
    reader.readTimeInfo(*this);
    reader.readInverseInfo(*this);
    reader.readDataInfo(*this);
    break;
  }
  default:
    throw std::runtime_error("Unknown Application");
    break;
  }
}
