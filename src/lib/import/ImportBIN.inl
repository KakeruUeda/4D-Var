/**
 * @file ImportBIN.inl
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef IMPORT_BIN_INL_H
#define IMPORT_BIN_INL_H

#include "Import.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

template <typename T> 
void IMPORT::importScalarDataBIN(const std::string &file, std::vector<T> &vec)
{
  std::ifstream inFile(file, std::ios::binary);

  if(!inFile) {
    std::cerr << "Error opening file: " << file << std::endl;
    return;
  }

  size_t dataSize;
  inFile.read(reinterpret_cast<char *>(&dataSize), sizeof(size_t));

  vec.resize(dataSize);
  inFile.read(reinterpret_cast<char *>(vec.data()), dataSize * sizeof(T));

  if(!inFile.good()) {
    std::cerr << "Error occurred while reading from file: " << file << std::endl;
    vec.clear();
  }

  inFile.close();
}

template <typename T> 
void IMPORT::importVectorDataBIN(const std::string &file, std::vector<std::vector<T>> &vec)
{
  std::ifstream inFile(file, std::ios::binary);

  if(!inFile) {
    std::cerr << "Error opening file: " << file << std::endl;
    return;
  }

  size_t outerSize;
  inFile.read(reinterpret_cast<char *>(&outerSize), sizeof(size_t));

  vec.resize(outerSize);

  for(size_t i = 0; i < outerSize; ++i) {
    size_t innerSize;
    inFile.read(reinterpret_cast<char *>(&innerSize), sizeof(size_t));

    if(innerSize > 0) {
      vec[i].resize(innerSize);
      inFile.read(reinterpret_cast<char *>(vec[i].data()), innerSize * sizeof(T));
    }
  }

  if(!inFile.good()) {
    std::cerr << "Error occurred while reading from file: " << file << std::endl;
    vec.clear();
  }

  inFile.close();
}

#endif  // IMPORT_BIN_INL_H