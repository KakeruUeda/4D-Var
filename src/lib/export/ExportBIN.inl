/**
 * @file ExportBIN.inl
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef EXPORT_BIN_INL_H
#define EXPORT_BIN_INL_H

#include "Export.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

template <typename T> 
void EXPORT::exportScalarDataBIN(const std::string &file, const std::vector<T> &vec)
{
  std::ofstream outFile(file, std::ios::binary);

  if(!outFile) {
    std::cerr << "Error opening file: " << file << std::endl;
    return;
  }

  size_t dataSize = vec.size();

  outFile.write(reinterpret_cast<const char *>(&dataSize), sizeof(size_t));
  outFile.write(reinterpret_cast<const char *>(vec.data()), dataSize * sizeof(T));

  outFile.close();

  if(!outFile.good()) {
    std::cerr << "Error occurred while writing to file: " << file << std::endl;
  }
}

template <typename T> 
void EXPORT::exportVectorDataBIN(const std::string &file, const std::vector<std::vector<T>> &vec)
{
  std::ofstream outFile(file, std::ios::binary);

  if(!outFile) {
    std::cerr << "Error opening file: " << file << std::endl;
    return;
  }

  size_t outerSize = vec.size();
  outFile.write(reinterpret_cast<const char *>(&outerSize), sizeof(size_t));

  for(const auto &innerVec : vec) {
    size_t innerSize = innerVec.size();
    outFile.write(reinterpret_cast<const char *>(&innerSize), sizeof(size_t));
    if(innerSize > 0) {
      outFile.write(reinterpret_cast<const char *>(innerVec.data()), innerSize * sizeof(T));
    }
  }

  outFile.close();

  if(!outFile.good()) {
    std::cerr << "Error occurred while writing to file: " << file << std::endl;
  }
}

#endif  // EXPORT_BIN_INL_H