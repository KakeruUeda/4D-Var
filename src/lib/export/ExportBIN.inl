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
void EXPORT::exportScalarDataBIN(const std::string &file, std::vector<T> &vec)
{
  std::ofstream ofs(file, std::ios::binary);
  if(!ofs) {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  ofs.write(reinterpret_cast<const char *>(vec.data()), vec.size() * sizeof(T));
  ofs.close();
  if(!ofs.good()) {
    std::cerr << "Error occurred at writing time." << std::endl;
  }
}

template <typename T>
void EXPORT::exportVectorDataBIN(const std::string &file, std::vector<std::vector<T>> &vec)
{
  std::ofstream ofs(file, std::ios::binary);
  if(!ofs) {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  int rows = vec.size();
  ofs.write(reinterpret_cast<const char *>(&rows), sizeof(rows));

  for(const auto &row : vec) {
    int cols = row.size();
    ofs.write(reinterpret_cast<const char *>(&cols), sizeof(cols));
    ofs.write(reinterpret_cast<const char *>(row.data()), cols * sizeof(T));
  }

  ofs.close();
  if(!ofs.good()) {
    std::cerr << "Error occurred at writing time." << std::endl;
  }
}

#endif  // EXPORT_BIN_INL_H