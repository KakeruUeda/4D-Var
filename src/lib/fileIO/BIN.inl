/**
 * @file BIN.inl
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef BIN_INL_H
#define BIN_INL_H

#include "FileIO.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

template <typename T>
void BINTMP::exportScalarDataBIN(const std::string &file, std::vector<T> &vec)
{
  std::ofstream ofs(file, std::ios::binary);
  if (!ofs)
  {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  ofs.write(reinterpret_cast<const char *>(vec.data()), vec.size() * sizeof(T));
  ofs.close();
  if (!ofs.good())
  {
    std::cerr << "Error occurred at writing time." << std::endl;
  }
}

template <typename T>
void BINTMP::exportVectorDataBIN(const std::string &file, std::vector<std::vector<T>> &vec)
{
  std::ofstream ofs(file, std::ios::binary);
  if (!ofs)
  {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  int rows = vec.size();
  ofs.write(reinterpret_cast<const char *>(&rows), sizeof(rows));

  for (const auto &row : vec)
  {
    int cols = row.size();
    ofs.write(reinterpret_cast<const char *>(&cols), sizeof(cols));
    ofs.write(reinterpret_cast<const char *>(row.data()), cols * sizeof(T));
  }

  ofs.close();
  if (!ofs.good())
  {
    std::cerr << "Error occurred at writing time." << std::endl;
  }
}

template <typename T>
void BINTMP::importScalarDataBIN(const std::string &file, std::vector<T> &vec)
{
  std::ifstream ifs(file, std::ios::binary);
  if (!ifs)
  {
    throw std::runtime_error("Couldn't open file: " + file);
  }

  ifs.seekg(0, std::ios::end);
  std::streamsize size = ifs.tellg();
  ifs.seekg(0, std::ios::beg);

  if (size % sizeof(T) != 0)
  {
    throw std::runtime_error("File size error: " + file);
  }

  vec.resize(size / sizeof(T));

  if (!ifs.read(reinterpret_cast<char *>(vec.data()), size))
  {
    throw std::runtime_error("Couldn't read file: " + file);
  }

  ifs.close();
}

template <typename T>
void BINTMP::importVectorDataBIN(const std::string &file, std::vector<std::vector<T>> &vec)
{
  std::ifstream ifs(file, std::ios::binary);
  if (!ifs)
  {
    throw std::runtime_error("Couldn't open file: " + file);
  }

  int rows;
  if (!ifs.read(reinterpret_cast<char *>(&rows), sizeof(rows)))
  {
    throw std::runtime_error("Failed to read rows: " + file);
  }

  vec.resize(rows);

  for (int j = 0; j < rows; j++)
  {
    int cols;
    if (!ifs.read(reinterpret_cast<char *>(&cols), sizeof(cols)))
    {
      throw std::runtime_error("Failed to read cols: " + file);
    }
    vec[j].resize(cols);
    if (!ifs.read(reinterpret_cast<char *>(vec[j].data()), cols * sizeof(T)))
    {
      throw std::runtime_error("Couldn't read file: " + file);
    }
  }

  ifs.close();
}

#endif // BIN_INL_H