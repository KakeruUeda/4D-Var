/**
 * @file DAT.inl
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef DAT_INL_H
#define DAT_INL_H

#include "FileIO.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

template <typename T>
void DAT::exportScalarDataDAT(const std::string &file, const std::vector<T> &vec)
{
  std::ofstream ofs(file);
  if (!ofs)
  {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  for (const T &val : vec)
  {
    ofs << val << "\n";
  }

  ofs.close();
  if (!ofs.good())
  {
    std::cerr << "Error occurred at writing time." << std::endl;
  }
}

template <typename T>
void DAT::exportVectorDataDAT(const std::string &file, const std::vector<std::vector<T>> &vec)
{
  std::ofstream ofs(file);
  if (!ofs)
  {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  ofs << vec.size() << "\n";
  for (const auto &row : vec)
  {
    ofs << row.size() << " ";
    for (const T &val : row)
    {
      ofs << val << " ";
    }
    ofs << "\n";
  }

  ofs.close();
  if (!ofs.good())
  {
    std::cerr << "Error occurred at writing time." << std::endl;
  }
}

template <typename T>
void DAT::importScalarDataDAT(const std::string &file, std::vector<T> &vec)
{
  std::ifstream ifs(file);
  if (!ifs)
  {
    throw std::runtime_error("Couldn't open file: " + file);
  }

  vec.clear();
  T value;
  while (ifs >> value)
  {
    vec.push_back(value);
  }

  if (ifs.bad())
  {
    throw std::runtime_error("Error occurred while reading file: " + file);
  }

  ifs.close();
}

template <typename T>
void DAT::importVectorDataDAT(const std::string &file, std::vector<std::vector<T>> &vec)
{
  std::ifstream ifs(file);
  if (!ifs)
  {
    throw std::runtime_error("Couldn't open file: " + file);
  }

  vec.clear();
  std::string line;
  if (std::getline(ifs, line))
  {
    std::istringstream iss(line);
    int rows;
    if (!(iss >> rows))
    {
      throw std::runtime_error("Failed to read number of rows: " + file);
    }

    vec.resize(rows);

    for (int i = 0; i < rows; ++i)
    {
      if (!std::getline(ifs, line))
      {
        throw std::runtime_error("Failed to read row data: " + file);
      }
      std::istringstream rowStream(line);
      int cols;
      if (!(rowStream >> cols))
      {
        throw std::runtime_error("Failed to read number of columns in row " + std::to_string(i) + ": " + file);
      }

      vec[i].resize(cols);
      for (int j = 0; j < cols; ++j)
      {
        if (!(rowStream >> vec[i][j]))
        {
          throw std::runtime_error("Failed to read data in row " + std::to_string(i) + ", column " + std::to_string(j) + ": " + file);
        }
      }
    }
  }

  if (ifs.bad())
  {
    throw std::runtime_error("Error occurred while reading file: " + file);
  }

  ifs.close();
}

template <typename T>
void DAT::exportMapDataDAT(const std::string &file, const std::map<int, std::vector<T>> &dataMap)
{
  std::ofstream ofs(file);
  if (!ofs)
  {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  for (const auto &entry : dataMap)
  {
    ofs << entry.first;
    for (const auto &val : entry.second)
    {
      ofs << " " << val;
    }
    ofs << "\n";
  }

  ofs.close();
  if (!ofs.good())
  {
    std::cerr << "Error occurred at writing time." << std::endl;
  }
}

template <typename T>
void DAT::exportMapDataDAT(const std::string &file, const std::map<int, T> &dataMap)
{
  std::ofstream ofs(file);
  if (!ofs)
  {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  for (const auto &entry : dataMap)
  {
    ofs << entry.first << " " << entry.second << "\n";
  }

  ofs.close();
  if (!ofs.good())
  {
    std::cerr << "Error occurred at writing time." << std::endl;
  }
}

#endif // DAT_INL_H
