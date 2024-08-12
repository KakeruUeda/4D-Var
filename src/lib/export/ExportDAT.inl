/**
 * @file ExportDAT.inl
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef DAT_INL_H
#define DAT_INL_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "Export.h"

template <typename T>
void EXPORT::exportScalarDataDAT(const std::string &file, const std::vector<T> &vec)
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
void EXPORT::exportVectorDataDAT(const std::string &file, const std::vector<std::vector<T>> &vec)
{
  std::ofstream ofs(file);
  if (!ofs)
  {
    std::cerr << "Could not open " << file << std::endl;
    return;
  }

  for (const auto &row : vec)
  {
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
void EXPORT::exportMapDataDAT(const std::string &file, const std::map<int, std::vector<T>> &dataMap)
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
void EXPORT::exportMapDataDAT(const std::string &file, const std::map<int, T> &dataMap)
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
