/**
 * @file ImportDAT.inl
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef IMPORT_DAT_INL_H
#define IMPORT_DAT_INL_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "Import.h"

template <typename T>
void IMPORT::importScalarDataDAT(const std::string &file, std::vector<T> &vec)
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
void IMPORT::importVectorDataDAT(const std::string &file, std::vector<std::vector<T>> &vec)
{
  std::ifstream ifs(file);
  if (!ifs)
  {
    throw std::runtime_error("Couldn't open : " + file);
  }

  vec.clear();
  std::string line;

  while (std::getline(ifs, line))
  {
    std::istringstream rowStream(line);
    std::vector<T> row;
    T value;
    while (rowStream >> value)
    {
      row.push_back(value);
    }
    vec.push_back(row);
  }

  if (ifs.bad())
  {
    throw std::runtime_error("Error while reading: " + file);
  }

  ifs.close();
}

#endif // DAT_INL_H
