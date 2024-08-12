/**
 * @file Import.h
 * @author K.Ueda
 * @date Jun, 2024
 */

#ifndef IMPORT_H
#define IMPORT_H

#include <iostream>
#include <vector>
#include "MyMPI.h"

extern MyMPI mpi;

class IMPORT
{
public:
  // DAT
  template <typename T>
  static void importScalarDataDAT(const std::string &file, std::vector<T> &vec);
  template <typename T>
  static void importVectorDataDAT(const std::string &file, std::vector<std::vector<T>> &vec);

};


#include "ImportDAT.inl"

#endif // IMPORT_H