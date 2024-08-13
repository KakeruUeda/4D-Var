/**
 * @file Import.h
 * @author K.Ueda
 * @date Jun, 2024
 */

#ifndef IMPORT_H
#define IMPORT_H

#include <iostream>
#include <vector>
#include <set>
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
  template <typename T>
  static void importScalarDataDAT(const std::string &file, std::set<T> &set);
  template <typename K, typename V>
  static void importScalarMapDataDAT(const std::string &file, std::vector<K> &keys, std::vector<V> &values);
  template <typename K, typename V>
  static void importVectorMapDataDAT(const std::string &file, std::vector<K> &keys, std::vector<std::vector<V>> &values);
  // BIN
  template <typename T>
  static void importScalarDataBIN(const std::string &file, std::vector<T> &vec);
  template <typename T>
  static void importVectorDataBIN(const std::string &file, std::vector<std::vector<T>> &vec);
};

#include "ImportDAT.inl"
#include "ImportBIN.inl"

#endif // IMPORT_H