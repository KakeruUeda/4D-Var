/**
 * @file ImportDAT.inl
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef IMPORT_DAT_INL_H
#define IMPORT_DAT_INL_H

#include "Import.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

template <typename T> 
void IMPORT::importScalarDataDAT(const std::string &file, std::vector<T> &vec)
{
  std::ifstream ifs(file);
  if(!ifs) {
    throw std::runtime_error("Couldn't open file: " + file);
  }

  vec.clear();
  T value;
  while(ifs >> value) {
    vec.push_back(value);
  }

  if(ifs.bad()) {
    throw std::runtime_error("Error occurred while reading file: " + file);
  }

  ifs.close();
}

template <typename T> 
void IMPORT::importVectorDataDAT(const std::string &file, std::vector<std::vector<T>> &vec)
{
  std::ifstream ifs(file);
  if(!ifs) {
    throw std::runtime_error("Couldn't open : " + file);
  }

  vec.clear();
  std::string line;

  while(std::getline(ifs, line)) {
    std::istringstream rowStream(line);
    std::vector<T> row;
    T value;
    while(rowStream >> value) {
      row.push_back(value);
    }
    vec.push_back(row);
  }

  if(ifs.bad()) {
    throw std::runtime_error("Error while reading: " + file);
  }

  ifs.close();
}

template <typename K, typename V>
void IMPORT::importScalarMapDataDAT(const std::string &file, std::vector<K> &keys, std::vector<V> &values)
{
  std::ifstream ifs(file);
  if(!ifs) {
    throw std::runtime_error("Couldn't open file: " + file);
  }

  keys.clear();
  values.clear();
  K key;
  V value;
  while(ifs >> key >> value) {
    keys.push_back(key);
    values.push_back(value);
  }

  if(ifs.bad()) {
    throw std::runtime_error("Error occurred while reading file: " + file);
  }

  ifs.close();
}

template <typename K, typename V>
void IMPORT::importVectorMapDataDAT(const std::string &file, std::vector<K> &keys, std::vector<std::vector<V>> &values)
{
  std::ifstream ifs(file);
  if(!ifs) {
    throw std::runtime_error("Couldn't open file: " + file);
  }

  keys.clear();
  values.clear();
  K key;
  V value;
  std::string line;

  while(std::getline(ifs, line)) {
    std::istringstream iss(line);
    if(!(iss >> key)) {
      throw std::runtime_error("Error reading key from line: " + line);
    }
    keys.push_back(key);

    std::vector<V> value_vector;
    while(iss >> value) {
      value_vector.push_back(value);
    }
    values.push_back(value_vector);
  }

  if(ifs.bad()) {
    throw std::runtime_error("Error occurred while reading file: " + file);
  }

  ifs.close();
}

#endif  // IMPORT_DAT_INL_H
