/**
 * @file   Array.h
 * @autor  K.Ueda
 * @date   May, 2024
 */

#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <memory>
#include <algorithm>

template <class T>
class Array1D
{
public:
  Array1D(int width) : width(width), data(std::make_unique<T[]>(width)) {}
  Array1D() : width(0) {}

  inline T &operator()(int x)
  {
    return data[x];
  }

  inline int size() const
  {
    return width;
  }

  inline void resize(int newWidth)
  {
    data = std::make_unique<T[]>(newWidth);
    width = newWidth;
  }

  inline void zeroFill()
  {
    std::fill(data.get(), data.get() + width, 0);
  }

private:
  int width;
  std::unique_ptr<T[]> data;
};

template <class T>
class Array2D
{
public:
  Array2D(int height, int width) : width(width), height(height), data(std::make_unique<T[]>(width * height)) {}
  Array2D() : width(0), height(0) {}

  inline T &operator()(int y, int x)
  {
    return data[y * width + x];
  }

  inline int size() const
  {
    return width * height;
  }

  inline void resize(int newHeight, int newWidth)
  {
    data = std::make_unique<T[]>(newHeight * newWidth);
    height = newHeight;
    width = newWidth;
  }

  inline void zeroFill()
  {
    std::fill(data.get(), data.get() + height * width, 0);
  }

private:
  int width, height;
  std::unique_ptr<T[]> data;
};

template <class T>
class Array3D
{
public:
  Array3D(int depth, int height, int width) : width(width), height(height), depth(depth), data(std::make_unique<T[]>(depth * height * width)) {}
  Array3D() : width(0), height(0), depth(0) {}

  inline T &operator()(int z, int y, int x)
  {
    return data[z * height * width + y * width + x];
  }

  inline int size() const
  {
    return depth * height * width;
  }

  inline void resize(int newDepth, int newHeight, int newWidth)
  {
    data = std::make_unique<T[]>(newDepth * newHeight * newWidth);
    depth = newDepth;
    height = newHeight;
    width = newWidth;
  }

  inline void zeroFill()
  {
    std::fill(data.get(), data.get() + depth * height * width, 0);
  }

private:
  int width, height, depth;
  std::unique_ptr<T[]> data;
};

template <class T>
class Array4D
{
public:
  Array4D(int time, int depth, int height, int width) : time(time), depth(depth), height(height), width(width), data(std::make_unique<T[]>(time * depth * height * width)) {}
  Array4D() : time(0), depth(0), height(0), width(0) {}

  inline T &operator()(int t, int z, int y, int x)
  {
    return data[t * depth * height * width + z * height * width + y * width + x];
  }

  inline int size() const
  {
    return time * depth * height * width;
  }

  inline void resize(int newTime, int newDepth, int newHeight, int newWidth)
  {
    data = std::make_unique<T[]>(newTime * newDepth * newHeight * newWidth);
    time = newTime;
    depth = newDepth;
    height = newHeight;
    width = newWidth;
  }

  inline void zeroFill()
  {
    std::fill(data.get(), data.get() + time * depth * height * width, 0);
  }

private:
  int time, depth, height, width;
  std::unique_ptr<T[]> data;
};

#endif
