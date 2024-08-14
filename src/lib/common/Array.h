/**
 * @file   Array.h
 * @autor  K.Ueda
 * @date   May, 2024
 */

#ifndef ARRAY_H
#define ARRAY_H

#include <algorithm>
#include <iostream>
#include <memory>

// Array1D class definition
template <class T> class Array1D
{
public:
  Array1D(int width) : width(width), data(std::make_unique<T[]>(width)) {}
  Array1D() : width(0) {}

  inline T &operator()(int x)
  {
    return data[x];
  }

  inline const T &operator()(int x) const
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

  inline void fillZero()
  {
    std::fill(data.get(), data.get() + width, 0);
  }

  // Scalar multiplication
  Array1D<T> operator*(const T& scalar) const
  {
    Array1D<T> result(width);
    for (int i = 0; i < width; ++i)
    {
      result(i) = (*this)(i) * scalar;
    }
    return result;
  }

  // Scalar division
  Array1D<T> operator/(const T& scalar) const
  {
    Array1D<T> result(width);
    for (int i = 0; i < width; ++i)
    {
      result(i) = (*this)(i) / scalar;
    }
    return result;
  }

  // Element-wise addition
  Array1D<T> operator+(const Array1D<T>& other) const
  {
    Array1D<T> result(width);
    for (int i = 0; i < width; ++i)
    {
      result(i) = (*this)(i) + other(i);
    }
    return result;
  }

  // Element-wise subtraction
  Array1D<T> operator-(const Array1D<T>& other) const
  {
    Array1D<T> result(width);
    for (int i = 0; i < width; ++i)
    {
      result(i) = (*this)(i) - other(i);
    }
    return result;
  }

private:
  int width;
  std::unique_ptr<T[]> data;
};

// Scalar multiplication (commutative) for Array1D
template <class T>
Array1D<T> operator*(const T& scalar, const Array1D<T>& array)
{
  return array * scalar;
}

// Array2D class definition
template <class T> class Array2D
{
public:
  Array2D(int height, int width) : width(width), height(height), data(std::make_unique<T[]>(width * height)) {}
  Array2D() : width(0), height(0) {}

  inline T &operator()(int y, int x)
  {
    return data[y * width + x];
  }

  inline const T &operator()(int y, int x) const
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

  inline void fillZero()
  {
    std::fill(data.get(), data.get() + height * width, 0);
  }

  // Scalar multiplication
  Array2D<T> operator*(const T& scalar) const
  {
    Array2D<T> result(height, width);
    for (int i = 0; i < height; ++i)
    {
      for (int j = 0; j < width; ++j)
      {
        result(i, j) = (*this)(i, j) * scalar;
      }
    }
    return result;
  }

  // Scalar division
  Array2D<T> operator/(const T& scalar) const
  {
    Array2D<T> result(height, width);
    for (int i = 0; i < height; ++i)
    {
      for (int j = 0; j < width; ++j)
      {
        result(i, j) = (*this)(i, j) / scalar;
      }
    }
    return result;
  }

  // Element-wise addition
  Array2D<T> operator+(const Array2D<T>& other) const
  {
    Array2D<T> result(height, width);
    for (int i = 0; i < height; ++i)
    {
      for (int j = 0; j < width; ++j)
      {
        result(i, j) = (*this)(i, j) + other(i, j);
      }
    }
    return result;
  }

  // Element-wise subtraction
  Array2D<T> operator-(const Array2D<T>& other) const
  {
    Array2D<T> result(height, width);
    for (int i = 0; i < height; ++i)
    {
      for (int j = 0; j < width; ++j)
      {
        result(i, j) = (*this)(i, j) - other(i, j);
      }
    }
    return result;
  }

private:
  int width, height;
  std::unique_ptr<T[]> data;
};

// Scalar multiplication (commutative) for Array2D
template <class T>
Array2D<T> operator*(const T& scalar, const Array2D<T>& array)
{
  return array * scalar;
}

// Array3D class definition
template <class T> class Array3D
{
public:
  Array3D(int depth, int height, int width)
      : width(width), height(height), depth(depth), data(std::make_unique<T[]>(depth * height * width)) {}
  Array3D() : width(0), height(0), depth(0) {}

  inline T &operator()(int z, int y, int x)
  {
    return data[z * height * width + y * width + x];
  }

  inline const T &operator()(int z, int y, int x) const
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

  inline void fillZero()
  {
    std::fill(data.get(), data.get() + depth * height * width, 0);
  }

  // Scalar multiplication
  Array3D<T> operator*(const T& scalar) const
  {
    Array3D<T> result(depth, height, width);
    for (int i = 0; i < depth; ++i)
    {
      for (int j = 0; j < height; ++j)
      {
        for (int k = 0; k < width; ++k)
        {
          result(i, j, k) = (*this)(i, j, k) * scalar;
        }
      }
    }
    return result;
  }

  // Scalar division
  Array3D<T> operator/(const T& scalar) const
  {
    Array3D<T> result(depth, height, width);
    for (int i = 0; i < depth; ++i)
    {
      for (int j = 0; j < height; ++j)
      {
        for (int k = 0; k < width; ++k)
        {
          result(i, j, k) = (*this)(i, j, k) / scalar;
        }
      }
    }
    return result;
  }

  // Element-wise addition
  Array3D<T> operator+(const Array3D<T>& other) const
  {
    Array3D<T> result(depth, height, width);
    for (int i = 0; i < depth; ++i)
    {
      for (int j = 0; j < height; ++j)
      {
        for (int k = 0; k < width; ++k)
        {
          result(i, j, k) = (*this)(i, j, k) + other(i, j, k);
        }
      }
    }
    return result;
  }

  // Element-wise subtraction
  Array3D<T> operator-(const Array3D<T>& other) const
  {
    Array3D<T> result(depth, height, width);
    for (int i = 0; i < depth; ++i)
    {
      for (int j = 0; j < height; ++j)
      {
        for (int k = 0; k < width; ++k)
        {
          result(i, j, k) = (*this)(i, j, k) - other(i, j, k);
        }
      }
    }
    return result;
  }

private:
  int width, height, depth;
  std::unique_ptr<T[]> data;
};

// Scalar multiplication (commutative) for Array3D
template <class T>
Array3D<T> operator*(const T& scalar, const Array3D<T>& array)
{
  return array * scalar;
}

// Array4D class definition
template <class T> class Array4D
{
public:
  Array4D(int time, int depth, int height, int width)
      : time(time), depth(depth), height(height), width(width),
        data(std::make_unique<T[]>(time * depth * height * width)) {}
  Array4D() : time(0), depth(0), height(0), width(0) {}

  inline T &operator()(int t, int z, int y, int x)
  {
    return data[t * depth * height * width + z * height * width + y * width + x];
  }

  inline const T &operator()(int t, int z, int y, int x) const
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

  inline void fillZero()
  {
    std::fill(data.get(), data.get() + time * depth * height * width, 0);
  }

  // Scalar multiplication
  Array4D<T> operator*(const T& scalar) const
  {
    Array4D<T> result(time, depth, height, width);
    for (int i = 0; i < time; ++i)
    {
      for (int j = 0; j < depth; ++j)
      {
        for (int k = 0; k < height; ++k)
        {
          for (int l = 0; l < width; ++l)
          {
            result(i, j, k, l) = (*this)(i, j, k, l) * scalar;
          }
        }
      }
    }
    return result;
  }

  // Scalar division
  Array4D<T> operator/(const T& scalar) const
  {
    Array4D<T> result(time, depth, height, width);
    for (int i = 0; i < time; ++i)
    {
      for (int j = 0; j < depth; ++j)
      {
        for (int k = 0; k < height; ++k)
        {
          for (int l = 0; l < width; ++l)
          {
            result(i, j, k, l) = (*this)(i, j, k, l) / scalar;
          }
        }
      }
    }
    return result;
  }

  // Element-wise addition
  Array4D<T> operator+(const Array4D<T>& other) const
  {
    Array4D<T> result(time, depth, height, width);
    for (int i = 0; i < time; ++i)
    {
      for (int j = 0; j < depth; ++j)
      {
        for (int k = 0; k < height; ++k)
        {
          for (int l = 0; l < width; ++l)
          {
            result(i, j, k, l) = (*this)(i, j, k, l) + other(i, j, k, l);
          }
        }
      }
    }
    return result;
  }

  // Element-wise subtraction
  Array4D<T> operator-(const Array4D<T>& other) const
  {
    Array4D<T> result(time, depth, height, width);
    for (int i = 0; i < time; ++i)
    {
      for (int j = 0; j < depth; ++j)
      {
        for (int k = 0; k < height; ++k)
        {
          for (int l = 0; l < width; ++l)
          {
            result(i, j, k, l) = (*this)(i, j, k, l) - other(i, j, k, l);
          }
        }
      }
    }
    return result;
  }

private:
  int time, depth, height, width;
  std::unique_ptr<T[]> data;
};

// Scalar multiplication (commutative) for Array4D
template <class T>
Array4D<T> operator*(const T& scalar, const Array4D<T>& array)
{
  return array * scalar;
}

#endif