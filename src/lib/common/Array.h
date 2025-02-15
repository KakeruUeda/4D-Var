/**
 * @file   Array.h
 * @autor  K.Ueda
 * @date   May, 2024
 */

#ifndef ARRAY_H
#define ARRAY_H

#include <algorithm>
#include <fstream>
#include <iostream>

template <class T> class Array1D
{
public:
  Array1D(int width) : width_(width), data_(new T[width]())
  {
  }
  Array1D() : width_(0), data_(nullptr)
  {
  }
  ~Array1D()
  {
    delete[] data_;
  }

  Array1D(const Array1D &other) : width_(other.width_), data_(new T[other.width_])
  {
    std::copy(other.data_, other.data_ + other.width_, data_);
  }

  Array1D &operator=(const Array1D &other)
  {
    if(this != &other) {
      delete[] data_;
      width_ = other.width_;
      data_ = new T[other.width_];
      std::copy(other.data_, other.data_ + other.width_, data_);
    }
    return *this;
  }

  Array1D(Array1D &&other) noexcept : width_(other.width_), data_(other.data_)
  {
    other.width_ = 0;
    other.data_ = nullptr;
  }

  Array1D &operator=(Array1D &&other) noexcept
  {
    if(this != &other) {
      delete[] data_;
      width_ = other.width_;
      data_ = other.data_;
      other.width_ = 0;
      other.data_ = nullptr;
    }
    return *this;
  }

  inline T &operator()(int x)
  {
    return data_[x];
  }

  inline const T &operator()(int x) const
  {
    return data_[x];
  }

  inline int size() const
  {
    return width_;
  }

  inline void allocate(int newWidth)
  {
    if(data_) {
      delete[] data_;
    }
    data_ = new T[newWidth]();
    width_ = newWidth;
  }

  inline void fillZero()
  {
    std::fill(data_, data_ + width_, 0);
  }

  T *data() const
  {
    return data_;
  }

  /// @brief Export to dat file
  /// @param filename
  void exportDAT(const std::string &filename) const
  {
    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return;
    }
    ofs.write(reinterpret_cast<const char *>(data_), sizeof(T) * width_);
  }

  /// @brief import from dat file
  /// @param filename
  void importDAT(const std::string &filename)
  {
    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs) {
      std::cerr << "Could not open file for reading: " << filename << std::endl;
      return;
    }
    ifs.read(reinterpret_cast<char *>(data_), sizeof(T) * width_);
  }

  /// @brief export to binary file
  /// @param filename
  void exportBIN(const std::string &filename) const
  {
    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return;
    }
    ofs.write(reinterpret_cast<const char *>(&width_), sizeof(width_));
    ofs.write(reinterpret_cast<const char *>(data_), sizeof(T) * width_);
  }

  /// @brief import from binary file
  /// @param filename
  void importBIN(const std::string &filename)
  {
    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs) {
      std::cerr << "Could not open file for reading: " << filename << std::endl;
      return;
    }
    ifs.read(reinterpret_cast<char *>(&width_), sizeof(width_));
    delete[] data_;
    data_ = new T[width_];
    ifs.read(reinterpret_cast<char *>(data_), sizeof(T) * width_);
  }

private:
  int width_;
  T *data_;
};

template <class T> class Array2D
{
public:
  Array2D(int height, int width) : width_(width), height_(height), data_(new T[width * height]())
  {
  }
  Array2D() : width_(0), height_(0), data_(nullptr)
  {
  }
  ~Array2D()
  {
    delete[] data_;
  }

  Array2D(const Array2D &other)
      : width_(other.width_), height_(other.height_), data_(new T[other.width_ * other.height_])
  {
    std::copy(other.data_, other.data_ + (other.width_ * other.height_), data_);
  }

  Array2D &operator=(const Array2D &other)
  {
    if(this != &other) {
      delete[] data_;
      width_ = other.width_;
      height_ = other.height_;
      data_ = new T[other.width_ * other.height_];
      std::copy(other.data_, other.data_ + (other.width_ * other.height_), data_);
    }
    return *this;
  }

  Array2D(Array2D &&other) noexcept : width_(other.width_), height_(other.height_), data_(other.data_)
  {
    other.width_ = 0;
    other.height_ = 0;
    other.data_ = nullptr;
  }

  Array2D &operator=(Array2D &&other) noexcept
  {
    if(this != &other) {
      delete[] data_;
      width_ = other.width_;
      height_ = other.height_;
      data_ = other.data_;
      other.width_ = 0;
      other.height_ = 0;
      other.data_ = nullptr;
    }
    return *this;
  }

  Array2D operator*(T scalar) const
  {
    Array2D result(height_, width_);
    for(int y = 0; y < height_; ++y) {
      for(int x = 0; x < width_; ++x) {
        result(y, x) = (*this)(y, x) * scalar;
      }
    }
    return result;
  }

  Array2D operator-() const
  {
    Array2D result(height_, width_);
    for(int y = 0; y < height_; ++y) {
      for(int x = 0; x < width_; ++x) {
        result(y, x) = -(*this)(y, x);
      }
    }
    return result;
  }

  Array2D operator+(const Array2D &other) const
  {
    if(width_ != other.width_ || height_ != other.height_) {
      throw std::invalid_argument("Array2D dimensions must match for addition");
    }
    Array2D result(height_, width_);
    for(int y = 0; y < height_; ++y) {
      for(int x = 0; x < width_; ++x) {
        result(y, x) = (*this)(y, x) + other(y, x);
      }
    }
    return result;
  }

  inline T &operator()(int x)
  {
    return data_[x];
  }

  inline const T &operator()(int x) const
  {
    return data_[x];
  }

  inline T &operator()(int y, int x)
  {
    return data_[y * width_ + x];
  }

  inline const T &operator()(int y, int x) const
  {
    return data_[y * width_ + x];
  }

  inline int size() const
  {
    return width_ * height_;
  }

  inline void allocate(int newHeight, int newWidth)
  {
    if(data_) {
      delete[] data_;
    }
    data_ = new T[newHeight * newWidth]();
    height_ = newHeight;
    width_ = newWidth;
  }

  inline void fillZero()
  {
    std::fill(data_, data_ + height_ * width_, 0);
  }

  int height() const
  {
    return height_;
  }

  int width() const
  {
    return width_;
  }

  T *data() const
  {
    return data_;
  }

  /// @brief export to dat file
  /// @param filename
  void exportDAT(const std::string &filename) const
  {
    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return;
    }
    for(int y = 0; y < height_; ++y) {
      ofs.write(reinterpret_cast<const char *>(&data_[y * width_]), sizeof(T) * width_);
    }
  }

  /// @brief import from dat file
  /// @param filename File name.
  void importDAT(const std::string &filename)
  {
    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs) {
      std::cerr << "Could not open file for reading: " << filename << std::endl;
      return;
    }
    for(int y = 0; y < height_; ++y) {
      ifs.read(reinterpret_cast<char *>(&data_[y * width_]), sizeof(T) * width_);
    }
  }

  /// @brief Export to binary file
  /// @param filename File name.
  void exportBIN(const std::string &filename) const
  {
    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return;
    }
    ofs.write(reinterpret_cast<const char *>(&width_), sizeof(width_));
    ofs.write(reinterpret_cast<const char *>(&height_), sizeof(height_));

    ofs.write(reinterpret_cast<const char *>(data_), sizeof(T) * width_ * height_);
  }

  /// @brief Import from binary file
  /// @param filename file name.
  void importBIN(const std::string &filename)
  {
    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs) {
      std::cerr << "Could not open file for reading: " << filename << std::endl;
      return;
    }
    ifs.read(reinterpret_cast<char *>(&width_), sizeof(width_));
    ifs.read(reinterpret_cast<char *>(&height_), sizeof(height_));
    delete[] data_;
    data_ = new T[width_ * height_];
    ifs.read(reinterpret_cast<char *>(data_), sizeof(T) * width_ * height_);
  }

  /// @brief Export specific row data to dat file
  /// @param filename Export file name.
  /// @param y Row index.
  void exportDAT(const std::string &filename, int y) const
  {
    if(y < 0 || y >= height_) {
      std::cerr << "Invalid height specified: " << y << std::endl;
      return;
    }

    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return;
    }

    ofs.write(reinterpret_cast<const char *>(&data_[y * width_]), sizeof(T) * width_);
  }

  /// @brief Import specific row data from dat file.
  /// @param filename File name.
  /// @param y Row index.
  void importDAT(const std::string &filename, int y)
  {
    if(y < 0 || y >= height_) {
      std::cerr << "Invalid height specified: " << y << std::endl;
      return;
    }

    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs) {
      std::cerr << "Could not open file for reading: " << filename << std::endl;
      return;
    }

    ifs.read(reinterpret_cast<char *>(&data_[y * width_]), sizeof(T) * width_);
  }

  /// @brief Export specific row data to binary file.
  /// @param filename File name.
  /// @param y Row index.
  void exportBIN(const std::string &filename, int y) const
  {
    if(y < 0 || y >= height_) {
      std::cerr << "Invalid height specified: " << y << std::endl;
      return;
    }

    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return;
    }

    ofs.write(reinterpret_cast<const char *>(&width_), sizeof(width_));
    ofs.write(reinterpret_cast<const char *>(&data_[y * width_]), sizeof(T) * width_);
  }

  /// @brief Import specific row data from binary file.
  /// @param filename File name.
  /// @param y Row index.
  void importBIN(const std::string &filename, int y)
  {
    if(y < 0 || y >= height_) {
      std::cerr << "Invalid height specified: " << y << std::endl;
      return;
    }

    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs) {
      std::cerr << "Could not open file for reading: " << filename << std::endl;
      return;
    }

    long long int fileWidth, fileHeight;
    ifs.read(reinterpret_cast<char *>(&fileWidth), sizeof(fileWidth));

    if(fileWidth != width_) {
      std::cerr << "File dimensions do not match the current Array2D dimensions" << std::endl;
      return;
    }

    ifs.read(reinterpret_cast<char *>(&data_[y * width_]), sizeof(T) * width_);
  }

private:
 long int width_, height_;
  T *data_;
};

template <typename T> Array2D<T> operator*(T scalar, const Array2D<T> &array)
{
  Array2D<T> result(array.height(), array.width());
  for(int y = 0; y < array.height(); ++y) {
    for(int x = 0; x < array.width(); ++x) {
      result(y, x) = scalar * array(y, x);
    }
  }
  return result;
}

template <typename T> Array2D<T> operator*(const Array2D<T> &array, T scalar)
{
  return scalar * array;
}

template <class T> class Array3D
{
public:
  Array3D(long long int depth, long long int height, long long int width)
      : width_(width), height_(height), depth_(depth), data_(new T[depth * height * width]())
  {
  }
  Array3D() : width_(0), height_(0), depth_(0), data_(nullptr)
  {
  }
  ~Array3D()
  {
    delete[] data_;
  }

  Array3D(const Array3D &other)
      : width_(other.width_), height_(other.height_), depth_(other.depth_),
        data_(new T[other.depth_ * other.height_ * other.width_])
  {
    std::copy(other.data_, other.data_ + (other.depth_ * other.height_ * other.width_), data_);
  }

  Array3D &operator=(const Array3D &other)
  {
    if(this != &other) {
      delete[] data_;
      width_ = other.width_;
      height_ = other.height_;
      depth_ = other.depth_;
      data_ = new T[other.depth_ * other.height_ * other.width_];
      std::copy(other.data_, other.data_ + (other.depth_ * other.height_ * other.width_), data_);
    }
    return *this;
  }

  Array3D(Array3D &&other) noexcept
      : width_(other.width_), height_(other.height_), depth_(other.depth_), data_(other.data_)
  {
    other.width_ = 0;
    other.height_ = 0;
    other.depth_ = 0;
    other.data_ = nullptr;
  }

  Array3D &operator=(Array3D &&other) noexcept
  {
    if(this != &other) {
      delete[] data_;
      width_ = other.width_;
      height_ = other.height_;
      depth_ = other.depth_;
      data_ = other.data_;
      other.width_ = 0;
      other.height_ = 0;
      other.depth_ = 0;
      other.data_ = nullptr;
    }
    return *this;
  }

  Array3D operator*(T scalar) const
  {
    Array3D result(depth_, height_, width_);
    for(long long int z = 0; z < depth_; ++z) {
      for(long long int y = 0; y < height_; ++y) {
        for(long long int x = 0; x < width_; ++x) {
          result(z, y, x) = (*this)(z, y, x) * scalar;
        }
      }
    }
    return result;
  }

  Array3D operator-() const
  {
    Array3D result(depth_, height_, width_);
    for(long long int z = 0; z < depth_; ++z) {
      for(long long int y = 0; y < height_; ++y) {
        for(long long int x = 0; x < width_; ++x) {
          result(z, y, x) = -(*this)(z, y, x);
        }
      }
    }
    return result;
  }

  Array3D operator+(const Array3D &other) const
  {
    if(width_ != other.width_ || height_ != other.height_ || depth_ != other.depth_) {
      throw std::invalid_argument("Array3D dimensions must match for addition");
    }
    Array3D result(depth_, height_, width_);
    for(long long int z = 0; z < depth_; ++z) {
      for(long long int y = 0; y < height_; ++y) {
        for(long long int x = 0; x < width_; ++x) {
          result(z, y, x) = (*this)(z, y, x) + other(z, y, x);
        }
      }
    }
    return result;
  }

  inline T &operator()(long long int x)
  {
    return data_[x];
  }

  inline const T &operator()(long long int x) const
  {
    return data_[x];
  }

  inline T &operator()(long long int z, long long int y, long long int x)
  {
    return data_[z * height_ * width_ + y * width_ + x];
  }

  inline const T &operator()(long long int z, long long int y, long long int x) const
  {
    return data_[z * height_ * width_ + y * width_ + x];
  }

  inline long long int size() const
  {
    return depth_ * height_ * width_;
  }

  inline void allocate(long long int newDepth, long long int newHeight, long long int newWidth)
  {
    if(data_) {
      delete[] data_;
    }
    data_ = new T[newDepth * newHeight * newWidth]();
    depth_ = newDepth;
    height_ = newHeight;
    width_ = newWidth;
  }

  inline void fillZero()
  {
    std::fill(data_, data_ + depth_ * height_ * width_, 0);
  }

  long long int width() const
  {
    return width_;
  }

  long long int height() const
  {
    return height_;
  }

  long long int depth() const
  {
    return depth_;
  }

  T *data() const
  {
    return data_;
  }

  /// @brief Export from dat file
  /// @param filename File name.
  void exportDAT(const std::string &filename) const
  {
    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return;
    }
    for(long long int z = 0; z < depth_; ++z) {
      for(long long int y = 0; y < height_; ++y) {
        ofs.write(reinterpret_cast<const char *>(&data_[z * height_ * width_ + y * width_]), sizeof(T) * width_);
      }
    }
  }

  /// @brief Import from dat file
  /// @param filename File name.
  void importDAT(const std::string &filename)
  {
    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs) {
      std::cerr << "Could not open file for reading: " << filename << std::endl;
      return;
    }
    for(long long int z = 0; z < depth_; ++z) {
      for(long long int y = 0; y < height_; ++y) {
        ifs.read(reinterpret_cast<char *>(&data_[z * height_ * width_ + y * width_]), sizeof(T) * width_);
      }
    }
  }

  /// @brief Export to binary file
  /// @param filename File name.
  void exportBIN(const std::string &filename) const
  {
    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return;
    }
    ofs.write(reinterpret_cast<const char *>(&width_), sizeof(width_));
    ofs.write(reinterpret_cast<const char *>(&height_), sizeof(height_));
    ofs.write(reinterpret_cast<const char *>(&depth_), sizeof(depth_));
    ofs.write(reinterpret_cast<const char *>(data_), sizeof(T) * size());
  }

  /// @brief Import data from binary file
  /// @param filename File name.
  void importBIN(const std::string &filename)
  {
    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs) {
      std::cerr << "Could not open file for reading: " << filename << std::endl;
      return;
    }
    ifs.read(reinterpret_cast<char *>(&width_), sizeof(width_));
    ifs.read(reinterpret_cast<char *>(&height_), sizeof(height_));
    ifs.read(reinterpret_cast<char *>(&depth_), sizeof(depth_));
    delete[] data_;
    data_ = new T[size()];
    ifs.read(reinterpret_cast<char *>(data_), sizeof(T) * size());
  }

  /// @brief Export specific slice data to dat file binary file.
  /// @param filename Export file name.
  /// @param z Slice index.
  void exportDAT(const std::string &filename, long long int z) const
  {
    if(z < 0 || z >= depth_) {
      std::cerr << "Invalid depth specified: " << z << std::endl;
      return;
    }

    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return;
    }

    for(long long int y = 0; y < height_; ++y) {
      ofs.write(reinterpret_cast<const char *>(&data_[z * height_ * width_ + y * width_]), sizeof(T) * width_);
    }
  }

  /// @brief Import specific slice data to dat file.
  /// @param filename Export file name.
  /// @param z Slice index.
  void importDAT(const std::string &filename, long long int z)
  {
    if(z < 0 || z >= depth_) {
      std::cerr << "Invalid depth specified: " << z << std::endl;
      return;
    }

    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs) {
      std::cerr << "Could not open file for reading: " << filename << std::endl;
      return;
    }

    for(long long int y = 0; y < height_; ++y) {
      ifs.read(reinterpret_cast<char *>(&data_[z * height_ * width_ + y * width_]), sizeof(T) * width_);
    }
  }

  /// @brief Import specific slice data to binary file.
  /// @param filename File name.
  /// @param z Slice index.
  void exportBIN(const std::string &filename, long long int z) const
  {
    if(z < 0 || z >= depth_) {
      std::cerr << "Invalid depth specified: " << z << std::endl;
      return;
    }

    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs) {
      std::cerr << "Could not open file for writing: " << filename << std::endl;
      return;
    }

    ofs.write(reinterpret_cast<const char *>(&width_), sizeof(width_));
    ofs.write(reinterpret_cast<const char *>(&height_), sizeof(height_));
    ofs.write(reinterpret_cast<const char *>(&data_[z * height_ * width_]), sizeof(T) * height_ * width_);
  }

  /// @brief Import specific slice data from binary file.
  /// @param filename File name.
  /// @param z Slice index.
  void importBIN(const std::string &filename, long long int z)
  {
    if(z < 0 || z >= depth_) {
      std::cerr << "Invalid depth specified: " << z << std::endl;
      return;
    }
    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs) {
      std::cerr << "Could not open file for reading: " << filename << std::endl;
      return;
    }

    long int fileWidth, fileHeight, fileDepth;
    ifs.read(reinterpret_cast<char *>(&fileWidth), sizeof(fileWidth));
    ifs.read(reinterpret_cast<char *>(&fileHeight), sizeof(fileHeight));

    if(fileWidth != width_ || fileHeight != height_) {
      std::cerr << "File dimensions do not match the current Array3D dimensions" << std::endl;
      //std::cout << z * height_ * width_ << " " << fileWidth << " " << width_ << " " << fileHeight << " " << height_ << std::endl;
      return;
    }

    ifs.read(reinterpret_cast<char *>(&data_[z * height_ * width_]), sizeof(T) * height_ * width_);
  }

private:
  long long int width_, height_, depth_;
  T *data_;
};

template <typename T, typename Scalar> Array3D<T> operator*(Scalar scalar, const Array3D<T> &array)
{
  Array3D<T> result(array.depth(), array.height(), array.width());
  for(long long int z = 0; z < array.depth(); ++z) {
    for(long long int y = 0; y < array.height(); ++y) {
      for(long long int x = 0; x < array.width(); ++x) {
        result(z, y, x) = scalar * array(z, y, x);
      }
    }
  }
  return result;
}

template <typename T, typename Scalar> Array3D<T> operator*(const Array3D<T> &array, Scalar scalar)
{
  return scalar * array;
}

#endif