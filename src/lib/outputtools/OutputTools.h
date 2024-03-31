#ifndef OUTPUTTOOLS_H
#define OUTPUTTOOLS_H

#include <iostream>
#include "Array.h"
#include "MyMPI.h"
extern MyMPI mpi;

class OutputVTI
{
    public:
        template <typename T> 
        void exportNodeDataVTI(const std::string file, Array2D<T> &node, const size_t dim,
                               const size_t nx, const size_t ny, const size_t nz, 
                               const double dx, const double dy, const double dz);
        template <typename T>
        void exportCellDataVTI(const std::string file, Array2D<T> &cell, const size_t dim,
                               const size_t nx, const size_t ny, const size_t nz, 
                               const double dx, const double dy, const double dz);
};

class OutputDat
{
    public:
};

#endif