#include "OutputTools.h"
#include <iostream>

template <typename T> 
void OutputVTI::exportNodeDataVTI(const std::string file, Array2D<T> &node, const size_t dim,
                                    const size_t nx, const size_t ny, const size_t nz, 
                                    const double dx, const double dy, const double dz)
{
    FILE *fp;
    fp=fopen(file.c_str(),"w");

    if(fp == NULL)
    {
        if(mpi.myId == 0) std::cout << file << " open error" << std::endl;
        if(mpi.myId == 0) exit(1);
    }

    if(dim == 2)
    {
        fprintf(fp,"<?xml version=\"1.0\"?>\n");
        fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, 1, 0.0, 0.0, 0.0, dx, dy, 0.0);  
        fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, 1);  
        fprintf(fp,"<PointData>\n");
        fprintf(fp,"<DataArray type=\"Float64\" Name=\"node\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(size_t i=0; i<(nx+1)*(ny+1); i++){
            fprintf(fp,"%e %e 0.0\n", node(i, 0), node(i, 1));
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"</PointData>\n");
        fprintf(fp,"</Piece>\n");
        fprintf(fp,"</ImageData>\n");
        fprintf(fp,"</VTKFile>\n");
        fclose(fp);
    }
    else if(dim == 3)
    {
        fprintf(fp,"<?xml version=\"1.0\"?>\n");
        fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0,nz, 0.0, 0.0, 0.0, dx, dy, dz);  
        fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);  
        fprintf(fp,"<PointData>\n");
        fprintf(fp,"<DataArray type=\"Float64\" Name=\"node\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(size_t i=0; i<(nx+1)*(ny+1)*(nz+1); i++){
            fprintf(fp,"%e %e %e\n", node(i, 0), node(i, 1), node(i, 2));
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"</PointData>\n");
        fprintf(fp,"</Piece>\n");
        fprintf(fp,"</ImageData>\n");
        fprintf(fp,"</VTKFile>\n");
        fclose(fp);
    }
}

template <typename T>
void OutputVTI::exportCellDataVTI(const std::string file, Array2D<T> &cell, const size_t dim,
                                    const size_t nx, const size_t ny, const size_t nz, 
                                    const double dx, const double dy, const double dz)
{
    FILE *fp;
    fp=fopen(file.c_str(),"w");

    if(fp == NULL)
    {
        if(mpi.myId == 0) std::cout << file << " open error" << std::endl;
        if(mpi.myId == 0) exit(1);
    }

    if(dim == 2)
    {
        fprintf(fp,"<?xml version=\"1.0\"?>\n");
        fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, 1, 0.0, 0.0, 0.0, dx, dy, 0.0);  
        fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, 1);  
        fprintf(fp,"<CellData>\n");
        fprintf(fp,"<DataArray type=\"Float64\" Name=\"cell\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(size_t i=0; i<nx*ny; i++){
            fprintf(fp,"%e %e 0.0\n", cell(i, 0), cell(i, 1));
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"</CellData>\n");
        fprintf(fp,"</Piece>\n");
        fprintf(fp,"</ImageData>\n");
        fprintf(fp,"</VTKFile>\n");
        fclose(fp);
    }
    else if(dim == 3)
    {
        fprintf(fp,"<?xml version=\"1.0\"?>\n");
        fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0,nz, 0.0, 0.0, 0.0, dx, dy, dz);  
        fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);  
        fprintf(fp,"<CellData>\n");
        fprintf(fp,"<DataArray type=\"Float64\" Name=\"cell\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(size_t i=0; i<nx*ny*nz; i++){
            fprintf(fp,"%e %e %e\n", cell(i, 0), cell(i, 1), cell(i, 2));
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"</PointData>\n");
        fprintf(fp,"</Piece>\n");
        fprintf(fp,"</ImageData>\n");
        fprintf(fp,"</VTKFile>\n");
        fclose(fp);
    }
}