/**
 * @file VTI.h
 * @author K.Ueda
 * @date Jun, 2024
*/

#include "VTK.h"

void VTK::exportDataVTI(const std::string file, DataGrid &data,
                              const int &t, const int &dim)
{
    FILE *fp;
    fp=fopen(file.c_str(),"w");
    if(fp == NULL){
        if(mpi.myId == 0) std::cout << file << " open error" << std::endl;
        if(mpi.myId == 0) exit(1);
    }

    if(dim == 2){
        fprintf(fp,"<?xml version=\"1.0\"?>\n");
        fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, data.nx, 0, data.ny, 0, 1, 0.0, 0.0, 0.0, data.dx, data.dy, 0.0);  
        fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, data.nx, 0, data.ny, 0, 1);  
        fprintf(fp,"<CellData>\n");
        fprintf(fp,"<DataArray type=\"Float64\" Name=\"cell\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(int j=0; j<data.ny; j++){
            for(int i=0; i<data.nx; i++){
                fprintf(fp,"%e %e 0e0\n", data(j, i).vCFD[t][0], data(j, i).vCFD[t][1]);
            }
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"</CellData>\n");
        fprintf(fp,"</Piece>\n");
        fprintf(fp,"</ImageData>\n");
        fprintf(fp,"</VTKFile>\n");
        fclose(fp);
    }
    else if(dim == 3){
        fprintf(fp,"<?xml version=\"1.0\"?>\n");
        fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, data.nx, 0, data.ny, 0,data.nz, 0.0, 0.0, 0.0, data.dx, data.dy, data.dz);  
        fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, data.nx, 0, data.ny, 0, data.nz);  
        fprintf(fp,"<CellData>\n");
        fprintf(fp,"<DataArray type=\"Float64\" Name=\"data\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(int k=0; k<data.nz; k++){
            for(int j=0; j<data.ny; j++){
                for(int i=0; i<data.nx; i++){
                    fprintf(fp,"%e %e %e\n", data(k, j, i).vCFD[t][0], data(k, j, i).vCFD[t][1], data(k, j, i).vCFD[t][2]);
                }
            }
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"</CellData>\n");
        fprintf(fp,"</Piece>\n");
        fprintf(fp,"</ImageData>\n");
        fprintf(fp,"</VTKFile>\n");
        fclose(fp);
    }
}


void VTK::exportVelocityDataVTI(const std::string file, DataGrid &data,
                                   const int &t, const int &dim)
{
    FILE *fp;
    fp=fopen(file.c_str(),"w");
    if(fp == NULL){
        if(mpi.myId == 0) std::cout << file << " open error" << std::endl;
        if(mpi.myId == 0) exit(1);
    }

    if(dim == 2){
        fprintf(fp,"<?xml version=\"1.0\"?>\n");
        fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, data.nx, 0, data.ny, 0, 1, 0.0, 0.0, 0.0, data.dx, data.dy, 0.0);  
        fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, data.nx, 0, data.ny, 0, 1);  
        fprintf(fp,"<CellData>\n");
        fprintf(fp,"<DataArray type=\"Float64\" Name=\"cell\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(int j=0; j<data.ny; j++){
            for(int i=0; i<data.nx; i++){
                fprintf(fp,"%e %e 0e0\n", data(j, i).vCFD[t][0], data(j, i).vCFD[t][1]);
            }
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"</CellData>\n");
        fprintf(fp,"</Piece>\n");
        fprintf(fp,"</ImageData>\n");
        fprintf(fp,"</VTKFile>\n");
        fclose(fp);
    }
    else if(dim == 3){
        fprintf(fp,"<?xml version=\"1.0\"?>\n");
        fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, data.nx, 0, data.ny, 0,data.nz, 0.0, 0.0, 0.0, data.dx, data.dy, data.dz);  
        fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, data.nx, 0, data.ny, 0, data.nz);  
        fprintf(fp,"<CellData>\n");
        fprintf(fp,"<DataArray type=\"Float64\" Name=\"vCFD\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(int k=0; k<data.nz; k++){
            for(int j=0; j<data.ny; j++){
                for(int i=0; i<data.nx; i++){
                    fprintf(fp,"%e %e %e\n", data(k, j, i).vCFD[t][0], data(k, j, i).vCFD[t][1], data(k, j, i).vCFD[t][2]);
                }
            }
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"<DataArray type=\"Float64\" Name=\"vMRI\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(int k=0; k<data.nz; k++){
            for(int j=0; j<data.ny; j++){
                for(int i=0; i<data.nx; i++){
                    fprintf(fp,"%e %e %e\n", data(k, j, i).vMRI[t][0], data(k, j, i).vMRI[t][1], data(k, j, i).vMRI[t][2]);
                }
            }
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"<DataArray type=\"Float64\" Name=\"ve\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(int k=0; k<data.nz; k++){
            for(int j=0; j<data.ny; j++){
                for(int i=0; i<data.nx; i++){
                    fprintf(fp,"%e %e %e\n", data(k, j, i).ve[t][0], data(k, j, i).ve[t][1], data(k, j, i).ve[t][2]);
                }
            }
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"</CellData>\n");
        fprintf(fp,"</Piece>\n");
        fprintf(fp,"</ImageData>\n");
        fprintf(fp,"</VTKFile>\n");
        fclose(fp);
    }
}

