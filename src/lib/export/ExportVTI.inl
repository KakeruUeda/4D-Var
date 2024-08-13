/**
 * @file ExportVTI.inl
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef EXPORT_VTI_INL_H
#define EXPORT_VTI_INL_H

#include <fstream>
#include <iostream>
#include <cstdio>
#include "Export.h"

template<typename T>
void EXPORT::exportScalarPointDataVTI(const std::string &file, const char *dataName, std::vector<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
	FILE *fp;
	fp=fopen(file.c_str(), "w");

	if(fp==NULL)
	{
		std::cerr<<file<<" open error"<<std::endl;
		exit(1);
	}
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp, "<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);
	fprintf(fp, "<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
	fprintf(fp, "<PointData>\n");
	fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"0\">\n", dataName);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</PointData>\n");
	fprintf(fp, "<CellData>\n");
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</Piece>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fclose(fp);

	unsigned long allsize=(nx+1)*(ny+1)*(nz+1)*sizeof(float);

	float *data=new float[(nx+1)*(ny+1)*(nz+1)];
	for(int k=0; k<nz+1; k++)
	{
		for(int j=0; j<ny+1; j++)
		{
			for(int i=0; i<nx+1; i++)
			{
				int n=i+j*(nx+1)+k*(nx+1)*(ny+1);
				data[n]=static_cast<float>(p[n]);
			}
		}
	}

	std::fstream ofs;
	ofs.open(file.c_str(), std::ios::out|std::ios::app|std::ios_base::binary);
	ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
	ofs.write(reinterpret_cast<char *>(data), allsize);
	ofs.close();

	delete[] data;

	fp=fopen(file.c_str(), "a");
	fprintf(fp, "\n");
	fprintf(fp, "</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

template<typename T>
void EXPORT::exportVectorPointDataVTI(const std::string &file, const char *dataName, std::vector<std::vector<T> > &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
	FILE *fp;
	fp=fopen(file.c_str(), "w");

	if(fp==NULL)
	{
		std::cerr<<file<<" open error"<<std::endl;
		exit(1);
	}
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp, "<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);
	fprintf(fp, "<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
	fprintf(fp, "<PointData>\n");
	fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\">\n", dataName);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</PointData>\n");
	fprintf(fp, "<CellData>\n");
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</Piece>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fclose(fp);

	unsigned long allsize=(nx+1)*(ny+1)*(nz+1)*3*sizeof(float);

	float *data=new float[(nx+1)*(ny+1)*(nz+1)*3];
	for(int k=0; k<nz+1; k++)
	{
		for(int j=0; j<ny+1; j++)
		{
			for(int i=0; i<nx+1; i++)
			{
				int n=i+j*(nx+1)+k*(nx+1)*(ny+1);
				data[0+i*3+j*(nx+1)*3+k*(nx+1)*(ny+1)*3]=static_cast<float>(p[n][0]);
				data[1+i*3+j*(nx+1)*3+k*(nx+1)*(ny+1)*3]=static_cast<float>(p[n][1]);
				data[2+i*3+j*(nx+1)*3+k*(nx+1)*(ny+1)*3]=static_cast<float>(p[n][2]);
			}
		}
	}

	std::fstream ofs;
	ofs.open(file.c_str(), std::ios::out|std::ios::app|std::ios_base::binary);
	ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
	ofs.write(reinterpret_cast<char *>(data), allsize);
	ofs.close();

	delete[] data;

	fp=fopen(file.c_str(), "a");
	fprintf(fp, "\n");
	fprintf(fp, "</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

template<typename T>
void EXPORT::exportScalarCellDataVTI(const std::string &file, const char *dataName, std::vector<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
	FILE *fp;
	fp=fopen(file.c_str(), "w");

	if(fp==NULL)
	{
		std::cerr<<file<<" open error"<<std::endl;
		exit(1);
	}
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp, "<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);
	fprintf(fp, "<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
	fprintf(fp, "<PointData>\n");
	fprintf(fp, "</PointData>\n");
	fprintf(fp, "<CellData>\n");
	fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"0\">\n", dataName);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</Piece>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fclose(fp);

	float *data=new float[nx*ny*nz];
	unsigned long allsize=nx*ny*nz*sizeof(float);
	for(int k=0; k<nz; k++)
	{
		for(int j=0; j<ny; j++)
		{
			for(int i=0; i<nx; i++)
			{
				int n=i+j*nx+k*nx*ny;
				data[n]=static_cast<float>(c[n]);
			}
		}
	}

	std::fstream ofs;
	ofs.open(file.c_str(), std::ios::out|std::ios::app|std::ios_base::binary);
	ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
	ofs.write(reinterpret_cast<char *>(data), allsize);
	ofs.close();

	delete[] data;

	fp=fopen(file.c_str(), "a");
	fprintf(fp, "\n");
	fprintf(fp, "</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

template<typename T>
void EXPORT::exportVectorCellDataVTI(const std::string &file, const char *dataName, std::vector<std::vector<T> > &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
	FILE *fp;
	fp=fopen(file.c_str(), "w");

	if(fp==NULL)
	{
		std::cerr<<file<<" open error"<<std::endl;
		exit(1);
	}
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp, "<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);
	fprintf(fp, "<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
	fprintf(fp, "<PointData>\n");
	fprintf(fp, "</PointData>\n");
	fprintf(fp, "<CellData>\n");
	fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\">\n", dataName);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</Piece>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fclose(fp);

	float *data=new float[nx*ny*nz*3];
	unsigned long allsize=nx*ny*nz*3*sizeof(float);
	for(int k=0; k<nz; k++)
	{
		for(int j=0; j<ny; j++)
		{
			for(int i=0; i<nx; i++)
			{
				int n=i+j*nx+k*nx*ny;
				data[0+i*3+j*nx*3+k*nx*ny*3]=static_cast<float>(c[n][0]);
				data[1+i*3+j*nx*3+k*nx*ny*3]=static_cast<float>(c[n][1]);
				data[2+i*3+j*nx*3+k*nx*ny*3]=static_cast<float>(c[n][2]);
			}
		}
	}

	std::fstream ofs;
	ofs.open(file.c_str(), std::ios::out|std::ios::app|std::ios_base::binary);
	ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
	ofs.write(reinterpret_cast<char *>(data), allsize);
	ofs.close();

	delete[] data;

	fp=fopen(file.c_str(), "a");
	fprintf(fp, "\n");
	fprintf(fp, "</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

template<typename T>
void EXPORT::exportScalarPointDataVTI(const std::string &file, const char *dataName, Array1D<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
	FILE *fp;
	fp=fopen(file.c_str(), "w");

	if(fp==NULL)
	{
		std::cerr<<file<<" open error"<<std::endl;
		exit(1);
	}
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp, "<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);
	fprintf(fp, "<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
	fprintf(fp, "<PointData>\n");
	fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"0\">\n", dataName);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</PointData>\n");
	fprintf(fp, "<CellData>\n");
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</Piece>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fclose(fp);

	unsigned long allsize=(nx+1)*(ny+1)*(nz+1)*sizeof(float);

	float *data=new float[(nx+1)*(ny+1)*(nz+1)];
	for(int k=0; k<nz+1; k++)
	{
		for(int j=0; j<ny+1; j++)
		{
			for(int i=0; i<nx+1; i++)
			{
				int n=i+j*(nx+1)+k*(nx+1)*(ny+1);
				data[n]=static_cast<float>(p(n));
			}
		}
	}

	std::fstream ofs;
	ofs.open(file.c_str(), std::ios::out|std::ios::app|std::ios_base::binary);
	ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
	ofs.write(reinterpret_cast<char *>(data), allsize);
	ofs.close();

	delete[] data;

	fp=fopen(file.c_str(), "a");
	fprintf(fp, "\n");
	fprintf(fp, "</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

template<typename T>
void EXPORT::exportVectorPointDataVTI(const std::string &file, const char *dataName, Array2D<T> &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
	FILE *fp;
	fp=fopen(file.c_str(), "w");

	if(fp==NULL)
	{
		std::cerr<<file<<" open error"<<std::endl;
		exit(1);
	}
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp, "<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);
	fprintf(fp, "<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
	fprintf(fp, "<PointData>\n");
	fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\">\n", dataName);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</PointData>\n");
	fprintf(fp, "<CellData>\n");
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</Piece>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fclose(fp);

	unsigned long allsize=(nx+1)*(ny+1)*(nz+1)*3*sizeof(float);

	float *data=new float[(nx+1)*(ny+1)*(nz+1)*3];
	for(int k=0; k<nz+1; k++)
	{
		for(int j=0; j<ny+1; j++)
		{
			for(int i=0; i<nx+1; i++)
			{
				int n=i+j*(nx+1)+k*(nx+1)*(ny+1);
				data[0+i*3+j*(nx+1)*3+k*(nx+1)*(ny+1)*3]=static_cast<float>(p(n, 0));
				data[1+i*3+j*(nx+1)*3+k*(nx+1)*(ny+1)*3]=static_cast<float>(p(n, 1));
				data[2+i*3+j*(nx+1)*3+k*(nx+1)*(ny+1)*3]=static_cast<float>(p(n, 2));
			}
		}
	}

	std::fstream ofs;
	ofs.open(file.c_str(), std::ios::out|std::ios::app|std::ios_base::binary);
	ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
	ofs.write(reinterpret_cast<char *>(data), allsize);
	ofs.close();

	delete[] data;

	fp=fopen(file.c_str(), "a");
	fprintf(fp, "\n");
	fprintf(fp, "</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

template<typename T>
void EXPORT::exportScalarCellDataVTI(const std::string &file, const char *dataName, Array1D<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
	FILE *fp;
	fp=fopen(file.c_str(), "w");

	if(fp==NULL)
	{
		std::cerr<<file<<" open error"<<std::endl;
		exit(1);
	}
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp, "<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);
	fprintf(fp, "<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
	fprintf(fp, "<PointData>\n");
	fprintf(fp, "</PointData>\n");
	fprintf(fp, "<CellData>\n");
	fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"0\">\n", dataName);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</Piece>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fclose(fp);

	float *data=new float[nx*ny*nz];
	unsigned long allsize=nx*ny*nz*sizeof(float);
	for(int k=0; k<nz; k++)
	{
		for(int j=0; j<ny; j++)
		{
			for(int i=0; i<nx; i++)
			{
				int n=i+j*nx+k*nx*ny;
				data[n]=static_cast<float>(c(n));
			}
		}
	}

	std::fstream ofs;
	ofs.open(file.c_str(), std::ios::out|std::ios::app|std::ios_base::binary);
	ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
	ofs.write(reinterpret_cast<char *>(data), allsize);
	ofs.close();

	delete[] data;

	fp=fopen(file.c_str(), "a");
	fprintf(fp, "\n");
	fprintf(fp, "</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

template<typename T>
void EXPORT::exportVectorCellDataVTI(const std::string &file, const char *dataName, Array2D<T> &c, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
	FILE *fp;
	fp=fopen(file.c_str(), "w");

	if(fp==NULL)
	{
		std::cerr<<file<<" open error"<<std::endl;
		exit(1);
	}
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fp, "<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);
	fprintf(fp, "<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
	fprintf(fp, "<PointData>\n");
	fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\">\n", dataName);
	fprintf(fp, "</DataArray>\n");
	fprintf(fp, "</PointData>\n");
	fprintf(fp, "<CellData>\n");
	fprintf(fp, "</CellData>\n");
	fprintf(fp, "</Piece>\n");
	fprintf(fp, "</ImageData>\n");
	fprintf(fp, "<AppendedData encoding=\"raw\">\n");
	fprintf(fp, "_");
	fclose(fp);

	float *data=new float[nx*ny*nz*3];
	unsigned long allsize=nx*ny*nz*3*sizeof(float);
	for(int k=0; k<nz; k++)
	{
		for(int j=0; j<ny; j++)
		{
			for(int i=0; i<nx; i++)
			{
				int n=i+j*nx+k*nx*ny;
				data[0+i*3+j*nx*3+k*nx*ny*3]=static_cast<float>(c(n, 0));
				data[1+i*3+j*nx*3+k*nx*ny*3]=static_cast<float>(c(n, 1));
				data[2+i*3+j*nx*3+k*nx*ny*3]=static_cast<float>(c(n, 2));
			}
		}
	}

	std::fstream ofs;
	ofs.open(file.c_str(), std::ios::out|std::ios::app|std::ios_base::binary);
	ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
	ofs.write(reinterpret_cast<char *>(data), allsize);
	ofs.close();

	delete[] data;

	fp=fopen(file.c_str(), "a");
	fprintf(fp, "\n");
	fprintf(fp, "</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");
	fclose(fp);
}

#endif // EXPORT_VTI_INL_H