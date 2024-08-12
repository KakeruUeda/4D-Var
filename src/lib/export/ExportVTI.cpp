/**
 * @file DAT.inl
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef EXPORT_VTI_H
#define EXPORT_VTI_H

#include "Export.h"

void EXPORT::exportVelocityDataVTI(const std::string &file, DataGrid &data, const int t)
{
  FILE *fp = fopen(file.c_str(), "w");

  if (fp == NULL)
  {
    std::cout << file << " open error" << std::endl;
    exit(1);
  }

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"%e %e %e\" Spacing=\"%e %e %e\">\n", 0, data.nx, 0, data.ny, 0, data.nz, 0e0, 0e0, 0e0, data.dx, data.dy, data.dz);
  fprintf(fp, "<Piece Extent=\"%d %d %d %d %d %d\">\n", 0, data.nx, 0, data.ny, 0, data.nz);
  fprintf(fp, "<PointData>\n");
  fprintf(fp, "</PointData>\n");
  fprintf(fp, "<CellData>\n");

  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"vCFD\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(unsigned long) + data.nx * data.ny * data.nz * 3 * sizeof(float);
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"vMRI\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(unsigned long) + data.nx * data.ny * data.nz * 3 * sizeof(float);
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"ve\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</ImageData>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");

  fclose(fp);

  std::fstream ofs;
  ofs.open(file.c_str(), std::ios::out | std::ios::app | std::ios::binary);
  unsigned long allsize = data.nx * data.ny * data.nz * 3 * sizeof(float);

  // Write vCFD data
  float *data1 = new float[data.nx * data.ny * data.nz * 3];
  for (int k = 0; k < data.nz; k++)
  {
    for (int j = 0; j < data.ny; j++)
    {
      for (int i = 0; i < data.nx; i++)
      {
        int index = i + j * data.nx + k * data.nx * data.ny;
        data1[0 + i * 3 + j * data.nx * 3 + k * data.nx * data.ny * 3] = static_cast<float>(data(k, j, i).vCFD[t][0]);
        data1[1 + i * 3 + j * data.nx * 3 + k * data.nx * data.ny * 3] = static_cast<float>(data(k, j, i).vCFD[t][1]);
        data1[2 + i * 3 + j * data.nx * 3 + k * data.nx * data.ny * 3] = static_cast<float>(data(k, j, i).vCFD[t][2]);
      }
    }
  }

  ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
  ofs.write(reinterpret_cast<char *>(data1), allsize);
  delete[] data1;

  // Write vMRI data
  float *data2 = new float[data.nx * data.ny * data.nz * 3];
  for (int k = 0; k < data.nz; k++)
  {
    for (int j = 0; j < data.ny; j++)
    {
      for (int i = 0; i < data.nx; i++)
      {
        int index = i + j * data.nx + k * data.nx * data.ny;
        data2[0 + i * 3 + j * data.nx * 3 + k * data.nx * data.ny * 3] = static_cast<float>(data(k, j, i).vMRI[t][0]);
        data2[1 + i * 3 + j * data.nx * 3 + k * data.nx * data.ny * 3] = static_cast<float>(data(k, j, i).vMRI[t][1]);
        data2[2 + i * 3 + j * data.nx * 3 + k * data.nx * data.ny * 3] = static_cast<float>(data(k, j, i).vMRI[t][2]);
      }
    }
  }

  ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
  ofs.write(reinterpret_cast<char *>(data2), allsize);
  delete[] data2;

  // Write ve data
  float *data3 = new float[data.nx * data.ny * data.nz * 3];
  for (int k = 0; k < data.nz; k++)
  {
    for (int j = 0; j < data.ny; j++)
    {
      for (int i = 0; i < data.nx; i++)
      {
        int index = i + j * data.nx + k * data.nx * data.ny;
        data3[0 + i * 3 + j * data.nx * 3 + k * data.nx * data.ny * 3] = static_cast<float>(data(k, j, i).ve[t][0]);
        data3[1 + i * 3 + j * data.nx * 3 + k * data.nx * data.ny * 3] = static_cast<float>(data(k, j, i).ve[t][1]);
        data3[2 + i * 3 + j * data.nx * 3 + k * data.nx * data.ny * 3] = static_cast<float>(data(k, j, i).ve[t][2]);
      }
    }
  }

  ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
  ofs.write(reinterpret_cast<char *>(data3), allsize);
  delete[] data3;

  ofs.close();

  fp = fopen(file.c_str(), "a");
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

#endif // EXPORT_VTI_INL_H