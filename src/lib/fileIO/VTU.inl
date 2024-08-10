/**
 * @file FileIO.h
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef VTK_INL_H
#define VTK_INL_H

#include "FileIO.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

template <typename T>
void VTKTMP::exportScalarPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<T> &p)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell(ic).node.size(); p++)
      fprintf(fp, "%d ", cell(ic).node[p]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    num += cell(ic).node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
    fprintf(fp, "%d\n", cell(ic).cellType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", dataName, offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal;
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  std::fstream ofs;
  ofs.open(file.c_str(), std::ios::out | std::ios::app | std::ios_base::binary);
  float *data_point_d3 = new float[node.nNodesGlobal * 3];
  float *data_point_d1 = new float[node.nNodesGlobal];

  num = 0;
  int size = 0;
  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d3[num] = (float)node.x[in][0];
    num++;
    data_point_d3[num] = (float)node.x[in][1];
    num++;
    data_point_d3[num] = (float)node.x[in][2];
    num++;
  }
  size = sizeof(float) * node.nNodesGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d3, size);

  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d1[in] = (float)p[in];
  }
  size = sizeof(float) * node.nNodesGlobal;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d1, size);

  delete[] data_point_d3;
  delete[] data_point_d1;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

template <typename T>
void VTKTMP::exportVectorPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<std::vector<T>> &p)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell(ic).node.size(); p++)
      fprintf(fp, "%d ", cell(ic).node[p]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    num += cell(ic).node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
    fprintf(fp, "%d\n", cell(ic).cellType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", dataName, offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal * 3;
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  std::fstream ofs;
  ofs.open(file.c_str(), std::ios::out | std::ios::app | std::ios_base::binary);
  float *data_point_d3 = new float[node.nNodesGlobal * 3];

  int size = 0;
  num = 0;
  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d3[num] = (float)node.x[in][0];
    num++;
    data_point_d3[num] = (float)node.x[in][1];
    num++;
    data_point_d3[num] = (float)node.x[in][2];
    num++;
  }
  size = sizeof(float) * node.nNodesGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d3, size);

  num = 0;
  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d3[num] = (float)p[in][0];
    num++;
    data_point_d3[num] = (float)p[in][1];
    num++;
    data_point_d3[num] = (float)p[in][2];
    num++;
  }
  size = sizeof(float) * node.nNodesGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d3, size);

  delete[] data_point_d3;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

template <typename T>
void VTKTMP::exportScalarCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<T> &c)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell(ic).node.size(); p++)
      fprintf(fp, "%d ", cell(ic).node[p]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    num += cell(ic).node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
    fprintf(fp, "%d\n", cell(ic).cellType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", dataName, offset);
  offset += sizeof(int) + sizeof(float) * cell.nCellsGlobal;
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  std::fstream ofs;
  ofs.open(file.c_str(), std::ios::out | std::ios::app | std::ios_base::binary);
  float *data_point_d3 = new float[node.nNodesGlobal * 3];
  float *data_cell_d1 = new float[cell.nCellsGlobal];

  num = 0;
  int size = 0;
  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d3[num] = (float)node.x[in][0];
    num++;
    data_point_d3[num] = (float)node.x[in][1];
    num++;
    data_point_d3[num] = (float)node.x[in][2];
    num++;
  }
  size = sizeof(float) * node.nNodesGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d3, size);

  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    data_cell_d1[ic] = (float)c[ic];
  }
  size = sizeof(float) * cell.nCellsGlobal;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_cell_d1, size);

  delete[] data_point_d3;
  delete[] data_cell_d1;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

template <typename T>
void VTKTMP::exportVectorCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, std::vector<std::vector<T>> &c)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell(ic).node.size(); p++)
      fprintf(fp, "%d ", cell(ic).node[p]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    num += cell(ic).node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
    fprintf(fp, "%d\n", cell(ic).cellType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", dataName, offset);
  offset += sizeof(int) + sizeof(float) * cell.nCellsGlobal * 3;
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  std::fstream ofs;
  ofs.open(file.c_str(), std::ios::out | std::ios::app | std::ios_base::binary);
  float *data_point_d3 = new float[node.nNodesGlobal * 3];
  float *data_cell_d3 = new float[cell.nCellsGlobal * 3];

  int size = 0;
  num = 0;
  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d3[num] = (float)node.x[in][0];
    num++;
    data_point_d3[num] = (float)node.x[in][1];
    num++;
    data_point_d3[num] = (float)node.x[in][2];
    num++;
  }
  size = sizeof(float) * node.nNodesGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d3, size);

  num = 0;
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    data_cell_d3[num] = (float)c[ic][0];
    num++;
    data_cell_d3[num] = (float)c[ic][1];
    num++;
    data_cell_d3[num] = (float)c[ic][2];
    num++;
  }
  size = sizeof(float) * cell.nCellsGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_cell_d3, size);

  delete[] data_point_d3;
  delete[] data_cell_d3;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

template <typename T>
void VTKTMP::exportScalarPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array1D<T> &p)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell(ic).node.size(); p++)
      fprintf(fp, "%d ", cell(ic).node[p]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    num += cell(ic).node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
    fprintf(fp, "%d\n", cell(ic).cellType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", dataName, offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal;
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  std::fstream ofs;
  ofs.open(file.c_str(), std::ios::out | std::ios::app | std::ios_base::binary);
  float *data_point_d3 = new float[node.nNodesGlobal * 3];
  float *data_point_d1 = new float[node.nNodesGlobal];

  num = 0;
  int size = 0;
  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d3[num] = (float)node.x[in][0];
    num++;
    data_point_d3[num] = (float)node.x[in][1];
    num++;
    data_point_d3[num] = (float)node.x[in][2];
    num++;
  }
  size = sizeof(float) * node.nNodesGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d3, size);

  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d1[in] = (float)p(in);
  }
  size = sizeof(float) * node.nNodesGlobal;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d1, size);

  delete[] data_point_d3;
  delete[] data_point_d1;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

template <typename T>
void VTKTMP::exportVectorPointDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array2D<T> &p)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell(ic).node.size(); p++)
      fprintf(fp, "%d ", cell(ic).node[p]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    num += cell(ic).node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
    fprintf(fp, "%d\n", cell(ic).cellType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", dataName, offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal * 3;
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  std::fstream ofs;
  ofs.open(file.c_str(), std::ios::out | std::ios::app | std::ios_base::binary);
  float *data_point_d3 = new float[node.nNodesGlobal * 3];

  int size = 0;
  num = 0;
  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d3[num] = (float)node.x[in][0];
    num++;
    data_point_d3[num] = (float)node.x[in][1];
    num++;
    data_point_d3[num] = (float)node.x[in][2];
    num++;
  }
  size = sizeof(float) * node.nNodesGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d3, size);

  num = 0;
  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d3[num] = (float)p(in, 0);
    num++;
    data_point_d3[num] = (float)p(in, 1);
    num++;
    data_point_d3[num] = (float)p(in, 2);
    num++;
  }
  size = sizeof(float) * node.nNodesGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d3, size);

  delete[] data_point_d3;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

template <typename T>
void VTKTMP::exportScalarCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array1D<T> &c)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell(ic).node.size(); p++)
      fprintf(fp, "%d ", cell(ic).node[p]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    num += cell(ic).node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
    fprintf(fp, "%d\n", cell(ic).cellType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", dataName, offset);
  offset += sizeof(int) + sizeof(float) * cell.nCellsGlobal;
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  std::fstream ofs;
  ofs.open(file.c_str(), std::ios::out | std::ios::app | std::ios_base::binary);
  float *data_point_d3 = new float[node.nNodesGlobal * 3];
  float *data_cell_d1 = new float[cell.nCellsGlobal];

  num = 0;
  int size = 0;
  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d3[num] = (float)node.x[in][0];
    num++;
    data_point_d3[num] = (float)node.x[in][1];
    num++;
    data_point_d3[num] = (float)node.x[in][2];
    num++;
  }
  size = sizeof(float) * node.nNodesGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d3, size);

  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    data_cell_d1[ic] = (float)c(ic);
  }
  size = sizeof(float) * cell.nCellsGlobal;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_cell_d1, size);

  delete[] data_point_d3;
  delete[] data_cell_d1;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

template <typename T>
void VTKTMP::exportVectorCellDataVTU(const std::string &file, const char *dataName, Node &node, Cell &cell, Array2D<T> &c)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(float) * node.nNodesGlobal * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    for (int p = 0; p < cell(ic).node.size(); p++)
      fprintf(fp, "%d ", cell(ic).node[p]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    num += cell(ic).node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
    fprintf(fp, "%d\n", cell(ic).cellType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", dataName, offset);
  offset += sizeof(int) + sizeof(float) * cell.nCellsGlobal * 3;
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  std::fstream ofs;
  ofs.open(file.c_str(), std::ios::out | std::ios::app | std::ios_base::binary);
  float *data_point_d3 = new float[node.nNodesGlobal * 3];
  float *data_cell_d3 = new float[cell.nCellsGlobal * 3];

  int size = 0;
  num = 0;
  for (int in = 0; in < node.nNodesGlobal; in++)
  {
    data_point_d3[num] = (float)node.x[in][0];
    num++;
    data_point_d3[num] = (float)node.x[in][1];
    num++;
    data_point_d3[num] = (float)node.x[in][2];
    num++;
  }
  size = sizeof(float) * node.nNodesGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_d3, size);

  num = 0;
  for (int ic = 0; ic < cell.nCellsGlobal; ic++)
  {
    data_cell_d3[num] = (float)c(ic, 0);
    num++;
    data_cell_d3[num] = (float)c(ic, 1);
    num++;
    data_cell_d3[num] = (float)c(ic, 2);
    num++;
  }
  size = sizeof(float) * cell.nCellsGlobal * 3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_cell_d3, size);

  delete[] data_point_d3;
  delete[] data_cell_d3;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    std::cerr << file << " open error" << std::endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

#endif // VTK_INL_H
