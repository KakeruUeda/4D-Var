/**
 * @file DAT.inl
 * @author K.Ueda
 * @date August, 2024
 */

#ifndef EXPORT_VTU_H
#define EXPORT_VTU_H

#include "Export.h"

void EXPORT::exportMeshPartitionVTU(const std::string &file, Node &node, Cell &cell)
{
  FILE *fp;
  if((fp = fopen(file.c_str(), "w")) == NULL) {
    cout << file << " open error" << endl;
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
  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    for(int p = 0; p < cell(ic).node.size(); p++)
      fprintf(fp, "%d ", cell(ic).node[p]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    num += cell(ic).node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(int ic = 0; ic < cell.nCellsGlobal; ic++)
    fprintf(fp, "%d\n", cell(ic).cellType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"subNodeId\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(int) * node.nNodesGlobal;
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  // This must be Int32.
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"subCellId\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(int) * cell.nCellsGlobal;
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  float *data_point_d3 = new float[node.nNodesGlobal * 3];
  int *data_point_i1 = new int[node.nNodesGlobal];
  int *data_cell_i1 = new int[cell.nCellsGlobal];

  num = 0;
  int size = 0;
  for(int in = 0; in < node.nNodesGlobal; in++) {
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

  for(int in = 0; in < node.nNodesGlobal; in++) {
    data_point_i1[in] = node.subId[in];
  }
  size = sizeof(int) * node.nNodesGlobal;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_point_i1, size);

  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    data_cell_i1[ic] = cell(ic).subId;
  }
  size = sizeof(int) * cell.nCellsGlobal;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_cell_i1, size);

  delete[] data_point_d3;
  delete[] data_point_i1;
  delete[] data_cell_i1;

  ofs.close();

  if((fp = fopen(file.c_str(), "a")) == NULL) {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

void EXPORT::exportPhiVTU(const std::string &file, Node &node, Cell &cell)
{
  FILE *fp;
  if((fp = fopen(file.c_str(), "w")) == NULL) {
    cout << file << " open error" << endl;
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
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    for(int p = 0; p < cell(ic).node.size(); p++)
      fprintf(fp, "%d ", cell(ic).node[p]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    num += cell(ic).node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(int ic = 0; ic < cell.nCellsGlobal; ic++)
    fprintf(fp, "%d\n", cell(ic).cellType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  // This must be Int32.
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"phi\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(float) * cell.nCellsGlobal;
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  float *data_point_d3 = new float[node.nNodesGlobal * 3];
  float *data_cell_d1 = new float[cell.nCellsGlobal];

  num = 0;
  int size = 0;
  for(int in = 0; in < node.nNodesGlobal; in++) {
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

  for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
    data_cell_d1[ic] = (float)cell(ic).phi;
  }
  size = sizeof(float) * cell.nCellsGlobal;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_cell_d1, size);

  delete[] data_point_d3;
  delete[] data_cell_d1;

  ofs.close();

  if((fp = fopen(file.c_str(), "a")) == NULL) {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

#endif  // EXPORT_VTU_INL_H
