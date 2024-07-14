/**
 * @file VTU.h
 * @author K.Ueda
 * @date Jun, 2024
*/

#include "VTK.h"

void VTK::exportMeshPartitionVTU(const std::string &file, Node &node, Cell &cell)
{
    FILE *fp;
    if((fp = fopen(file.c_str(), "w")) == NULL){
        std::cout << file << " open error" << std::endl;
        exit(1);
    }
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
    fprintf(fp, "<Points>\n");
    int offset = 0;
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int in=0; in<node.nNodesGlobal; in++){
        fprintf(fp, "%e %e %e\n", node.x[in][0], node.x[in][1], node.x[in][2]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
    for(int p=0; p<cell(ic).node.size(); p++) fprintf(fp, "%d ", cell(ic).node[p]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    int num = 0;
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        num += cell(ic).node.size();
        fprintf(fp, "%d\n", num);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++) fprintf(fp, "%d\n", cell(ic).cellType);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Cells>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"cellSubId\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        fprintf(fp,"%d \n", cell(ic).subId);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</CellData>\n");
    fprintf(fp, "<PointData>\n");
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"nodeSubId\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for(int in=0; in<node.nNodesGlobal; in++){
        fprintf(fp,"%d \n", node.subId[in]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</PointData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}


void VTK::exportPhiVTU(const std::string &file, Node &node, Cell &cell)
{
    FILE *fp;
    if((fp = fopen(file.c_str(), "w")) == NULL){
        std::cout << file << " open error" << std::endl;
        exit(1);
    }
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
    fprintf(fp, "<Points>\n");
    int offset = 0;
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int in=0; in<node.nNodesGlobal; in++){
        fprintf(fp, "%e %e %e\n", node.x[in][0], node.x[in][1], node.x[in][2]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
    for(int p=0; p<cell(ic).node.size(); p++) fprintf(fp, "%d ", cell(ic).node[p]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    int num = 0;
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        num += cell(ic).node.size();
        fprintf(fp, "%d\n", num);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++) fprintf(fp, "%d\n", cell(ic).cellType);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Cells>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        fprintf(fp,"%e \n", cell(ic).phi);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</CellData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}

void VTK::exportSolutionVTU(const std::string &file, Node &node, Cell &cell, DataType dataType)
{
    FILE *fp;
    if((fp = fopen(file.c_str(), "w")) == NULL){
        std::cout << file << " open error" << std::endl;
        exit(1);
    }
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
    fprintf(fp, "<Points>\n");
    int offset = 0;
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int in=0; in<node.nNodesGlobal; in++){
        fprintf(fp, "%e %e %e\n", node.x[in][0], node.x[in][1], node.x[in][2]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
    for(int p=0; p<cell(ic).node.size(); p++) fprintf(fp, "%d ", cell(ic).node[p]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    int num = 0;
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        num += cell(ic).node.size();
        fprintf(fp, "%d\n", num);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++) fprintf(fp, "%d\n", cell(ic).cellType);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Cells>\n");
    if(dataType == DataType::VELOCITY){
        fprintf(fp, "<PointData>\n");
        fprintf(fp, "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(int in=0; in<node.nNodesGlobal; in++){
            fprintf(fp,"%e %e %e \n", node.v[in][0], node.v[in][1], node.v[in][2]);
        }
        fprintf(fp, "</DataArray>\n");
        fprintf(fp, "</PointData>\n");
    }else if(dataType == DataType::PRESSURE){
        fprintf(fp, "<PointData>\n");
        fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for(int in=0; in<node.nNodesGlobal; in++){
            fprintf(fp,"%e \n", node.p[in]);
        }
        fprintf(fp, "</DataArray>\n");
        fprintf(fp, "</PointData>\n");
    }
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}


void VTK::exportMainVariablesVTU(const std::string &file, Node &node, Cell &cell, const int t, DataType dataType)
{
    FILE *fp;
    if((fp = fopen(file.c_str(), "w")) == NULL){
        std::cout << file << " open error" << std::endl;
        exit(1);
    }
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
    fprintf(fp, "<Points>\n");
    int offset = 0;
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int in=0; in<node.nNodesGlobal; in++){
        fprintf(fp, "%e %e %e\n", node.x[in][0], node.x[in][1], node.x[in][2]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
    for(int p=0; p<cell(ic).node.size(); p++) fprintf(fp, "%d ", cell(ic).node[p]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    int num = 0;
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        num += cell(ic).node.size();
        fprintf(fp, "%d\n", num);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++) fprintf(fp, "%d\n", cell(ic).cellType);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Cells>\n");
    if(dataType == DataType::VELOCITY){
        fprintf(fp, "<PointData>\n");
        fprintf(fp, "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(int in=0; in<node.nNodesGlobal; in++){
            fprintf(fp,"%e %e %e \n", node.vt[t][in][0], node.vt[t][in][1], node.vt[t][in][2]);
        }
        fprintf(fp, "</DataArray>\n");
        fprintf(fp, "</PointData>\n");
    }else if(dataType == DataType::PRESSURE){
        fprintf(fp, "<PointData>\n");
        fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for(int in=0; in<node.nNodesGlobal; in++){
            fprintf(fp,"%e \n", node.pt[t][in]);
        }
        fprintf(fp, "</DataArray>\n");
        fprintf(fp, "</PointData>\n");
    }
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}



void VTK::exportAdjointSolutionVTU(const std::string &file, Node &node, Cell &cell,
                                      DataType dataType)
{
    FILE *fp;
    if((fp = fopen(file.c_str(), "w")) == NULL){
        std::cout << file << " open error" << std::endl;
        exit(1);
    }
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
    fprintf(fp, "<Points>\n");
    int offset = 0;
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int in=0; in<node.nNodesGlobal; in++){
        fprintf(fp, "%e %e %e\n", node.x[in][0], node.x[in][1], node.x[in][2]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
    for(int p=0; p<cell(ic).node.size(); p++) fprintf(fp, "%d ", cell(ic).node[p]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    int num = 0;
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        num += cell(ic).node.size();
        fprintf(fp, "%d\n", num);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++) fprintf(fp, "%d\n", cell(ic).cellType);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Cells>\n");
    if(dataType == DataType::ADJOINT_W){
        fprintf(fp, "<PointData>\n");
        fprintf(fp, "<DataArray type=\"Float64\" Name=\"w\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(int in=0; in<node.nNodesGlobal; in++){
            fprintf(fp,"%e %e %e \n", node.w[in][0], node.w[in][1], node.w[in][2]);
        }
        fprintf(fp, "</DataArray>\n");
        fprintf(fp, "</PointData>\n");
    }else if(dataType == DataType::ADJOINT_Q){
        fprintf(fp, "<PointData>\n");
        fprintf(fp, "<DataArray type=\"Float64\" Name=\"q\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for(int in=0; in<node.nNodesGlobal; in++){
            fprintf(fp,"%e \n", node.q[in]);
        }
        fprintf(fp, "</DataArray>\n");
        fprintf(fp, "</PointData>\n");
    }else if(dataType == DataType::ADJOINT_L){
        fprintf(fp, "<PointData>\n");
        fprintf(fp, "<DataArray type=\"Float64\" Name=\"l\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for(int in=0; in<node.nNodesGlobal; in++){
            fprintf(fp,"%e %e %e \n", node.l[in][0], node.l[in][1], node.l[in][2]);
        }
        fprintf(fp, "</DataArray>\n");
        fprintf(fp, "</PointData>\n");
    }
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}

void VTK::exportFeedbackForceVTU(const std::string &file, Node &node, Cell &cell,
                                    std::vector<std::vector<double>> &feedback)
{
    FILE *fp;
    if((fp = fopen(file.c_str(), "w")) == NULL){
        std::cout << file << " open error" << std::endl;
        exit(1);
    }
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
    fprintf(fp, "<Points>\n");
    int offset = 0;
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int in=0; in<node.nNodesGlobal; in++){
        fprintf(fp, "%e %e %e\n", node.x[in][0], node.x[in][1], node.x[in][2]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
    for(int p=0; p<cell(ic).node.size(); p++) fprintf(fp, "%d ", cell(ic).node[p]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    int num = 0;
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        num += cell(ic).node.size();
        fprintf(fp, "%d\n", num);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++) fprintf(fp, "%d\n", cell(ic).cellType);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Cells>\n");
    fprintf(fp, "<PointData>\n");
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"feedback\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int in=0; in<node.nNodesGlobal; in++){
        fprintf(fp,"%e %e %e \n", feedback[in][0], feedback[in][1], feedback[in][2]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</PointData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}

void VTK::exportSnapShotVTU(const std::string &file, Node &node, Cell &cell, SnapShot &snap, const int &snapCount)
{
    FILE *fp;
    if((fp = fopen(file.c_str(), "w")) == NULL){
        std::cout << file << " open error" << std::endl;
        exit(1);
    }
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.nNodesGlobal, cell.nCellsGlobal);
    fprintf(fp, "<Points>\n");
    int offset = 0;
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int in=0; in<node.nNodesGlobal; in++){
        fprintf(fp, "%e %e %e\n", node.x[in][0], node.x[in][1], node.x[in][2]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Points>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
    for(int p=0; p<cell(ic).node.size(); p++) fprintf(fp, "%d ", cell(ic).node[p]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    int num = 0;
    for(int ic=0; ic<cell.nCellsGlobal; ic++){
        num += cell(ic).node.size();
        fprintf(fp, "%d\n", num);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for(int ic=0; ic<cell.nCellsGlobal; ic++) fprintf(fp, "%d\n", cell(ic).cellType);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Cells>\n");
    fprintf(fp, "<PointData>\n");
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int in=0; in<node.nNodesGlobal; in++){
        fprintf(fp,"%e %e %e \n", snap.v[snapCount][in][0], snap.v[snapCount][in][1], snap.v[snapCount][in][2]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</PointData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}