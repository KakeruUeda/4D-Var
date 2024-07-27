/**
 * @file VTI.cpp
 * @author K.Ueda
 * @date Jun, 2024
*/

#include "VTK.h"


void VTK::exportScalarPointDataVTI(const std::string &file, const char *dataName,
                                   std::vector<double> &p,
                                   const int nx,const int ny,const int nz,
                                   const double dx,const double dy,const double dz)
{
    FILE *fp;
    fp = fopen(file.c_str(), "w");

    if(fp == NULL){
        cout << file << " open error" << endl;
        exit(1);
    }
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n",0,nx,0,ny,0,nz,0e0,0e0,0e0,dx,dy,dz);  
    fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n",0,nx,0,ny,0,nz);  
    fprintf(fp, "<PointData>\n");
    fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"0\">\n",dataName);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</PointData>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "</CellData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</ImageData>\n");
    fprintf(fp,"<AppendedData encoding=\"raw\">\n");
    fprintf(fp,"_");
    fclose(fp);

    unsigned long allsize=(nx+1)*(ny+1)*(nz+1)*sizeof(float);

    float *data = new float[(nx+1)*(ny+1)*(nz+1)];
    for(int k=0; k<nz+1; k++){
        for(int j=0; j<ny+1; j++){
            for(int i=0; i<nx+1; i++){
                int n = i + j * (nx+1) + k * (nx+1) * (ny+1);
                data[n] = (float)p[n];  
            }
        }
    }

  fstream ofs;
  ofs.open(file.c_str(),ios::out | ios::app | ios_base::binary);
  ofs.write((char *)&allsize,sizeof(allsize));
  ofs.write((char *)data,allsize);
  ofs.close();

  free(data);

  fp = fopen(file.c_str(), "a");
  fprintf(fp,"\n");
  fprintf(fp,"</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

}

void VTK::exportVectorPointDataVTI(const std::string &file, const char *dataName, 
                                   std::vector<std::vector<double>> &p,
                                   const int nx, const int ny, const int nz,
                                   const double dx,const double dy,const double dz)
{
    FILE *fp;
    fp = fopen(file.c_str(), "w");

    if(fp == NULL){
        cout << file << " open error" << endl;
        exit(1);
    }
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n",0,nx,0,ny,0,nz,0e0,0e0,0e0,dx,dy,dz);  
    fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n",0,nx,0,ny,0,nz);  
    fprintf(fp, "<PointData>\n");
    fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\">\n", dataName);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</PointData>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "</CellData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</ImageData>\n");
    fprintf(fp,"<AppendedData encoding=\"raw\">\n");
    fprintf(fp,"_");
    fclose(fp);

    unsigned long allsize = (nx+1)*(ny+1)*(nz+1)*3*sizeof(float);

    float *data = new float[(nx+1)*(ny+1)*(nz+1)*3];
    for(int k=0; k<nz+1; k++){
        for(int j=0; j<ny+1; j++){
            for(int i=0; i<nx+1; i++){
                int n = i + j * (nx+1) + k * (nx+1) * (ny+1);
                data[0+i*3+j*(nx+1)*3+k*(nx+1)*(ny+1)*3] = (float)p[n][0];
                data[1+i*3+j*(nx+1)*3+k*(nx+1)*(ny+1)*3] = (float)p[n][1];
                data[2+i*3+j*(nx+1)*3+k*(nx+1)*(ny+1)*3] = (float)p[n][2];
            }
        }
    }

    fstream ofs;
    ofs.open(file.c_str(),ios::out | ios::app | ios_base::binary);
    ofs.write((char *)&allsize,sizeof(allsize));
    ofs.write((char *)data,allsize);
    ofs.close();

    free(data);

    fp = fopen(file.c_str(), "a");
    fprintf(fp,"\n");
    fprintf(fp,"</AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);

}

void VTK::exportScalarCellDataVTI(const std::string &file, const char *dataName, std::vector<double> &c,
                                  const int nx, const int ny, const int nz,
                                  const double dx, const double dy, const double dz)
{
    FILE *fp;
    fp = fopen(file.c_str(), "w");

    if(fp == NULL){
        cout << file << " open error" << endl;
        exit(1);
    }
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n",0,nx,0,ny,0,nz,0e0,0e0,0e0,dx,dy,dz);  
    fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n",0,nx,0,ny,0,nz);  
    fprintf(fp, "<PointData>\n");
    fprintf(fp, "</PointData>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"0\">\n",dataName);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</CellData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</ImageData>\n");
    fprintf(fp,"<AppendedData encoding=\"raw\">\n");
    fprintf(fp,"_");
    fclose(fp);

    float *data = new float[nx*ny*nz];
    unsigned long allsize = nx*ny*nz*sizeof(float);
    for(int k=0; k<nz; k++){
        for(int j=0; j<ny; j++){
            for(int i=0; i<nx; i++){ 
                int n = i + j * nx + k * nx * ny;
                data[n] = (float)c[n];
            }
        }
    }

    fstream ofs;
    ofs.open(file.c_str(),ios::out | ios::app | ios_base::binary);
    ofs.write((char *)&allsize,sizeof(allsize));
    ofs.write((char *)data,allsize);
    ofs.close();

    free(data);

    fp = fopen(file.c_str(), "a");
    fprintf(fp,"\n");
    fprintf(fp,"</AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);

}

void VTK::exportVectorCellDataVTI(const std::string &file, const char *dataName, 
                                 std::vector<std::vector<double>> &c,
                                 const int nx,const int ny,const int nz,
                                 const double dx,const double dy,const double dz)
{
    FILE *fp;
    fp = fopen(file.c_str(), "w");

    if(fp == NULL){
        cout << file << " open error" << endl;
        exit(1);
    }
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n",0,nx,0,ny,0,nz,0e0,0e0,0e0,dx,dy,dz);  
    fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n",0,nx,0,ny,0,nz);  
    fprintf(fp, "<PointData>\n");
    fprintf(fp, "</PointData>\n");
    fprintf(fp, "<CellData>\n");
    fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\">\n",dataName);
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</CellData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</ImageData>\n");
    fprintf(fp,"<AppendedData encoding=\"raw\">\n");
    fprintf(fp,"_");
    fclose(fp);

    float *data = new float[nx*ny*nz*3];
    unsigned long allsize = nx*ny*nz*3*sizeof(float);
    for(int k=0; k<nz; k++){
        for(int j=0; j<ny; j++){
            for(int i=0; i<nx; i++){
                int n = i + j * nx + k * nx * ny;
                data[0+i*3+j*nx*3+k*nx*ny*3] = (float)c[n][0];
                data[1+i*3+j*nx*3+k*nx*ny*3] = (float)c[n][1];
                data[2+i*3+j*nx*3+k*nx*ny*3] = (float)c[n][2];
            }
        }
    }

    fstream ofs;
    ofs.open(file.c_str(),ios::out | ios::app | ios_base::binary);
    ofs.write((char *)&allsize,sizeof(allsize));
    ofs.write((char *)data,allsize);
    ofs.close();

    free(data);

    fp = fopen(file.c_str(), "a");
    fprintf(fp,"\n");
    fprintf(fp,"</AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);

}

void VTK::exportVelocityDataVTI(const std::string &file, DataGrid &data, const int t)
{
    FILE *fp = fopen(file.c_str(), "w");

    if (fp == NULL) {
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
    for (int k=0; k<data.nz; k++) {
        for (int j=0; j<data.ny; j++) {
            for (int i=0; i<data.nx; i++) {
                int index = i + j * data.nx + k * data.nx * data.ny;
                data1[0+i*3+j*data.nx*3+k*data.nx*data.ny*3] = static_cast<float>(data(k, j, i).vCFD[t][0]);
                data1[1+i*3+j*data.nx*3+k*data.nx*data.ny*3] = static_cast<float>(data(k, j, i).vCFD[t][1]);
                data1[2+i*3+j*data.nx*3+k*data.nx*data.ny*3] = static_cast<float>(data(k, j, i).vCFD[t][2]);
            }
        }
    }

    ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
    ofs.write(reinterpret_cast<char *>(data1), allsize);
    delete[] data1;

    // Write vMRI data
    float *data2 = new float[data.nx * data.ny * data.nz * 3];
    for (int k=0; k<data.nz; k++) {
        for (int j=0; j<data.ny; j++) {
            for (int i=0; i<data.nx; i++) {
                int index = i + j * data.nx + k * data.nx * data.ny;
                data2[0+i*3+j*data.nx*3+k*data.nx*data.ny*3] = static_cast<float>(data(k, j, i).vMRI[t][0]);
                data2[1+i*3+j*data.nx*3+k*data.nx*data.ny*3] = static_cast<float>(data(k, j, i).vMRI[t][1]);
                data2[2+i*3+j*data.nx*3+k*data.nx*data.ny*3] = static_cast<float>(data(k, j, i).vMRI[t][2]);
            }
        }
    }

    ofs.write(reinterpret_cast<char *>(&allsize), sizeof(allsize));
    ofs.write(reinterpret_cast<char *>(data2), allsize);
    delete[] data2;

    // Write ve data
    float *data3 = new float[data.nx * data.ny * data.nz * 3];
    for (int k=0; k<data.nz; k++) {
        for (int j=0; j<data.ny; j++) {
            for (int i=0; i<data.nx; i++) {
                int index = i + j * data.nx + k * data.nx * data.ny;
                data3[0+i*3+j*data.nx*3+k*data.nx*data.ny*3] = static_cast<float>(data(k, j, i).ve[t][0]);
                data3[1+i*3+j*data.nx*3+k*data.nx*data.ny*3] = static_cast<float>(data(k, j, i).ve[t][1]);
                data3[2+i*3+j*data.nx*3+k*data.nx*data.ny*3] = static_cast<float>(data(k, j, i).ve[t][2]);
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

