#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

void addGaussianNoise(vector<double> &u, vector<double> &v, vector<double> &w, vector<double> &uNoise,
                      vector<double> &vNoise, vector<double> &wNoise, double snr_u, double snr_v, double snr_w, double sigma, vector<double> &mask)
{
  auto maxVal = [](const vector<double> &vec) { return *max_element(vec.begin(), vec.end()); };

  auto aveVal = [](const vector<double> &vec, const vector<double> &mask) {
    double value = 0e0;
    size_t count = 0;

    for(int k = 0; k < 27; k++) {
      for(int j = 0; j < 27; j++) {
        for(int i = 0; i < 27; i++) {
          int ic = i + j * 27 + k * 27 * 27;
          if(mask[ic] < 1e-8) continue;
          value += fabs(vec[ic]);
          count++;
        }
      }
    }
    return value / count;
    //return accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
  };

  double uSigma = aveVal(u, mask) / snr_u;
  double vSigma = aveVal(v, mask) / snr_v;
  double wSigma = aveVal(w, mask) / snr_w;

  double sigma_ave = (uSigma + vSigma + wSigma) / 3e0;

  std::cout << uSigma << " " << vSigma << " " << wSigma << " " << sigma_ave << std::endl;

  random_device rd;
  mt19937 gen(rd());

  normal_distribution<> du(0, uSigma);
  normal_distribution<> dv(0, vSigma);
  normal_distribution<> dw(0, wSigma);

  for(size_t i = 0; i < u.size(); ++i) {
    uNoise[i] = u[i] + du(gen);
    vNoise[i] = v[i] + dv(gen);
    wNoise[i] = w[i] + dw(gen);
  }

  for(int k = 0; k < 27; k++) {
    for(int j = 0; j < 27; j++) {
      for(int i = 0; i < 27; i++) {
        int ic = i + j * 27 + k * 27 * 27;
        if((fabs(u[ic]) < 1e-8) && (fabs(v[ic]) < 1e-8) && (fabs(w[ic]) < 1e-8)) {
            uNoise[ic] = vNoise[ic] = wNoise[ic] = 0e0;
        }
      }
    }
  }
}

void inputData(const string &inputFile, vector<double> &u, vector<double> &v, vector<double> &w, int &nx, int &ny,
               int &nz, double &lx, double &ly, double &lz, double &dx, double &dy, double &dz)
{
  ifstream in(inputFile);
  if(!in) {
    cerr << "Error opening file: " << inputFile << endl;
    exit(1);
  }

  size_t totalSize = nx * ny * nz;

  u.resize(totalSize, 0.0);
  v.resize(totalSize, 0.0);
  w.resize(totalSize, 0.0);

  for(size_t i = 0; i < totalSize; ++i) {
    in >> u[i] >> v[i] >> w[i];
  }
}

void export_vti(const string &outputFile, const vector<double> &u, const vector<double> &v, const vector<double> &w,
                const vector<double> &uNoise, const vector<double> &vNoise, const vector<double> &wNoise, int nx,
                int ny, int nz, double dx, double dy, double dz)
{
  ofstream out(outputFile);
  if(!out) {
    cerr << "Error opening file: " << outputFile << endl;
    exit(1);
  }

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  out << "<ImageData WholeExtent=\"0 " << nx << " 0 " << ny << " 0 " << nz << "\" Origin=\"0 0 0\" Spacing=\"" << dx
      << " " << dy << " " << dz << "\">\n";
  out << "<Piece Extent=\"0 " << nx << " 0 " << ny << " 0 " << nz << "\">\n";
  out << "<PointData/>";
  out << "<CellData>\n";

  auto writeDataArray = [&](const string &name, const vector<double> &x, const vector<double> &y,
                            const vector<double> &z) {
    out << "<DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(size_t i = 0; i < x.size(); ++i) {
      out << x[i] << " " << y[i] << " " << z[i] << "\n";
    }
    out << "</DataArray>\n";
  };

  writeDataArray("original velocity", u, v, w);
  writeDataArray("noisy velocity", uNoise, vNoise, wNoise);

  out << "</CellData>\n";
  out << "</Piece>\n";
  out << "</ImageData>\n";
  out << "</VTKFile>\n";
}

void processTimeSeries(const string &dataDir, const string &snrFile_u, const string &snrFile_v, const string &snrFile_w,
                       const string &sigmaFile, const string &maskFile, const string &outputDir_vtk, const string &outputDir_dat)
{
  ifstream snrStream_u(snrFile_u);
  if(!snrStream_u) {
    cerr << "Error opening SNR file: " << snrFile_u << endl;
    exit(1);
  }

  ifstream snrStream_v(snrFile_v);
  if(!snrStream_v) {
    cerr << "Error opening SNR file: " << snrFile_v << endl;
    exit(1);
  }

  ifstream snrStream_w(snrFile_w);
  if(!snrStream_w) {
    cerr << "Error opening SNR file: " << snrFile_w << endl;
    exit(1);
  }

  ifstream sigmaStream(sigmaFile);
  if(!sigmaStream) {
    cerr << "Error opening sigma file: " << sigmaFile << endl;
    exit(1);
  }

  ifstream maskStream(maskFile);
  if(!maskStream) {
    cerr << "Error opening sigma file: " << maskFile << endl;
    exit(1);
  }

  double value;

  vector<double> snr_u;
  while(snrStream_u >> value) {
    snr_u.push_back(value);
  }
  vector<double> snr_v;
  while(snrStream_v >> value) {
    snr_v.push_back(value);
  }
  vector<double> snr_w;
  while(snrStream_w >> value) {
    snr_w.push_back(value);
  }
  vector<double> sigma;
  while(sigmaStream >> value) {
    sigma.push_back(value);
  }
  vector<double> mask;
  while(maskStream >> value) {
    mask.push_back(value);
  }
  vector<string> dataFiles;
  for(int t = 0; t < snr_u.size(); t++) {
    std::string file = dataDir + "/data_" + std::to_string(t) + ".dat";
    dataFiles.push_back(file);
  }

  //sort(dataFiles.begin(), dataFiles.end());

  if(!(snr_u.size() == snr_v.size() && snr_v.size() == snr_w.size() && snr_w.size() == dataFiles.size())) {
    cerr << "Mismatch between number of data files and SNR values." << endl;
    cerr << snr_u.size() << " " << snr_v.size() << " " << snr_w.size() << " " << dataFiles.size() << std::endl;
    exit(1);
  }

  for(size_t t = 0; t < dataFiles.size(); ++t) {
    string inputFile = dataFiles[t];
    stringstream vtiFile, noisyDataFile;
    vtiFile << outputDir_vtk << "/data_" << t << ".vti";
    noisyDataFile << outputDir_dat << "/data_" << t << ".dat";

    int nx, ny, nz;
    double lx, ly, lz, dx, dy, dz;

    vector<double> u, v, w;
    vector<double> uNoise, vNoise, wNoise;

    nx = ny = nz = 27;
    lx = ly = 0.050625;
    lz = 0.0513;

    dx = lx / nx;
    dy = ly / ny;
    dz = lz / nz;

    inputData(inputFile, u, v, w, nx, ny, nz, lx, ly, lz, dx, dy, dz);

    size_t totalSize = nx * ny * nz;
    uNoise.resize(totalSize, 0.0);
    vNoise.resize(totalSize, 0.0);
    wNoise.resize(totalSize, 0.0);

    addGaussianNoise(u, v, w, uNoise, vNoise, wNoise, snr_u[t], snr_v[t], snr_w[t], sigma[t], mask);
    export_vti(vtiFile.str(), u, v, w, uNoise, vNoise, wNoise, nx, ny, nz, dx, dy, dz);

    ofstream out(noisyDataFile.str());
    if(!out) {
      cerr << "Error opening file for writing: " << noisyDataFile.str() << endl;
      exit(1);
    }

    for(int k = 0; k < nz; ++k) {
      for(int j = 0; j < ny; ++j) {
        for(int i = 0; i < nx; ++i) {
          size_t idx = k * ny * nx + j * nx + i;
          out << uNoise[idx] << " " << vNoise[idx] << " " << wNoise[idx] << "\n";
        }
      }
    }
  }
}

int main()
{
  string dataDir = "../../../direct/voxelDataCreation/output/Ubend_test/data";
  string snrFile_u = "snr_u.dat";
  string snrFile_v = "snr_v.dat";
  string snrFile_w = "snr_w.dat";
  string sigmaFile = "sigma.dat";
  string maskFile = "mask.dat";
  string outputDir_vtk = "output_vtk";
  string outputDir_dat = "output_dat";

  fs::create_directory(outputDir_vtk);
  fs::create_directory(outputDir_dat);

  processTimeSeries(dataDir, snrFile_u, snrFile_v, snrFile_w, sigmaFile, maskFile, outputDir_vtk, outputDir_dat);

  return 0;
}
