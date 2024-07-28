/**
 * @file BIN.cpp
 * @author K.Ueda
 * @date July, 2024
*/

#include "VTK.h"

void BIN::exportScalarDataBIN(const std::string &file, std::vector<double> &vec)
{
    std::ofstream ofs(file, std::ios::binary);
    if(!ofs) {
        std::cerr << "Could not open " << file << std::endl;
        return;
    }

    ofs.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(double));

    ofs.close();
    if(!ofs.good()) {
        std::cerr << "Error occurred at writing time." << std::endl;
    }
}

void BIN::exportVectorDataBIN(const std::string &file, std::vector<std::vector<double>> &vec)
{
    std::ofstream ofs(file, std::ios::binary);
    if(!ofs){
        std::cerr << "Could not open " << file << std::endl;
        return;
    }

    int rows = vec.size();
    ofs.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    
    for(const auto& row : vec){
        int cols = row.size();
        ofs.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
        ofs.write(reinterpret_cast<const char*>(row.data()), cols * sizeof(double));
    }

    ofs.close();
    if(!ofs.good()){
        std::cerr << "Error occurred at writing time." << std::endl;
    }
}

void BIN::importScalarDataBIN(const std::string &file, std::vector<double> &vec)
{
    std::ifstream ifs(file, std::ios::binary);
    if(!ifs){
        throw std::runtime_error("Couldn't open file: " + file);
    }

    ifs.seekg(0, std::ios::end);
    std::streamsize size = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    if(size % sizeof(double) != 0){
        throw std::runtime_error("File size error: " + file);
    }

    vec.resize(size / sizeof(double));

    if(!ifs.read(reinterpret_cast<char*>(vec.data()), size)){
        throw std::runtime_error("Couldn't read file: " + file);
    }

    ifs.close();
}

void BIN::importVectorDataBIN(const std::string &file, std::vector<std::vector<double>> &vec)
{
    std::ifstream ifs(file, std::ios::binary);
    if(!ifs){
        throw std::runtime_error("Couldn't open file: " + file);
    }

    int rows;
    if(!ifs.read(reinterpret_cast<char*>(&rows), sizeof(rows))) {
        throw std::runtime_error("Failed to read rows: " + file);
    }
 
    vec.resize(rows);

    for(int j=0; j<rows; j++){
        int cols;
        if(!ifs.read(reinterpret_cast<char*>(&cols), sizeof(cols))) {
            throw std::runtime_error("Failed to read cols: " + file);
        }
        vec[j].resize(cols);
        if(!ifs.read(reinterpret_cast<char*>(vec[j].data()), cols * sizeof(double))){
            throw std::runtime_error("Couldn't read file: " + file);
        }
    }

    ifs.close();
}