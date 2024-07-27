/**
 * @file BIN.cpp
 * @author K.Ueda
 * @date July, 2024
*/

#include "VTK.h"

void BIN::exportScalarDataBIN(const std::string &file, std::vector<double> &vec)
{
    std::ofstream ofs(file, std::ios::binary);
    if (!ofs) {
        std::cerr << "Could not open " << file << std::endl;
        return;
    }

    ofs.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(double));

    ofs.close();
    if (!ofs.good()) {
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

    size_t rows = vec.size();
    ofs.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    
    for(const auto& row : vec){
        size_t cols = row.size();
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
        std::cerr << "Couldn't open file: " << file << std::endl;
        return;
    }

    ifs.seekg(0, std::ios::end);
    std::streamsize size = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    if(size % sizeof(double) != 0){
        std::cerr << "File size error: " << file << std::endl;
        return;
    }

    vec.resize(size / sizeof(double));

    if(!ifs.read(reinterpret_cast<char*>(vec.data()), size)){
        std::cerr << "Couldn't read file: " << file << std::endl;
    }

    ifs.close();
}

void BIN::importVectorDataBIN(const std::string &file, std::vector<std::vector<double>> &vec)
{
    std::ifstream ifs(file, std::ios::binary);
    if(!ifs){
        std::cerr << "Couldn't open file: " << file << std::endl;
        return;
    }

    size_t rows;
    if(!ifs.read(reinterpret_cast<char*>(&rows), sizeof(rows))) {
        std::cerr << "Failed to read rows: " << file << std::endl;
        return;
    }

    vec.resize(rows);

    for(size_t j=0; j<rows; j++){
        size_t cols;
        if(!ifs.read(reinterpret_cast<char*>(&cols), sizeof(cols))) {
            std::cerr << "Failed to read cols: " << file << std::endl;
            return;
        }
        vec[j].resize(cols);
        if(!ifs.read(reinterpret_cast<char*>(vec[j].data()), cols * sizeof(double))){
            std::cerr << "Couldn't read file: " << file << std::endl;
            return;
        }
    }

    ifs.close();
}