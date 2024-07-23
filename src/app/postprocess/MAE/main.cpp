#include "PostInverseProblem.h"
MyMPI mpi;

int main(int argc, char* argv[])
{
    std::string inputFile = argv[1]; 
    std::string appName = "MAE";
    Config conf(inputFile, appName);
    if(conf.isReadingError) return EXIT_FAILURE;

    std::string output = "output";
    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    conf.outputDir = "output/" + conf.outputDir;
    mkdir(conf.outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}