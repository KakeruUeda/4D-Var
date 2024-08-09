#include "PostInverseProblem.h"
MyMPI mpi;

int main(int argc, char *argv[])
{
    std::string inputFile = argv[1];
    std::string appName = "FLOWRATE";
    Config conf(inputFile, appName);
    if (conf.isReadingError)
        return EXIT_FAILURE;

    PostInverseProblem post;
    post.initialize(conf);

    std::string output = "output";
    mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    post.outputDir = "output/" + post.outputDir;
    mkdir(post.outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    std::vector<double> flowRateRefVec(post.nRef, 0e0);
    std::vector<double> flowRateOptVec(post.nRef, 0e0);

    ofstream outFlowRateVelRef(post.outputDir + "/flowRateVelRef.dat");
    for (int t = 0; t < post.nRef; t++)
    {
        double flowRate = post.compFlowRate(post.velRef[t]);
        outFlowRateVelRef << t << " " << flowRate << std::endl;
        flowRateRefVec[t] = flowRate;
    }
    outFlowRateVelRef.close();

    ofstream outFlowRateVelOpt(post.outputDir + "/flowRateVelOpt.dat");
    for (int t = 0; t < post.nRef; t++)
    {
        double flowRate = post.compFlowRate(post.velOpt[t]);
        outFlowRateVelOpt << t << " " << flowRate << std::endl;
        flowRateOptVec[t] = flowRate;
    }
    outFlowRateVelOpt.close();

    double meanFlowRateError;

    for (int t = 0; t < post.nRef; t++)
    {
        meanFlowRateError += post.compFlowRateError(flowRateRefVec[t], flowRateOptVec[t]);
    }
    meanFlowRateError /= post.nRef;

    ofstream outMeanFlowRateError(post.outputDir + "/meanFlowRateError.dat");
    outMeanFlowRateError << meanFlowRateError << std::endl;
    outMeanFlowRateError.close();

    return EXIT_SUCCESS;
}