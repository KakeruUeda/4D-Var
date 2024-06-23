#include "Config.h"

Config::Config(std::string inputFile, std::string appName)
{
    setApplication(appName);
    tryOpenConfigFile(inputFile);
    tryReadConfigFile();
}

void Config::setApplication(std::string appName)
{   
    if     (appName == "STRGRID")  app = Application::STRGRID;
    else if(appName == "SNS")      app = Application::SNS;
    else if(appName == "USNS")     app = Application::USNS;
    else if(appName == "TDVAR")    app = Application::TDVAR;
    else if(appName == "FDVAR")    app = Application::FDVAR;
    else
        if(mpi.myId == 0)
            std::cout << "Unknown appName" << std::endl;

    return;
}

void Config::tryOpenConfigFile(std::string inputFile)
{
    try{
        int error;
        if ((error = tp.read(inputFile)) != TP_NO_ERROR)
            throw std::runtime_error("Open error");
    }
    catch (const std::runtime_error& e){
		if(mpi.myId == 0) 
            std::cout << e.what() << std::endl;
		if(mpi.myId == 0) 
            isReadingError = true;
    }
}

void Config::tryReadConfigFile()
{
    try{
        readConfigFile();
    }
    catch(const std::runtime_error& e){
        if(mpi.myId == 0)
            std::cout << e.what() << std::endl;
        if(mpi.myId == 0)
            isReadingError = true;
    }
}

void Config::readConfigFile()
{
    switch(app){
        case Application::STRGRID:
            readBasicParameter(); 
            readStructuredGridParameter(); 
            readStructuredBoundaryParameter();
            break;

        case Application::USNS:
            readBasicParameter();
            readGridParameter();
            readTimeParameter();
            readPysicalParameter();
            readDarcyParameter();
            readPostprocessParameter();
            break;

        default:
            if(mpi.myId == 0)
                throw std::runtime_error("Unknown Application");
            break;
    }
}

void Config::setSolidBoundary()
{
    std::vector<double> vecTmp;
    vecTmp.resize(3, 0e0);

    for(int ic=0; ic<nCellsGlobal; ic++){
        if(phi[ic] < 1e-12){
            for(int p=0; p<nNodesInCell; p++){
                std::vector<double> vecTmp(dim, 0e0);
                vDirichlet[cell[ic][p]] = vecTmp;
            }
        }
    }
}

void Config::setFluidDomain()
{
    int nCellsGlobalTmp = nCellsGlobal;
    std::vector<std::vector<int>> cellTmp = cell;

    nCellsGlobal = 0;
    cell.erase(cell.begin(), cell.end());
    for(int ic=0; ic<nCellsGlobalTmp; ic++){
        if(phi[ic] < 1e-12) continue;
        sortCell.push_back(ic);
        cell.push_back(cellTmp[ic]);
        nCellsGlobal++;
    }

    std::vector<double> phiTmp = phi;
    phi.erase(phi.begin(), phi.end());

    int count = 0;
    for(int ic=0; ic<nCellsGlobalTmp; ic++){
        if(phiTmp[ic] < 1e-12) continue;
        phi.push_back(phiTmp[ic]);
    }

    sortNode.resize(nCellsGlobal * nNodesInCell);

    count = 0;
    for(int ic=0; ic<nCellsGlobal; ic++)
        for(int p=0; p<nNodesInCell; p++)
            sortNode[count++] = cell[ic][p];

    sort(sortNode.begin(), sortNode.end());
    sortNode.erase(unique(sortNode.begin(), sortNode.end()), sortNode.end());

    int nNodesGlobalTmp = nNodesGlobal;
    nNodesGlobal = 0;
    std::vector<int> sortNodeNew(sortNode.size(), 0);
    std::vector<int> convertNodeOldToNew(nNodesGlobalTmp, 0);
    for(int in=0; in<sortNode.size(); in++){
        sortNodeNew[in] = nNodesGlobal++;
        convertNodeOldToNew[sortNode[in]] = sortNodeNew[in];
    }

    for(int ic=0; ic<nCellsGlobal; ic++)
        for(int p=0; p<nNodesInCell; p++)
            cell[ic][p] = convertNodeOldToNew[cell[ic][p]];

    std::vector<std::vector<double>> nodeTmp = node;
    node.erase(node.begin(), node.end());

    for(int in=0; in<sortNode.size(); in++)
        node.push_back(nodeTmp[sortNode[in]]);

    std::map<int, std::vector<double>> vDirichletTmp = vDirichlet;
    std::map<int, double> pDirichletTmp = pDirichlet;

    vDirichlet.clear();
    pDirichlet.clear();

    for(int in=0; in<sortNode.size(); in++){
        if(vDirichletTmp.size() == 0) break;
        int key = sortNode[in];
        auto it = vDirichletTmp.find(key);
        if(it != vDirichletTmp.end()){
            int keyNew = convertNodeOldToNew[key];
            vDirichlet[keyNew] = vDirichletTmp[key];
        }
    }

    for(int in=0; in<sortNode.size(); in++){
        if(pDirichletTmp.size() == 0) break;
        int key = sortNode[in];
        auto it = pDirichletTmp.find(key);
        if(it != pDirichletTmp.end()){
            int keyNew = convertNodeOldToNew[key];
            pDirichlet[keyNew] = pDirichletTmp[key];
        }
    }
}
