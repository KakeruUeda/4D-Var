#include "Config.h"

void Config::readGridTypeParameter()
{
    std::string str, base_label, label;

    base_label = "/Grid";

    label = base_label + "/type";
    if(!tp.getInspectedValue(label, gridTypeString))
        throw std::runtime_error(label + " is not set");
    
    if(gridTypeString != "Structured" && gridTypeString != "Unstructured")
        throw std::runtime_error("Unknown GridType");  

    return;
}

void Config::readBasicParameter()
{
    std::string str, base_label, label;

    base_label = "/Base";
    label = base_label + "/dim";

    if(!tp.getInspectedValue(label, dim))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/outputDir";
    if(!tp.getInspectedValue(label, outputDir))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/numOfOMP";
    if(!tp.getInspectedValue(label, nOMP))
        throw std::runtime_error(label + " is not set");

    return;
}

void Config::readPysicalParameter()
{
    std::string str, base_label, label;

    base_label = "/PysicalParameter";
    label = base_label + "/rho";
    if (!tp.getInspectedValue(label, rho))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/mu";
    if (!tp.getInspectedValue(label, mu))
        throw std::runtime_error(label + " is not set");

    return;
}


void Config::readTimeParameter()
{
    /*
    string str,base_label,label;

    base_label = "/TimeParameter";
    label = base_label + "/dt";
    if (!tp.getInspectedValue(label, dt))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/timeMax";
    if (!tp.getInspectedValue(label, timeMax))
        throw std::runtime_error(label + " is not set");

    string ON_OFF;

    label = base_label + "/pulsatile_flow";
    if (!tp.getInspectedValue(label, ON_OFF))
        throw std::runtime_error(label + " is not set");
 
    if(ON_OFF == "ON")
        pulsatile_flow = ON;
    else if(ON_OFF == "OFF")
        pulsatile_flow = OFF;
    else 
        throw std::runtime_error("ON or OFF is not set");
  
    label = base_label + "/pulse_begin_itr";
    if (!tp.getInspectedValue(label, pulse_begin_itr))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/T";
    if (!tp.getInspectedValue(label, T))
        throw std::runtime_error(label + " is not set");
    */
    return;
}

void Config::readImageData()
{
    std::string str, base_label, label;
    std::string imageFile;

    base_label = "/Domain";
    label = base_label + "/image";

    if(!tp.getInspectedValue(label, imageFile))
        throw std::runtime_error(label + " is not set");
    
    std::ifstream ifsImage(imageFile);
    if(!ifsImage)
        throw std::runtime_error(imageFile + " open error");

    while(getline(ifsImage, str)){
        std::istringstream iss(str);
        double line;

        getline(iss, str, ' ');
        line = stod(str);

        phi.push_back(line);
    }
    ifsImage.close();
}

void Config::readGridParameter()
{
    std::string str, base_label, label;
    std::string nodeFile;
    label = base_label + "/node";

    if(!tp.getInspectedValue(label, nodeFile))
        throw std::runtime_error(label + " is not set");

    std::ifstream ifsNode(nodeFile);
    if(!ifsNode)
        throw std::runtime_error(nodeFile + " open error");

    while(getline(ifsNode, str)){
        std::istringstream iss(str);
        std::vector<double> nodeTmp;

        for(int d=0; d<dim; d++){
            getline(iss, str, ' ');
            nodeTmp.push_back(stod(str));
        }
        node.push_back(nodeTmp);
    }
    ifsNode.close();

    std::string cellFile;
    label = base_label + "/cell";

    if(!tp.getInspectedValue(label, cellFile))
        throw std::runtime_error(label + " is not set");

    std::ifstream ifsCell(cellFile);

    while(getline(ifsCell, str)){
        std::istringstream iss(str);
        std::vector<double> cellTmp;

        for(int d=0; d<nNodesInCell; d++){
            getline(iss, str, ' ');
            cellTmp.push_back(stod(str));
        }
        cell.push_back(cellTmp);
    }
    ifsCell.close();

}

void Config::readStructuredGridParameter()
{
    std::string str, base_label, label;
    std::string imageFile;
    int tmpInt[dim];
    double tmpDouble[dim];
  
    base_label = "/UnstructuredGrid";
    label = base_label + "/nx";
    if (!tp.getInspectedVector(label, tmpInt, dim))
        throw std::runtime_error(label + " is not set");

    nx = tmpInt[0];
    ny = tmpInt[1];
    nz = tmpInt[2];

    if(dim == 2){
        nxNodes = nx + 1; 
        nyNodes = ny + 1;
        nzNodes = 1;
        nxCells = nx;
        nyCells = ny;
        nzCells = 1;
    }
    if(dim == 3){
        nxNodes = nx + 1; 
        nyNodes = ny + 1;
        nzNodes = nz + 1;
        nxCells = nx;
        nyCells = ny;
        nzCells = nz;
    }

    label = base_label + "/lx";
    if(!tp.getInspectedVector(label, tmpDouble, dim))
        throw std::runtime_error(label + " is not set");

    lx = tmpDouble[0];
    ly = tmpDouble[1];
    if(dim == 2)      lz = 0.0;
    else if(dim == 3) lz = tmpDouble[2];

    dx = lx / (double)nx;
    dy = ly / (double)ny;
    dz = lz / (double)nz;

    nCellsGlobal = nxCells * nyCells * nzCells;
    nNodesGlobal = nxNodes * nyNodes * nzNodes;

    label = base_label + "/nNodesInCell";
    if (!tp.getInspectedValue(label, nNodesInCell))
        throw std::runtime_error(label + " is not set"); 

    if(nNodesInCell == 4 && dim != 2)
        throw std::runtime_error("nNodesInCell is not consistent with dim");
    
    if(nNodesInCell == 8 && dim != 3)
        throw std::runtime_error("nNodesInCell is not consistent with dim");

    return;
}

void Config::readStructuredBoundaryParameter()
{
    std::string str, base_label;
    std::string labelType, labelValue;
    std::string bdTypeTmp;
    int tmp = 0;

    base_label = "/Boundary";

    if(dim == 2){
        bdStr.push_back("bottom");
        labelType = base_label + "/bottom/type";
        labelValue = base_label + "/bottom/value";
        readBoundaryTypeAndValue(labelType, labelValue, tmp);

        bdStr.push_back("top");
        labelType = base_label + "/top/type";
        labelValue = base_label + "/top/value";
        readBoundaryTypeAndValue(labelType, labelValue, tmp);

        bdStr.push_back("left");
        labelType = base_label + "/left/type";
        labelValue = base_label + "/left/value";
        readBoundaryTypeAndValue(labelType, labelValue, tmp);

        bdStr.push_back("right");
        labelType = base_label + "/right/type";
        labelValue = base_label + "/right/value";
        readBoundaryTypeAndValue(labelType, labelValue, tmp);
    }else if(dim == 3){
        bdStr.push_back("bottom");
        labelType = base_label + "/bottom/type";
        labelValue = base_label + "/bottom/value";
        readBoundaryTypeAndValue(labelType, labelValue, tmp);

        bdStr.push_back("top");
        labelType = base_label + "/top/type";
        labelValue = base_label + "/top/value";
        readBoundaryTypeAndValue(labelType, labelValue, tmp);

        bdStr.push_back("left");
        labelType = base_label + "/left/type";
        labelValue = base_label + "/left/value";
        readBoundaryTypeAndValue(labelType, labelValue, tmp);

        bdStr.push_back("right");
        labelType = base_label + "/right/type";
        labelValue = base_label + "/right/value";
        readBoundaryTypeAndValue(labelType, labelValue, tmp);

        bdStr.push_back("front");
        labelType = base_label + "/front/type";
        labelValue = base_label + "/front/value";
        readBoundaryTypeAndValue(labelType, labelValue, tmp);

        bdStr.push_back("back");
        labelType = base_label + "/back/type";
        labelValue = base_label + "/back/value";
        readBoundaryTypeAndValue(labelType, labelValue, tmp);
    }else{
        throw std::runtime_error("Undefined dim");
    }

    return;
}

void Config::readBoundaryTypeAndValue(std::string labelType, std::string labelValue, int &tmp)
{
    std::string bdTypeTmp;
    if (!tp.getInspectedValue(labelType, bdTypeTmp))
        throw std::runtime_error(labelType + " is not set");

    bdType.push_back(bdTypeTmp);

    if(bdTypeTmp == "v"){
        double value[dim];
        if (!tp.getInspectedVector(labelValue, value, dim))
            throw std::runtime_error(labelValue + " is not set");  

        bdValue.emplace_back();
        for(int k=0; k<dim; k++)
            bdValue[tmp].push_back(value[k]);

    }else if(bdTypeTmp == "p"){
        double value;
        if(!tp.getInspectedValue(labelValue, value))
            throw std::runtime_error(labelValue + " is not set");

        bdValue.emplace_back();
        bdValue[tmp].push_back(value);

    }else if(bdTypeTmp == "free"){
        bdValue.emplace_back();
    }else if(bdTypeTmp == "file"){
    }
    else{
        throw std::runtime_error("label " + bdTypeTmp + " undefined");
    }

    tmp++;

    return;
}
