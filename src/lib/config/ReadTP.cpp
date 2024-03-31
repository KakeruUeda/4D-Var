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
    int tmp_dim, tmp_nOMP;

    base_label = "/Base";
    label = base_label + "/dim";

    if(!tp.getInspectedValue(label, tmp_dim))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/outputDir";
    if(!tp.getInspectedValue(label, outputDir))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/numOfOMP";
    if(!tp.getInspectedValue(label, tmp_nOMP))
        throw std::runtime_error(label + " is not set");

    dim  = static_cast<size_t>(tmp_dim);
    nOMP = static_cast<size_t>(tmp_nOMP);

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

void Config::readBoundaryMethodParameter()
{
    /*
    string str,base_label,label;
    string method;

    base_label = "/BoundaryMethod";
    label = base_label + "/boundary";
    if (!tp.getInspectedValue(label, method))
        throw std::runtime_error(label + " is not set");

    if(method == "XFEM")
        bd = BOUNDARY::XFEM;
    else if(method == "Darcy")
        bd = BOUNDARY::DARCY;
    else
        throw std::runtime_error("Boundary method is not defined");

    return;
    */
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

void Config::readGridParameter()
{
    std::string str, base_label, label;
    int tmpInt[dim];
    double tmpDouble[dim];
    int nNodesInCellTmp;
  
    base_label = "/Grid";
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
    if (!tp.getInspectedVector(label, tmpDouble, dim))
        throw std::runtime_error(label + " is not set");

    lx = tmpDouble[0];
    ly = tmpDouble[1];
    if(dim == 2) lz = 0.0;
    else if(dim == 3) lz = tmpDouble[2];

    dx = lx / (double)nx;
    dy = ly / (double)ny;
    dz = lz / (double)nz;

    nCellsGlobal = nx * ny * nz;
    if(dim == 2) nNodesGlobal = (nx+1) * (ny+1);
    else if(dim == 3) nNodesGlobal = (nx+1) * (ny+1) * (nz+1);

    label = base_label + "/nNodesInCell";
    if (!tp.getInspectedValue(label, nNodesInCellTmp))
        throw std::runtime_error(label + " is not set"); 

    nNodesInCell  = static_cast<size_t>(nNodesInCellTmp);

    return;
}

void Config::readImageParameter()
{
    std::string str, base_label, label;
    std::string imageFile;

    phi.resize(nCellsGlobal);
    for(int i=0;i<nCellsGlobal;i++)
        phi.at(i) = 1.0;

    label = "/Domain/image";
    if(!tp.getInspectedValue(label, imageFile))
        throw std::runtime_error(label + " is not set");
    
    std::ifstream ifsImage(imageFile);
    if(!ifsImage)
        throw std::runtime_error(imageFile + "open error");
    
    while(getline(ifsImage, str))
    {
        std::istringstream iss(str);
        double line;

        for(int i=0; i<dim+1; i++)
        {
            if(i < dim) continue;
            getline(iss, str, ' ');
            line = stod(str);
        }
        phi.push_back(line);
    }
    ifsImage.close();

    return;
}

void Config::readBoundaryParameter()
{
    std::string str, base_label;
    std::string labelType, labelValue;
    std::string bdTypeTmp;
    size_t tmp = 0;

    base_label = "/Boundary";

    if(dim == 2)
    {
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
    }
    else if(dim == 3)
    {
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
    }
    else
    {
        throw std::runtime_error("Undefined dim");
    }

    return;
}

void Config::readBoundaryTypeAndValue(std::string labelType, std::string labelValue, size_t &tmp)
{
    std::string bdTypeTmp;
    if (!tp.getInspectedValue(labelType, bdTypeTmp))
        throw std::runtime_error(labelType + " is not set");

    bdType.push_back(bdTypeTmp);

     if(bdTypeTmp == "v")
    {
        double value[dim];
        if (!tp.getInspectedVector(labelValue, value, dim))
            throw std::runtime_error(labelValue + " is not set");  
        bdValue.emplace_back();
        for(int k=0; k<dim; k++)
            bdValue[tmp].push_back(value[k]);
    }
    else if(bdTypeTmp == "p")
    {
        double value;
        if(!tp.getInspectedValue(labelValue, value))
            throw std::runtime_error(labelValue + " is not set");
        bdValue.emplace_back();
        bdValue[tmp].push_back(value);
    }
    else if(bdTypeTmp == "free")
    {
        bdValue.emplace_back();
    }
    else if(bdTypeTmp == "file")
    {
    }
    else
    {
        throw std::runtime_error("label " + bdTypeTmp + " undefined");
    }

    tmp++;

    return;
}
