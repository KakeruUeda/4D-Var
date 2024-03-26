#include "Config.h"

void Config::readGridType()
{
    std::string str, base_label, label;

    base_label = "/Grid";

    std::string gridTypeString;
    label = base_label + "/type";
    if(!tp.getInspectedValue(label, gridTypeString))
        throw std::runtime_error(label + " is not set");
    
    if(gridTypeString == "Structured")
        gridType = GridType::STRUCTURED;
    else if(gridTypeString == "Unstructured")
        gridType = GridType::UNSTRUCTURED;
    else 
        throw std::runtime_error("Unknown grid type");

    return;
}

void Config::readBase()
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

void Config::readPysicalParam()
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

void Config::readBoundaryMethod()
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

void Config::readTimeParam()
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

void Config::readGrid()
{
    std::string str, base_label, label;
    int tmpInt[dim];
    double tmpDouble[dim];
  
    base_label = "/Grid";
    label = base_label + "/nx";
    if (!tp.getInspectedVector(label, tmpInt, dim))
        throw std::runtime_error(label + " is not set");

    nx.resize(dim);
    for(int i=0; i<dim; i++)
        nx(i) = tmpInt[i];

    label = base_label + "/lx";
    if (!tp.getInspectedVector(label, tmpDouble, dim))
        throw std::runtime_error(label + " is not set");

    lx.resize(dim);
    for(int i=0; i<dim; i++)
        lx(i) = tmpDouble[i];

    return;
}