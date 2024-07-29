/**
 * @file ReadTP.cpp
 * @author K.Ueda
 * @date May, 2024
*/

#include "Config.h"

/*****************************
 * @brief Read text parameter.
 */
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
}


/*****************************
 * @brief Read text parameter.
 */
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
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readNRParameter()
{
    std::string str, base_label, label;

    base_label = "/NRParameter";
    label = base_label + "/NRtolerance";
    if (!tp.getInspectedValue(label, NRtolerance))
        throw std::runtime_error(label + " is not set");
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readTimeParameter()
{
    std::string str, base_label, label;

    base_label = "/TimeParameter";
    label = base_label + "/dt";
    if (!tp.getInspectedValue(label, dt))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/timeMax";
    if (!tp.getInspectedValue(label, timeMax))
        throw std::runtime_error(label + " is not set");

    std::string ON_OFF;

    label = base_label + "/pulsatileFlow";
    if (!tp.getInspectedValue(label, ON_OFF))
        throw std::runtime_error(label + " is not set");
 
    if(ON_OFF == "ON")
        pulsatileFlow = ON;
    else if(ON_OFF == "OFF")
        pulsatileFlow = OFF;
    else 
        throw std::runtime_error("ON or OFF is not set");
  
    label = base_label + "/pulseBeginItr";
    if (!tp.getInspectedValue(label, pulseBeginItr))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/T";
    if (!tp.getInspectedValue(label, T))
        throw std::runtime_error(label + " is not set");
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readDarcyParameter()
{
    std::string str, base_label, label;
    base_label = "/DarcyParameter";

    label = base_label + "/alpha";
    if(!tp.getInspectedValue(label, alpha))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/resistance";
    if(!tp.getInspectedValue(label, resistance))
        throw std::runtime_error(label + " is not set");
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readGridParameter()
{
    std::string str, base_label, label;
    std::string nodeFile;

    base_label = "/Grid";

    label = base_label + "/type";
    std::string gridTypeString;

    if(!tp.getInspectedValue(label, gridTypeString))
        throw std::runtime_error(label + " is not set");
    
    if(gridTypeString != "Structured" && gridTypeString != "Unstructured")
        throw std::runtime_error("Unknown GridType");  

    if(gridTypeString == "Structured")
        gridType = GridType::STRUCTURED;
    else if(gridTypeString == "Unstructured")
        gridType = GridType::UNSTRUCTURED;

    if(gridType == GridType::STRUCTURED){
        int tmpInt[dim];
        double tmpDouble[dim];

        label = base_label + "/nx";
        if(!tp.getInspectedVector(label, tmpInt, dim))
            throw std::runtime_error(label + " is not set.");

        nx = tmpInt[0];
        ny = tmpInt[1];
        nz = tmpInt[2];

        label = base_label + "/lx";
        if(!tp.getInspectedVector(label, tmpDouble, dim))
            throw std::runtime_error(label + " is not set");

        lx = tmpDouble[0];
        ly = tmpDouble[1];
        lz = tmpDouble[2];

        dx = lx / (double)nx;
        dy = ly / (double)ny;
        dz = lz / (double)nz;
    }

    std::string ON_OFF;

    label = base_label + "/extractFluid";
    if (!tp.getInspectedValue(label, ON_OFF))
        throw std::runtime_error(label + " is not set");
 
    if(ON_OFF == "ON")
        extractFluid = ON;
    else if(ON_OFF == "OFF")
        extractFluid = OFF;
    else 
        throw std::runtime_error("ON or OFF is not set");

    label = base_label + "/nNodesInCell";
    if (!tp.getInspectedValue(label, nNodesInCell))
        throw std::runtime_error(label + " is not set"); 

    if(nNodesInCell == 4 && dim != 2)
        throw std::runtime_error("nNodesInCell is not consistent with dim");
    
    if(nNodesInCell == 8 && dim != 3)
        throw std::runtime_error("nNodesInCell is not consistent with dim");
    
    label = base_label + "/node";

    if(!tp.getInspectedValue(label, nodeFile))
        throw std::runtime_error(label + " is not set");

    std::ifstream ifsNode(nodeFile);
    if(!ifsNode.is_open()){
        throw std::runtime_error("Failed to open file: " + nodeFile);
    }
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

    nNodesGlobal = node.size();

    std::string cellFile;
    label = base_label + "/cell";

    if(!tp.getInspectedValue(label, cellFile))
        throw std::runtime_error(label + " is not set");

    std::ifstream ifsCell(cellFile);
    if(!ifsCell.is_open()){
        throw std::runtime_error("Failed to open file: " + cellFile);
    }
    while(getline(ifsCell, str)){
        std::istringstream iss(str);
        std::vector<int> cellTmp;

        for(int d=0; d<nNodesInCell; d++){
            getline(iss, str, ' ');
            cellTmp.push_back(stod(str));
        }
        cell.push_back(cellTmp);
    }
    ifsCell.close();

    nCellsGlobal = cell.size();

    std::string imageFile;
    label = base_label + "/image";

    if(!tp.getInspectedValue(label, imageFile))
        throw std::runtime_error(label + " is not set");
    
    std::ifstream ifsImage(imageFile);
    if(!ifsImage.is_open()){
        throw std::runtime_error("Failed to open file: " + imageFile);
    }
    while(getline(ifsImage, str)){
        std::istringstream iss(str);
        for(int d=0; d<4; d++){
            getline(iss, str, ' ');
            if(d == 3) phi.push_back(stod(str));
        }
    }
    ifsImage.close();
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readBoundaryParameter()
{
    std::string str, base_label, label;
    std::string nodeFile;

    base_label = "/Boundary";

    std::string velFile;
    label = base_label + "/velocityDirichlet";

    if(!tp.getInspectedValue(label, velFile))
        throw std::runtime_error(label + " is not set");

    std::ifstream ifsVel(velFile);
    if(!ifsVel.is_open()){
        throw std::runtime_error("Failed to open file: " + velFile);
    }
     while(getline(ifsVel, str)){
        int index;
        std::istringstream iss(str);
        std::vector<double> vecTmp;

        for(int d=0; d<dim+1; d++){
            getline(iss, str, ' ');
            if(d == 0) index = stoi(str);
            else vecTmp.push_back(stod(str));
        }
        vDirichlet[index] = vecTmp;
    }
    ifsVel.close();

    std::string preFile;
    label = base_label + "/pressureDirichlet";

    if(!tp.getInspectedValue(label, preFile))
        throw std::runtime_error(label + " is not set");

    std::ifstream ifsPre(preFile);
    if(!ifsPre.is_open()){
        throw std::runtime_error("Failed to open file: " + preFile);
    }
     while(getline(ifsPre, str)){
        int index;
        std::istringstream iss(str);
        std::vector<double> preTmp;

        for(int d=0; d<1+1; d++){
            getline(iss, str, ' ');
            if(d == 0) index = stoi(str);
            else pDirichlet[index] = stoi(str);
        }
    }
    ifsPre.close();
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readControlBoundaryParameter()
{
    std::string str, base_label, label;
    base_label = "/Boundary";
    label = base_label + "/control";

    if(!tp.getInspectedValue(label, str))
        throw std::runtime_error(label + " is not set");
    
    if(str == "left")
        controlBoundary = ControlBoundary::left;
    if(str == "right")
        controlBoundary = ControlBoundary::right;
    if(str == "top")
        controlBoundary = ControlBoundary::top;
    if(str == "bottom")
        controlBoundary = ControlBoundary::bottom;
    if(str == "front")
        controlBoundary = ControlBoundary::front;
    if(str == "back")
        controlBoundary = ControlBoundary::back;
}


/*****************************
 * @brief Read text parameter.
 */
void Config::readInverseParameter()
{
    std::string str, base_label, label;
    base_label = "/Inverse"; 
    
    label = base_label + "/controlBoundary";

    if(!tp.getInspectedValue(label, str))
        throw std::runtime_error(label + " is not set");

    planeDir.resize(2, 0);
    
    if(str == "left"){
        controlBoundary = ControlBoundary::left;
        planeDir[0] = 1, planeDir[1] = 2;
    }else if(str == "right"){
        controlBoundary = ControlBoundary::right;
        planeDir[0] = 1, planeDir[1] = 2;
    }else if(str == "top"){
        controlBoundary = ControlBoundary::top;
        planeDir[0] = 0, planeDir[1] = 2;
    }else if(str == "bottom"){
        controlBoundary = ControlBoundary::bottom;
        planeDir[0] = 0, planeDir[1] = 2;
    }else if(str == "front"){
        controlBoundary = ControlBoundary::front;
        planeDir[0] = 0, planeDir[1] = 1;
    }else if(str == "back"){
        controlBoundary = ControlBoundary::back;
        planeDir[0] = 0, planeDir[1] = 1;
    }

    label = base_label + "/voxelVelocity";

    if(!tp.getInspectedValue(label, str))
        throw std::runtime_error(label + " is not set");

    if(str == "average"){
        vvox = VoxelVelocity::AVERAGE;
    }else if(str == "interpolation"){
        vvox = VoxelVelocity::INTERPOLATION;
    }else{
        throw std::runtime_error("undefined voxelVelocity");
    }

    label = base_label + "/aCF";
    if(!tp.getInspectedValue(label, aCF))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/bCF1";
    if(!tp.getInspectedValue(label, bCF1))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/bCF2";
    if(!tp.getInspectedValue(label, bCF2))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/gCF";
    if(!tp.getInspectedValue(label, gCF))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/outputItr";
    if(!tp.getInspectedValue(label, outputItr))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/loopMax";
    if(!tp.getInspectedValue(label, loopMax))
        throw std::runtime_error(label + " is not set");
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readDataParameter()
{
    std::string str, base_label, label;
    base_label = "/Data"; 

    int tmpInt[dim];
    double tmpDouble[dim];

    label = base_label + "/nControlNodesInCell";
    if(!tp.getInspectedValue(label, nControlNodesInCell))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/nSnapShot";
    if(!tp.getInspectedValue(label, nSnapShot))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/snapInterval";
    if(!tp.getInspectedValue(label, snapInterval))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/snapTimeBeginItr";
    if(!tp.getInspectedValue(label, snapTimeBeginItr))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/nNodesInDataCell";
    if(!tp.getInspectedValue(label, nNodesInCellData))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/nxData";
    if (!tp.getInspectedVector(label, tmpInt, dim))
        throw std::runtime_error(label + " is not set");

    nxData = tmpInt[0];
    nyData = tmpInt[1];
    nzData = tmpInt[2];

    label = base_label + "/lxData";
    if(!tp.getInspectedVector(label, tmpDouble, dim))
        throw std::runtime_error(label + " is not set");

    lxData = tmpDouble[0];
    lyData = tmpDouble[1];
    lzData = tmpDouble[2];

    dxData = lxData / (double)nxData;
    dyData = lyData / (double)nyData;
    dzData = lzData / (double)nzData;

    nCellsDataGlobal = nxData * nxData * nzData;  

    std::string controlBoundaryFile;
    label = base_label + "/controlBoundary";

    if(!tp.getInspectedValue(label, controlBoundaryFile))
        throw std::runtime_error(label + " is not set");

    std::ifstream ifsControlBoundary(controlBoundaryFile);
    if(!ifsControlBoundary.is_open()){
        throw std::runtime_error("Failed to open file: " + controlBoundaryFile);
    }
    while(getline(ifsControlBoundary, str)){
        std::istringstream iss(str);
        getline(iss, str, ' ');
        controlBoundaryMap.push_back(stoi(str));
    }
    ifsControlBoundary.close();

    std::string controlCellMapFile;
    label = base_label + "/controlCellMap";

    if(!tp.getInspectedValue(label, controlCellMapFile))
        throw std::runtime_error(label + " is not set");

    std::ifstream ifscontrolCellMap(controlCellMapFile);
    if(!ifscontrolCellMap.is_open()){
        throw std::runtime_error("Failed to open file: " + controlCellMapFile);
    }
    while(getline(ifscontrolCellMap, str)){
        std::istringstream iss(str);
        getline(iss, str, ' ');
        controlCellMap.push_back(stoi(str));
    }
    ifscontrolCellMap.close();

    std::string controlNodeInCellFile;
    label = base_label + "/controlNodeInCell";

    if(!tp.getInspectedValue(label, controlNodeInCellFile))
        throw std::runtime_error(label + " is not set");

    std::ifstream ifsControlNodeInCell(controlNodeInCellFile);
    if(!ifsControlNodeInCell.is_open()){
        throw std::runtime_error("Failed to open file: " + controlNodeInCellFile);
    }
    while(getline(ifsControlNodeInCell, str)){
        std::istringstream iss(str);
        std::vector<int> ctrTmp;
        for(int p=0; p<nControlNodesInCell; p++){
            getline(iss, str, ' ');
            ctrTmp.push_back(stod(str));
        }
        controlNodeInCell.push_back(ctrTmp);
    }
    ifsControlNodeInCell.close();
    
    label = base_label + "/inputDir";
    if(!tp.getInspectedValue(label, inputDir))
        throw std::runtime_error(label + " is not set");
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readStructuredGridParameter()
{
    std::string str, base_label, label;
    std::string imageFile;
    int tmpInt[dim];
    double tmpDouble[dim];
  
    base_label = "/StructuredGrid";
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
    if(dim == 2)  lz = 0.0;
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

/*****************************
 * @brief Read text parameter.
 */
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

/*****************************
 * @brief Read text parameter.
 */
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

/*****************************
 * @brief Read text parameter.
 */
void Config::readPostInverseBasicParameter()
{
    std::string str, base_label, label;
    base_label = "/PostInverseBasic";

    label = base_label + "/gridType";
    std::string gridTypeString;

    if(!tp.getInspectedValue(label, gridTypeString))
        throw std::runtime_error(label + " is not set");
    
    if(gridTypeString != "Structured" && gridTypeString != "Unstructured")
        throw std::runtime_error("Unknown GridType");  

    if(gridTypeString == "Structured")
        gridType = GridType::STRUCTURED;
    else if(gridTypeString == "Unstructured")
        gridType = GridType::UNSTRUCTURED;

    if(gridType == GridType::STRUCTURED){
        int tmpInt[dim];
        double tmpDouble[dim];

        label = base_label + "/nx";
        if(!tp.getInspectedVector(label, tmpInt, dim))
            throw std::runtime_error(label + " is not set.");

        nx = tmpInt[0];
        ny = tmpInt[1];
        nz = tmpInt[2];

        label = base_label + "/lx";
        if(!tp.getInspectedVector(label, tmpDouble, dim))
            throw std::runtime_error(label + " is not set");

        lx = tmpDouble[0];
        ly = tmpDouble[1];
        lz = tmpDouble[2];

        dx = lx / (double)nx;
        dy = ly / (double)ny;
        dz = lz / (double)nz;
    }

    label = base_label + "/nNodesInCell";
    if(!tp.getInspectedValue(label, nNodesInCell))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/nRef";
    if(!tp.getInspectedValue(label, nRef))
        throw std::runtime_error(label + " is not set");
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readPostInverseVelocityParameter()
{
    std::string str, base_label, label;
    base_label = "/PostInverseVelocity";

    int num = 0;
    std::string velRefFile;
    velRef.resize(nRef);

    while(1){
        if(num >= nRef) break;
        label = base_label + "/velRef" + std::to_string(num);

        if(!tp.getInspectedValue(label, velRefFile))
            throw std::runtime_error(label + " is not set");
        
        std::ifstream ifsVelRef(velRefFile);
        if(!ifsVelRef.is_open()){
            throw std::runtime_error("Failed to open file: " + velRefFile);
        }
        std::vector<std::vector<double>> velRefTmp;

        while(getline(ifsVelRef, str)){
            std::istringstream iss(str);
            std::vector<double> tmp;
            for(int d=0; d<dim; d++){
                getline(iss, str, ' ');
                tmp.push_back(stod(str));
            }
            velRefTmp.push_back(tmp);
        } 
        ifsVelRef.close();
        velRef[num] = velRefTmp;
        num++;
    }

    num = 0;
    std::string velOptFile;
    velOpt.resize(nRef);

    while(1){
        if(num >= nRef) break;
        label = base_label + "/velOpt" + std::to_string(num);

        if(!tp.getInspectedValue(label, velOptFile))
            throw std::runtime_error(label + " is not set");
        
        std::ifstream ifsVelOpt(velOptFile);
        if(!ifsVelOpt.is_open()){
            throw std::runtime_error("Failed to open file: " + velOptFile);
        }
        std::vector<std::vector<double>> velOptTmp;

        while(getline(ifsVelOpt, str)){
            std::istringstream iss(str);
            std::vector<double> tmp;
            for(int d=0; d<dim; d++){
                getline(iss, str, ' ');
                tmp.push_back(stod(str));
            }
            velOptTmp.push_back(tmp);
        } 
        ifsVelOpt.close();
        velOpt[num] = velOptTmp;
        num++;
    }
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readPostInverseFlowRateParameter()
{
    std::string str, base_label, label;
    base_label = "/PostInverseFrowRate";
    
    label = base_label + "/crossSection";
    if(!tp.getInspectedValue(label, str))
        throw std::runtime_error(label + " is not set");
    
    label = base_label + "/flowRateVelDir";
    if(!tp.getInspectedValue(label, flowRateVelDir))
        throw std::runtime_error(label + " is not set");

    if(str == "xy"){
        crossSection = CrossSection::XY;
    }else if(str == "yz"){
        crossSection = CrossSection::YZ;
    }else if(str == "zx"){
        crossSection = CrossSection::ZX;
    }

    label = base_label + "/crossPoint";
    if(!tp.getInspectedValue(label, crossPoint))
        throw std::runtime_error(label + " is not set");
}

/*****************************
 * @brief Read text parameter.
 */
void Config::readVoxelCreationParameter()
{
    std::string str, base_label, label;

    int tmpInt[dim];
    double tmpDouble[dim];

    base_label = "/VoxelCreation";

    label = base_label + "/inputDir";
    if(!tp.getInspectedValue(label, inputDir))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/stepMax";
    if(!tp.getInspectedValue(label, stepMax))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/nSnapShot";
    if(!tp.getInspectedValue(label, nSnapShot))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/snapInterval";
    if(!tp.getInspectedValue(label, snapInterval))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/snapTimeBeginItr";
    if(!tp.getInspectedValue(label, snapTimeBeginItr))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/nNodesInDataCell";
    if(!tp.getInspectedValue(label, nNodesInCellData))
        throw std::runtime_error(label + " is not set");

    label = base_label + "/origin";
    if (!tp.getInspectedVector(label, tmpDouble, dim))
        throw std::runtime_error(label + " is not set");

    xOrigin = tmpDouble[0];
    yOrigin = tmpDouble[1];
    zOrigin = tmpDouble[2];

    label = base_label + "/nxData";
    if (!tp.getInspectedVector(label, tmpInt, dim))
        throw std::runtime_error(label + " is not set");

    nxData = tmpInt[0];
    nyData = tmpInt[1];
    nzData = tmpInt[2];

    label = base_label + "/lxData";
    if(!tp.getInspectedVector(label, tmpDouble, dim))
        throw std::runtime_error(label + " is not set");

    lxData = tmpDouble[0];
    lyData = tmpDouble[1];
    lzData = tmpDouble[2];

    dxData = lxData / (double)nxData;
    dyData = lyData / (double)nyData;
    dzData = lzData / (double)nzData;

    nCellsDataGlobal = nxData * nxData * nzData;

    label = base_label + "/nxOpt";
    if(!tp.getInspectedVector(label, tmpInt, dim))
        throw std::runtime_error(label + " is not set");

    nxOpt = tmpInt[0];
    nyOpt = tmpInt[1];
    nzOpt = tmpInt[2];

    label = base_label + "/lxOpt";
    if(!tp.getInspectedVector(label, tmpDouble, dim))
        throw std::runtime_error(label + " is not set");

    lxOpt = tmpDouble[0];
    lyOpt = tmpDouble[1];
    lzOpt = tmpDouble[2];

    dxOpt = lxOpt / (double)nxOpt;
    dyOpt = lyOpt / (double)nyOpt;
    dzOpt = lzOpt / (double)nzOpt;

    nCellsOptGlobal = nxOpt * nxOpt * nzOpt;
    nNodesOptGlobal = (nxOpt+1) * (nyOpt+1) * (nzOpt+1);
}

