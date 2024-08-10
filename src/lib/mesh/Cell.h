/**
 * @file Cell.h
 * @author K.Ueda
 * @date Jun, 2024
 */

#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <cassert>
#include "Array.h"
#include "Config.h"
#include "VTKCellType.h"

struct CellInfo
{
public:
    VTKCellType cellType;
    int subId;
    double phi;

    std::vector<int> node, nodeNew;
    std::vector<int> dofsMap, dofsBCsMap;
    std::vector<int> dofStart;
    std::vector<int> dofStartPlane;
    std::vector<std::vector<double>> x;

    /// add ///
    std::vector<int> dofsMapWall, dofsBCsMapWall;

    inline void setArrayZero(int n);
};

class Cell
{
public:
    Cell() {}
    Cell(Config &conf) : nNodesInCell(conf.nNodesInCell),
                         nCellsGlobal(conf.nCellsGlobal), data(conf.nCellsGlobal) {}
    virtual ~Cell() {}

    inline CellInfo &operator()(int n)
    {
        return data[n];
    }

    inline int size()
    {
        return data.size();
    }

    inline void resize(int n)
    {
        data.resize(n);
    }

    int nCellsGlobal;
    int nNodesInCell;
    int nCellsStrGlobal;

    void initialize(Config &conf);
    void initializeAdjoint(Config &conf);

private:
    std::vector<CellInfo> data;
};

#endif