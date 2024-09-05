/**
 * @file Node.h
 * @author K.Ueda
 * @date Jun, 2024
 */

#ifndef OCTREE_H
#define OCTREE_H

#include "Cell.h"
#include <iostream>
#include <vector>

class Octree
{
public:
  struct OctreeNode
  {
    double minX, maxX, minY, maxY, minZ, maxZ;
    std::vector<int> cellIndices;
    OctreeNode *children[8] = {nullptr};  
  };

  Octree(Cell &cell, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, int maxDepth)
  {
    std::vector<int> cellIdRoot;
    cellIdRoot.reserve(cell.nCellsGlobal); 
    for(int ic = 0; ic < cell.nCellsGlobal; ic++) {
      cellIdRoot.push_back(ic);
    }
    root = buildOctree(cell, cellIdRoot, minX, maxX, minY, maxY, minZ, maxZ, maxDepth);
  }

  std::vector<int> getCandidateCells(std::vector<double> &coord) const
  {
    return searchOctree(root, coord);
  }

private:
  OctreeNode *root;

  OctreeNode *buildOctree(Cell &cell, const std::vector<int> &cellIdParent, double minX, double maxX, double minY,
                          double maxY, double minZ, double maxZ, int depth)
  {
    OctreeNode *node = new OctreeNode();
    node->minX = minX;
    node->maxX = maxX;
    node->minY = minY;
    node->maxY = maxY;
    node->minZ = minZ;
    node->maxZ = maxZ;

    if(depth == 0 || cellIdParent.size() <= 1) {
      node->cellIndices = cellIdParent; 
      return node;
    }

    double midX = (minX + maxX) / 2.0;
    double midY = (minY + maxY) / 2.0;
    double midZ = (minZ + maxZ) / 2.0;

    auto isCellBoundaryIncludedInVoxel = [&](const int i, double minX, double maxX, double minY, double maxY,
                                             double minZ, double maxZ) {
      return !(maxX < cell(i).minX || minX > cell(i).maxX || maxY < cell(i).minY || minY > cell(i).maxY ||
               maxZ < cell(i).minZ || minZ > cell(i).maxZ);
    };

    std::vector<std::vector<int>> cellsIdChildren(8);

    for(int ic = 0; ic < 8; ic++) {
      double childMinX = (ic % 2 == 1) ? midX : minX;
      double childMaxX = (ic % 2 == 1) ? maxX : midX;
      double childMinY = ((ic / 2) % 2 == 1) ? midY : minY;
      double childMaxY = ((ic / 2) % 2 == 1) ? maxY : midY;
      double childMinZ = (ic / 4 == 1) ? midZ : minZ;
      double childMaxZ = (ic / 4 == 1) ? maxZ : midZ;

      for(int ip = 0; ip < cellIdParent.size(); ip++) {
        if(isCellBoundaryIncludedInVoxel(cellIdParent[ip], childMinX, childMaxX, childMinY, childMaxY, childMinZ,
                                         childMaxZ)) {
          cellsIdChildren[ic].push_back(cellIdParent[ip]);
        }
      }

      if(!cellsIdChildren[ic].empty()) {
        node->children[ic] = buildOctree(cell, cellsIdChildren[ic], childMinX, childMaxX, childMinY, childMaxY,
                                         childMinZ, childMaxZ, depth - 1);
      }
    }

    return node;
  }

  std::vector<int> searchOctree(OctreeNode *node, std::vector<double> &coord) const
  {
    if(!node) return {};

    if(coord[0] < node->minX || coord[0] > node->maxX || 
       coord[1] < node->minY || coord[1] > node->maxY ||
       coord[2] < node->minZ || coord[2] > node->maxZ) {
      return {};
    }

    if(isLeafNode(node)) {
      return node->cellIndices;
    }

    for(int i = 0; i < 8; i++) {
      auto result = searchOctree(node->children[i], coord);
      if(!result.empty()) {
        return result;
      }
    }

    return {};
  }

  bool isLeafNode(const OctreeNode *node) const
  {
    for(int i = 0; i < 8; i++) {
      if(node->children[i] != nullptr) {
        return false;
      }
    }
    return true;
  }
};

#endif
