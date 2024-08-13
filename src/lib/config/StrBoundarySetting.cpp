/**
 * @file SetStrBoundary.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

void Config::setBoundaryVelocityValue(std::string face, double value[3])
{
  std::vector<double> vecValue(value, value + 3);

  auto updateBoundary = [&](auto condition)
  {
    for(int k=0; k<nz+1; k++){
      for(int j=0; j<ny+1; j++){
        for(int i=0; i<nx+1; i++){
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
          if(condition(i, j, k)){
            vDirichlet[in] = vecValue;
          }
        }
      }
    }
  };

  if (face == "top")
    updateBoundary([&](int, int j, int){ return j == ny; });
  else if (face == "bottom")
    updateBoundary([&](int, int j, int){ return j == 0; });
  else if (face == "left")
    updateBoundary([&](int i, int, int){ return i == 0; });
  else if (face == "right")
    updateBoundary([&](int i, int, int){ return i == nx; });
  else if (face == "front")
    updateBoundary([&](int, int, int k){ return k == 0; });
  else if (face == "back")
    updateBoundary([&](int, int, int k){ return k == nz; });
}

void Config::setBoundaryPressureValue(std::string face, const double value)
{
  auto updateBoundary = [&](auto condition)
  {
    for(int k=0; k<nz+1; k++){
      for(int j=0; j<ny+1; j++){
        for(int i=0; i<nx+1; i++){
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
          if(condition(i, j, k)){
            pDirichlet[in] = value;
          }
        }
      }
    }
  };

  if(face == "top")
    updateBoundary([&](int, int j, int){ return j == ny; });
  else if(face == "bottom")
    updateBoundary([&](int, int j, int){ return j == 0; });
  else if(face == "left")
    updateBoundary([&](int i, int, int){ return i == 0; });
  else if(face == "right")
    updateBoundary([&](int i, int, int){ return i == nx; });
  else if(face == "front")
    updateBoundary([&](int, int, int k){ return k == 0; });
  else if(face == "back")
    updateBoundary([&](int, int, int k){ return k == nz; });
}

void Config::setBoundaryPoiseuilleValue(std::string face)
{
  auto addVelocity = [&](int i, int j, int k, double value)
  {
    int n = i + j * (nx + 1) + k * (nx + 1) * (ny + 1);

    if(face == "left" || face == "right"){
      vDirichlet[n] = {value, 0.0, 0.0};
    }else if(face == "top" || face == "bottom"){
      vDirichlet[n] = {0.0, value, 0.0};
    }else if(face == "front" || face == "back"){
      vDirichlet[n] = {0.0, 0.0, value};
    }
  };

  auto calculatePoiseuilleVelocity = [&](int i, int j, int k) -> std::tuple<double, bool>
  {
    double r = 0.0, value = 0.0;
    bool onBoundary = false;

    if((face == "left" && i == 0) || (face == "right" && i == nx)){
      r = sqrt(pow(center[2] - k * dz, 2) + pow(center[1] - j * dy, 2));
      onBoundary = true;
    }
    else if((face == "top" && j == ny) || (face == "bottom" && j == 0)){
      r = sqrt(pow(center[0] - i * dx, 2) + pow(center[2] - k * dz, 2));
      onBoundary = true;
    }
    else if((face == "front" && k == 0) || (face == "back" && k == nz)){
      r = sqrt(pow(center[0] - i * dx, 2) + pow(center[1] - j * dy, 2));
      onBoundary = true;
    }

    if (onBoundary){
      if (r < R){
        value = (2.0 * Q / (PI * R * R)) * (1 - (r * r) / (R * R));
      }else{
        value = 0.0;
      }
    }

    return {value, onBoundary};
  };

  for(int k=0; k<nz+1; k++){
    for(int j=0; j<ny+1; j++){
      for(int i=0; i<nx+1; i++){
        auto [value, onBoundary] = calculatePoiseuilleVelocity(i, j, k);
        if (onBoundary){
          addVelocity(i, j, k, value);
        }
      }
    }
  }
}

void Config::setControlBoundary()
{
  auto collectBoundaryNodes = [&](int i, int j, int k)
  {
    int in = i + j * (nx+1) + k * (nx+1) * (ny+1);
    mapCB.push_back(in);
  };

  auto collectBoundaryCells = [&](int i, int j, int k)
  {
    int ic = i + j * nx + k * nx * ny;
    mapCBCell.push_back(ic);
  };

  auto collectBoundaryCellNodes = [&](int i, int j, int k, const std::vector<int> &nodeIndices)
  {
    int ic = i + j * nx + k * nx * ny;
    std::vector<int> vecTmp;
    for (int p : nodeIndices){
      vecTmp.push_back(cell[ic][p]);
    }
    mapCBInCell.push_back(vecTmp);
  };

  auto collectBoundary = [&](auto conditionNodes, auto conditionCells, const std::vector<int> &nodeIndices)
  {
    for(int k=0; k<nz+1; k++){
      for(int j=0; j<ny+1; j++){
        for(int i=0; i<nx+1; i++){
          if(conditionNodes(i, j, k)){
            collectBoundaryNodes(i, j, k);
          }
        }
      }
    }

    for(int k=0; k<nz; k++){
      for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
          if(conditionCells(i, j, k)){
            collectBoundaryCells(i, j, k);
            collectBoundaryCellNodes(i, j, k, nodeIndices);
          }
        }
      }
    }
  };

  switch (inletCB)
  {
  case ControlBoundary::left:
    collectBoundary(
        [&](int i, int, int){ return i == 0; },
        [&](int i, int, int){ return i == 0; },
        {0, 3, 7, 4});
    break;
  case ControlBoundary::right:
    collectBoundary(
        [&](int i, int, int){ return i == nx; },
        [&](int i, int, int){ return i == nx-1; },
        {1, 5, 6, 2});
    break;
  case ControlBoundary::top:
    collectBoundary(
        [&](int, int j, int){ return j == ny; },
        [&](int, int j, int){ return j == ny-1; },
        {2, 6, 7, 3});
    break;
  case ControlBoundary::bottom:
    collectBoundary(
        [&](int, int j, int){ return j == 0; },
        [&](int, int j, int){ return j == 0; },
        {0, 4, 5, 1});
    break;
  case ControlBoundary::front:
    collectBoundary(
        [&](int, int, int k){ return k == 0; },
        [&](int, int, int k){ return k == 0; },
        {0, 1, 2, 3});
    break;
  case ControlBoundary::back:
    collectBoundary(
        [&](int, int, int k){ return k == nz; },
        [&](int, int, int k){ return k == nz-1; },
        {4, 7, 6, 5});
    break;
  default:
    break;
  }
}
