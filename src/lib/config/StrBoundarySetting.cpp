/**
 * @file SetStrBoudnary.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

void Config::setBoundaryVelocityValue(std::string face, double value[3])
{
  std::vector<double> vecValue(3);
  for (int i = 0; i < 3; i++)
  {
    vecValue[i] = value[i];
  }

  if (face == "top")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (j == ny)
          {
            vDirichlet[in] = vecValue;
          }
        }
      }
    }
  }
  if (face == "bottom")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (j == 0)
          {
            vDirichlet[in] = vecValue;
          }
        }
      }
    }
  }
  if (face == "left")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (i == 0)
          {
            vDirichlet[in] = vecValue;
          }
        }
      }
    }
  }
  if (face == "right")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (i == nx)
          {
            vDirichlet[in] = vecValue;
          }
        }
      }
    }
  }
  if (face == "front")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (k == 0)
          {
            vDirichlet[in] = vecValue;
          }
        }
      }
    }
  }
  if (face == "back")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (k == ny)
          {
            vDirichlet[in] = vecValue;
          }
        }
      }
    }
  }
}

void Config::setBoundaryPressureValue(std::string face, const double value)
{
  if (face == "top")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (j == ny)
          {
            pDirichlet[in] = value;
          }
        }
      }
    }
  }
  if (face == "bottom")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (j == 0)
          {
            pDirichlet[in] = value;
          }
        }
      }
    }
  }
  if (face == "left")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (i == 0)
          {
            pDirichlet[in] = value;
          }
        }
      }
    }
  }
  if (face == "right")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (i == nx)
          {
            pDirichlet[in] = value;
          }
        }
      }
    }
  }
  if (face == "front")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (k == 0)
          {
            pDirichlet[in] = value;
          }
        }
      }
    }
  }
  if (face == "back")
  {
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        for (int i = 0; i < nx + 1; i++)
        {
          int in = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
          if (k == ny)
          {
            pDirichlet[in] = value;
          }
        }
      }
    }
  }
}

void Config::setBoundaryPoiseuilleValue(std::string face)
{
  auto addVelocity = [&](int i, int j, int k, double value)
  {
    int n = i + j * (nx + 1) + k * (nx + 1) * (ny + 1);

    if (face == "left" || face == "right")
    {
      vDirichlet[n] = {value, 0e0, 0e0};
    }
    else if (face == "top" || face == "bottom")
    {
      vDirichlet[n] = {0e0, value, 0e0};
    }
    else if (face == "front" || face == "back")
    {
      vDirichlet[n] = {0e0, 0e0, value};
    }
  };

  for (int k = 0; k < nz + 1; k++)
  {
    for (int j = 0; j < ny + 1; j++)
    {
      for (int i = 0; i < nx + 1; i++)
      {
        double r, value = 0.0;
        bool onBoundary = false;

        if ((face == "left" && i == 0) || (face == "right" && i == nx))
        {
          r = sqrt(pow(center[2] - k * dz, 2) + pow(center[1] - j * dy, 2));
          onBoundary = true;
        }
        else if ((face == "top" && j == ny) || (face == "bottom" && j == 0))
        {
          r = sqrt(pow(center[0] - i * dx, 2) + pow(center[2] - k * dz, 2));
          onBoundary = true;
        }
        else if ((face == "front" && k == 0) || (face == "back" && k == nz))
        {
          r = sqrt(pow(center[0] - i * dx, 2) + pow(center[1] - j * dy, 2));
          onBoundary = true;
        }

        if (onBoundary)
        {
          if (r < R)
          {
            value = (2.0 * Q / (PI * R * R)) * (1 - (r * r) / (R * R));
          }
          else if (r >= R)
          {
            value = 0.0;
          }
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
    int in = i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
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
    for (int p : nodeIndices)
    {
      vecTmp.push_back(cell[ic][p]);
    }
    mapCBInCell.push_back(vecTmp);
  };

  switch (inletCB)
  {
  case ControlBoundary::left:
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        collectBoundaryNodes(0, j, k);
      }
    }
    for (int k = 0; k < nz; k++)
    {
      for (int j = 0; j < ny; j++)
      {
        collectBoundaryCells(0, j, k);
        collectBoundaryCellNodes(0, j, k, {0, 3, 7, 4});
      }
    }
    break;

  case ControlBoundary::right:
    for (int k = 0; k < nz + 1; k++)
    {
      for (int j = 0; j < ny + 1; j++)
      {
        collectBoundaryNodes(nx, j, k);
      }
    }
    for (int k = 0; k < nz; k++)
    {
      for (int j = 0; j < ny; j++)
      {
        collectBoundaryCells(nx - 1, j, k);
        collectBoundaryCellNodes(nx - 1, j, k, {1, 5, 6, 2});
      }
    }
    break;

  case ControlBoundary::top:
    for (int k = 0; k < nz + 1; k++)
    {
      for (int i = 0; i < nx + 1; i++)
      {
        collectBoundaryNodes(i, ny, k);
      }
    }
    for (int k = 0; k < nz; k++)
    {
      for (int i = 0; i < nx; i++)
      {
        collectBoundaryCells(i, ny - 1, k);
        collectBoundaryCellNodes(i, ny - 1, k, {2, 6, 7, 3});
      }
    }
    break;

  case ControlBoundary::bottom:
    for (int k = 0; k < nz + 1; k++)
    {
      for (int i = 0; i < nx + 1; i++)
      {
        collectBoundaryNodes(i, 0, k);
      }
    }
    for (int k = 0; k < nz; k++)
    {
      for (int i = 0; i < nx; i++)
      {
        collectBoundaryCells(i, 0, k);
        collectBoundaryCellNodes(i, 0, k, {0, 4, 5, 1});
      }
    }
    break;

  case ControlBoundary::front:
    for (int j = 0; j < ny + 1; j++)
    {
      for (int i = 0; i < nx + 1; i++)
      {
        collectBoundaryNodes(i, j, 0);
      }
    }
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        collectBoundaryCells(i, j, 0);
        collectBoundaryCellNodes(i, j, 0, {0, 1, 2, 3});
      }
    }
    break;

  case ControlBoundary::back:
    for (int j = 0; j < ny + 1; j++)
    {
      for (int i = 0; i < nx + 1; i++)
      {
        collectBoundaryNodes(i, j, nz);
      }
    }
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        collectBoundaryCells(i, j, nz - 1);
        collectBoundaryCellNodes(i, j, nz - 1, {4, 7, 6, 5});
      }
    }
    break;

  default:
    break;
  }
}