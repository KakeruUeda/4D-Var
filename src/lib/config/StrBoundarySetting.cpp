/**
 * @file SetStrBoudnary.cpp
 * @author K.Ueda
 * @date August, 2024
 */

#include "Config.h"

void Config::setBoundaryVelocityValue(std::string face, double value[3])
{
  std::vector<double> vecValue(3);
  for(int i = 0; i < 3; i++)
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
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
          int in = k * (ny+1) * (nx+1) + j * (nx+1) + i;
          if (k == ny)
          {
            pDirichlet[in] = value;
          }
        }
      }
    }
  }
}