#include <iostream>
#include "Boundary.h"
#include "LatticeSite.h"

Boundary::Boundary()
{
}

TopWall::TopWall(const int d[2])
{
  nbNodes = d[0];
  dims[0] = d[0]; dims[1] = d[1];
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  int cc = 0;

  for (int x=0;x<d[0];x++)
    {
      nodes[cc][0] = x; nodes[cc][1] = d[1]-1;
      cc += 1;
    }
}

void TopWall::BoundaryCondition()
{
}
void TopWall::BoundaryCondition(VelSite **sites, VelSite **_sites)
{
  int x,y;
  int op[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];
      for (int k=0;k<9;k++)
	{
	  _sites[x][y].f[k] = sites[x][y].f[op[k]];
	}
    }
}

void TopWall::FreeSlipBC(VelSite **sites, VelSite **_sites)
{
  int x,y; int xm, xp;
  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];
      xp = (x + 1 + dims[0])%dims[0];
      xm = (x - 1 + dims[0])%dims[0];
      _sites[x][y].f[7] = sites[xp][y].f[6];
      _sites[x][y].f[4] = sites[x][y].f[2];
      _sites[x][y].f[8] = sites[xm][y].f[5];
    }
}

void TopWall::TemperatureBC(VelSite **velSites, ThermalSite **thermalSites,
				double **T, double ***u)
{
  int c[4][2] = {{1,0}, {0,1}, {-1,0}, {0, -1}};
  int x,y; 
  double Tcold = 0.0;
  double Tinf, rho, ub[2], uinf[2];
  for (int i=0;i<nbNodes;i++)
    {
      double eu;
      x = nodes[i][0]; y = nodes[i][1];

      thermalSites[x][y-1].computeRhoAndU(Tinf);
      velSites[x][y-1].computeRhoAndU(rho, uinf);
      velSites[x][y].computeRhoAndU(rho, ub);
      
      for (int k=0;k<4;k++)
	{
	  eu = c[k][0]*ub[0] + c[k][1]*ub[1];
	  thermalSites[x][y].f[k] = (Tinf/4.)*(1.+2.*eu) + thermalSites[x][y-1].f[k]
	    - thermalSites[x][y-1].fEq(k,Tinf,uinf);
	}
    }
}

BottomWall::BottomWall(const int d[2])
{
  nbNodes = d[0];
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  int cc = 0;

  for (int x=0;x<d[0];x++)
    {
      nodes[cc][0] = x; nodes[cc][1] = 0;
      cc += 1;
    }
}

void BottomWall::BoundaryCondition()
{
}
void BottomWall::BoundaryCondition(VelSite **sites, VelSite **_sites)
{
  int x,y;
  int op[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];
      for (int k=0;k<9;k++)
	{
	  _sites[x][y].f[k] = sites[x][y].f[op[k]];
	}
    }
}

void BottomWall::FreeSlipBC(VelSite **sites, VelSite **_sites)
{
  int x,y;
  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];
      _sites[x][y].f[6] = _sites[x][y].f[7];
      _sites[x][y].f[2] = _sites[x][y].f[4];
      _sites[x][y].f[5] = _sites[x][y].f[8];
    }
}

void BottomWall::TemperatureBC(VelSite **velSites, ThermalSite **thermalSites,
				   double **T, double ***u)
{
  int c[4][2] = {{1,0}, {0,1}, {-1,0}, {0, -1}};
  int x,y; 
  double Tcold = 0.0;
  double Tsup, rho, ub[2], usup[2];
  for (int i=0;i<nbNodes;i++)
    {
      double eu;
      x = nodes[i][0]; y = nodes[i][1];

      thermalSites[x][y+1].computeRhoAndU(Tsup);
      velSites[x][y+1].computeRhoAndU(rho, usup);
      velSites[x][y].computeRhoAndU(rho, ub);
      
      for (int k=0;k<4;k++)
	{
	  eu = c[k][0]*ub[0] + c[k][1]*ub[1];
	  thermalSites[x][y].f[k] = (Tsup/4.)*(1.+ 2.*eu) + thermalSites[x][y+1].f[k]
	    - thermalSites[x][y+1].fEq(k,Tsup,usup);
	}
    }
}
