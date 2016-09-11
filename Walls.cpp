#include <iostream>
#include "Boundary.h"
#include "LatticeSite.h"

Boundary::Boundary()
{
}

TopWall::TopWall(const int d[2])
{
  nbNodes = d[0];
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

void TopWall::TemperatureBC(VelSite **velSites, ThermalSite **thermalSites,
				double **T, double ***u)
{
  int x,y; double u0[2] = {0.0, 0.0};
  double Tcold = 0.0;
  double TT, pp, uu[2];
  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];

      thermalSites[x][y-1].computeRhoAndU(TT);
      velSites[x][y-1].computeRhoAndU(pp, uu);
      
      for (int k=0;k<4;k++)
	{
	  thermalSites[x][y].f[k] = TT/4. + thermalSites[x][y-1].f[k]
	    - thermalSites[x][y-1].fEq(k,TT,uu);
	}
    }
}


