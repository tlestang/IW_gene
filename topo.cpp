#include <cmath>
#include <iostream>

#include "Boundary.h"
#include "LatticeSite.h"



Topography::Topography(const int d[2], int h, double ***u,
		       VelSite **sites, VelSite **_sites)
{
  nbNodes = d[0];

    nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  double a;
  int cc = 0;
  for (int x=0;x<d[0];x++)
    {
      nodes[cc][0] = x;
      a = 2.*M_PI / d[0];
      nodes[cc][1] = floor(h*sin(a*x)) + h;
      
      for(int y=0;y<nodes[cc][1];y++)
	{
	  sites[x][y].setType(VelSite::Solid);
	  _sites[x][y].setType(VelSite::Solid);
	}
            
      cc += 1;
    }

}

void Topography::BoundaryCondition()
{
}

void Topography::FreeSlipBC(VelSite **sites, VelSite **_sites)
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

void Topography::TemperatureBC(VelSite **velSites, ThermalSite **thermalSites,
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
