#include <cmath>
#include <iostream>

#include "Boundary.h"
#include "LatticeSite.h"



Topography::Topography(const int d[2], int h, int period,
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
      a = 2.*M_PI / period;
      nodes[cc][1] = floor(h*sin(a*x)) + h;
      
      for(int y=0;y<nodes[cc][1];y++)
	{
	  sites[x][y].setType(VelSite::Solid);
	  _sites[x][y].setType(VelSite::Solid);
	}
            
      cc += 1;
    }

}

void Topography::FreeSlipBC(VelSite **sites, VelSite **_sites)
{
  int x,y;
  int op[9] = {0, 1, 4, 3, 2, 6, 5, 8, 7};
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
  int x,y; double u0[2] = {0.0, 0.0};
  double Tcold = 0.0;
  double TT, pp, uu[2];
  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];

      thermalSites[x][y+1].computeRhoAndU(TT);
      velSites[x][y+1].computeRhoAndU(pp, uu);
      
      for (int k=0;k<4;k++)
	{
	  thermalSites[x][y].f[k] = TT/4. + thermalSites[x][y+1].f[k]
	    - thermalSites[x][y+1].fEq(k,TT,uu);
	}
    }
}
