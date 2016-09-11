#include <cmath>
#include <iostream>

#include "Boundary.h"
#include "LatticeSite.h"



Topography::Topography(const int d[2], int h, int period,
		       LatticeSite **sites, LatticeSite **_sites)
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
	  sites[x][y].setType(LatticeSite::Solid);
	  _sites[x][y].setType(LatticeSite::Solid);
	}
            
      cc += 1;
    }

}

void Topography::FreeSlipBC(LatticeSite **sites, LatticeSite **_sites)
{
  int x,y;
  int op[9] = {0, 1, 4, 3, 2, 6, 5, 8, 7};
  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];
      _sites[x][y].f[6] = _sites[x][y].f[7];
      _sites[x][y].f[2] = _sites[x][y].f[4];
      _sites[x][y].f[5] = _sites[x][y].f[8];
      // for (int k=0;k<9;k++)
      // 	{
      // 	  _sites[x][y].f[k] = sites[x][y].f[op[k]];
      // 	}
      //  double ux, uy;
      // if(x==25)
      // 	{
      // 	  ux = uy = 0.0;
      // 	  for(int k=0;k<9;k++)
      // 	    {
      // 	      ux += _sites[x][y].f[k]*LatticeSite::e[k][0];
      // 	      uy += _sites[x][y].f[k]*LatticeSite::e[k][1];
      // 	    }
      // 	  std::cout << "ux = " << ux << " | uy = " << uy << std::endl;
      // 	}
    }
  
}
