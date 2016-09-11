#include <iostream>

#include "Boundary.h"
#include "LatticeSite.h"
using namespace std;
Square::Square(int L, int point[2])
{
  
  xmin = point[0]; xmax = xmin + L;
  ymin = point[1]; ymax = ymin + L;

}

void Square::BoundaryCondition(LatticeSite **sites)
{
}

void Square::BoundaryCondition(LatticeSite **sites, LatticeSite **_sites)
{
  int op[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
  for(int k=0;k<9;k++)
    {
      for(int x=xmin;x<xmax+1;x++)
      	{
	  //cout << sites[x][ymax].f[op[k]] << endl;
      	  _sites[x][ymax].f[k] = sites[x][ymax].f[op[k]];
      	  _sites[x][ymin].f[k] = sites[x][ymin].f[op[k]];
      	}
            for(int y=ymin;y<ymax+1;y++)
      	{
      	  _sites[xmax][y].f[k] = sites[xmax][y].f[op[k]];
      	  _sites[xmin][y].f[k] = sites[xmin][y].f[op[k]];
      	}
    }
}
