#include <iostream>
#include "Boundary.h"
#include "LatticeSite.h"

Boundary::Boundary()
{
}

Walls::Walls(const int d[2])
{
  nbNodes = 2 * d[0];//  + 2* d[1];
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  int cc = 0;
  for (int x=0;x<d[0];x++)
    {
      nodes[cc][0] = x; nodes[cc][1] = 0;
      cc += 1;
    }

  for (int x=0;x<d[0];x++)
    {
      nodes[cc][0] = x; nodes[cc][1] = d[1]-1;
      cc += 1;
    }
  // for (int y=0;y<d[1];y++)
  //   {
  //     nodes[cc][0] = 0; nodes[cc][1] = y;
  //     cc +=1;
  //   }
  //   for (int y=0;y<d[1];y++)
  //   {
  //     nodes[cc][0] = d[0]-1; nodes[cc][1] = y;
  //     cc +=1;
  //   }
}
void Walls::BoundaryCondition()
{
}
void Walls::BoundaryCondition(VelSite **sites, VelSite **_sites)
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


