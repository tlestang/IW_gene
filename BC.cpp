#include <iostream>
#include <fstream>
#include <cmath>

#include "LatticeBoltzmann.h"
#include "LatticeSite.h"
#include "Boundary.h"

using namespace std;

void LatticeBoltzmann::BoundaryConditions()
{
  int op[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    for(int x=0;x<dims[0];x++)
    {
      for (int y=0;y<dims[1];y++)
	{
	  if(velSites[x][y].isFluidSolid())
	    {
	      nlink = velSites[x][y].getNormalLink();
	      for(int k=0;k<9;k++)
		{
		  if(k < 5 || k == nlink) // BB if link normal to surface
		    {
		      nx = x + LatticeNode::e[k][0];
		      ny = y + LatticeNode::e[k][1];
		      velSites_[x][y].f[op[k]] = velSites[nx][ny].f[k];
		    }
		  else
		    {
		      				switch (k) {
				case 5:
					if (nodes_s_[x][y + 1].isFluid())
						nodes_s_[x][y].f[op[k]] = nodes_s_[x + 1][y].f[8];
					else
						nodes_s_[x][y].f[op[k]] = nodes_s_[x][y + 1].f[6];
					break;
				case 6:
					if (nodes_s_[x][y + 1].isFluid())
						nodes_s_[x][y].f[op[k]] = nodes_s_[x - 1][y].f[7];
					else
						nodes_s_[x][y].f[op[k]] = nodes_s_[x][y + 1].f[5];
					break;
				case 7:
					if (nodes_s_[x][y - 1].isFluid())
						nodes_s_[x][y].f[op[k]] = nodes_s_[x - 1][y].f[6];
					else
						nodes_s_[x][y].f[op[k]] = nodes_s_[x][y - 1].f[8];
					break;
				case 8:
					if (nodes_s_[x][y - 1].isFluid())
						nodes_s_[x][y].f[op[k]] = nodes_s_[x + 1][y].f[5];
					else
						nodes_s_[x][y].f[op[k]] = nodes_s_[x][y - 1].f[7];
					break;
						}
		    }
		}
	    }
	}
    }
}
