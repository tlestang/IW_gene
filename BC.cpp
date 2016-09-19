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
  int c[4][2] = {{1,0}, {0,1}, {-1,0}, {0, -1}};
    for(int x=0;x<dims[0];x++)
    {
      for (int y=0;y<dims[1];y++)
	{
	  if(velSites[x][y].isFluidSolid())
	    {
	      
	      // -------- Velocity Boundary Conditions -------------
	      
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

	      // ----- Temperature Boundary Condition

	      int yp = (y + 1 + dims[1])%dims[1];
	      double Tcold = 0.0; double eu;
	      double Tinf, Tsup, rho, ub[2], uinf[2], usup[2];
	      
	      if(velSites[x][yp].isSolid())
		{
		  thermalSites[x][y-1].computeRhoAndU(Tinf);
		  velSites[x][y-1].computeRhoAndU(rho, uinf);
		  velSites[x][y].computeRhoAndU(rho, ub);
		  
		  for (int k=0;k<4;k++)
		    {
		      eu = c[k][0]*ub[0] + c[k][1]*ub[1];
		      thermalSites[x][y].f[k] = (Tinf/4.)*(1.+2.*eu)
			+ thermalSites[x][y-1].f[k]
			- thermalSites[x][y-1].fEq(k,Tinf,uinf);
		    }
		}
	      else //bottom wall (topography)
		{
		  thermalSites[x][y+1].computeRhoAndU(Tsup);
		  velSites[x][y+1].computeRhoAndU(rho, usup);
		  velSites[x][y].computeRhoAndU(rho, ub);
      
		  for (int k=0;k<4;k++)
		    {
		      eu = c[k][0]*ub[0] + c[k][1]*ub[1];
		      thermalSites[x][y].f[k] = (Tsup/4.)*(1.+ 2.*eu)
			+ thermalSites[x][y+1].f[k]
			- thermalSites[x][y+1].fEq(k,Tsup,usup);
		    }
	
		}
	    }
	}
    }
}
