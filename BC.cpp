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
  int nlink, nx, ny;

  for(int x=1;x<dims[0]-1;x++)
    {

      for (int y=1;y<dims[1]-1;y++)
	{
	  
	  if(velSites[x][y].isFluidSolid())
	    {
	      // -------- Velocity Boundary Conditions -------------

	      nlink = velSites[x][y].getNormalLink();
	      for(int k=0;k<9;k++)
		{
		  if(velSites[x][y].brLinks[k] == 1)
		    {
		      if(k < 5 || k == nlink) // BB if link normal to surface
			{
			  nx = x + VelSite::e[k][0];
			  ny = y + VelSite::e[k][1];
			  velSites_[x][y].f[op[k]] = velSites_[nx][ny].f[k];
			}
		      else
			{
			  switch (k) {
			  case 5:
			    if (velSites_[x][y+1].isFluid())
			      velSites_[x][y].f[op[k]] = velSites_[x+1][y].f[8];
			    else
			      velSites_[x][y].f[op[k]] = velSites_[x][y+1].f[6];
			    break;
			  case 6:
			    if (velSites_[x][y+1].isFluid())
			      velSites_[x][y].f[op[k]] = velSites_[x-1][y].f[7];
			    else
			      velSites_[x][y].f[op[k]] = velSites_[x][y+1].f[5];
			    break;
			  case 7:
			    if (velSites_[x][y-1].isFluid())
			      velSites_[x][y].f[op[k]] = velSites_[x-1][y].f[6];
			    else
			      velSites_[x][y].f[op[k]] = velSites_[x][y-1].f[8];
			    break;
			  case 8:
			    if (velSites_[x][y-1].isFluid())
			      velSites_[x][y].f[op[k]] = velSites_[x+1][y].f[5];
			    else
			      velSites_[x][y].f[op[k]] = velSites_[x][y-1].f[7];
			    break;
			  }
			}
		    }
		}

	      // ----- Temperature Boundary Condition


	      double Tcold = 0.0; double eu;
	      double Tinf, Tsup, rho, ub[2], uinf[2], usup[2];
	      
	      if(y == dims[1]-2)
	      	{
	      	  thermalSites_[x][y-1].computeRhoAndU(Tinf);
	      	  velSites_[x][y-1].computeRhoAndU(rho, uinf);
	      	  velSites_[x][y].computeRhoAndU(rho, ub);
		  
	      	  for (int k=0;k<4;k++)
	      	    {
	      	      eu = c[k][0]*ub[0] + c[k][1]*ub[1];
	      	      thermalSites_[x][y].f[k] = (Tinf/4.)*(1.+2.*eu)
	      		+ thermalSites_[x][y-1].f[k]
	      		- thermalSites_[x][y-1].fEq(k,Tinf,uinf);
	      	    }
	      	}
	      else //bottom wall (topography)
	      	{

	      	  thermalSites_[x][y+1].computeRhoAndU(Tsup);
	      	  velSites_[x][y+1].computeRhoAndU(rho, usup);
	      	  velSites_[x][y].computeRhoAndU(rho, ub);
      
	      	  for (int k=0;k<4;k++)
	      	    {
	      	      eu = c[k][0]*ub[0] + c[k][1]*ub[1];
	      	      thermalSites_[x][y].f[k] = (Tsup/4.)*(1.+ 2.*eu)
	      		+ thermalSites_[x][y+1].f[k]
	      		- thermalSites_[x][y+1].fEq(k,Tsup,usup);
	      	    }
	
	      	}
	    }
	}
    }

}
