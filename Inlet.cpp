#include <cmath>

#include "Boundary.h"
#include "LatticeSite.h"

using namespace std;

Inlet::Inlet(int nozzleLength, const int d[2], const double u0_[2])
{

  u0[0] = u0_[0]; u0[1] = u0_[1];
  
  nbNodes = nozzleLength + 1;
  int nozzleStart = floor((d[1] - (nozzleLength + 1))*0.5);
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  int cc = 0;

  for (int y=nozzleStart;y<nozzleStart+nbNodes;y++)
    {
      nodes[cc][0] = 0; nodes[cc][1] = y;
      cc += 1;
    }
}

void Inlet::BoundaryCondition(LatticeSite **sites)
{
    int x,y;

  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];

      double a0 = 1./(1.0-u0[0]);
      // Computes local density
      rho = 2.0*a0*(sites[x][y].f[6]+sites[x][y].f[3]+sites[x][y].f[7])
	    + a0*(sites[x][y].f[2]+sites[x][y].f[4]+sites[x][y].f[0]);
      
      // Computes out-of-equilibrium parts
      for  (int k=0;k<9;k++)
	{
	  fneq[k] = sites[x][y].f[k] - sites[x][y].fEq(k, rho, u0);
	}

      fneq[8] = fneq[6] + 0.5*(fneq[2]-fneq[4]);
      fneq[5] = fneq[7] + 0.5*(fneq[4]-fneq[2]);
      fneq[1] = fneq[3];

      sites[x][y].f[8] = sites[x][y].fEq(8, rho, u0) + fneq[8];
      sites[x][y].f[5] = sites[x][y].fEq(5, rho, u0) + fneq[5];
      sites[x][y].f[1] = sites[x][y].fEq(1, rho, u0) + fneq[1];
    }
}

      
