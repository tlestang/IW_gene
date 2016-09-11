#include "Boundary.h"
#include "LatticeSite.h"

Outlet::Outlet(const int d[2])
{
  nbNodes = d[1]-2;
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  int cc = 0;

  for (int y=1;y<d[1]-1;y++)
    {
      nodes[cc][0] = d[0]-1; nodes[cc][1] = y;
      cc += 1;
    }

}

void Outlet::BoundaryCondition(LatticeSite **sites, double ***u)
{
  int x,y;
    
  for (int i=0;i<nbNodes;i++)
    {

      x = nodes[i][0]; y = nodes[i][1];
      u0[0] = u[x-1][y][0];
      //u0[0] = (1./3.)*(4.0*u[x-1][y][0]-u[x-2][y][0]);
      //u0[1] = (1./3.)*(4.0*u[x-1][y][1]-u[x-2][y][1]);
      u0[1] = 0.0;
      
      double a0 = 1./(1.0+u0[0]);
      // Computes local density
      rho = 2.0*a0*(sites[x][y].f[5]+sites[x][y].f[1]+sites[x][y].f[8])
	    + a0*(sites[x][y].f[2]+sites[x][y].f[4]+sites[x][y].f[0]);
      // Computes out-of-equilibrium parts

      for  (int k=0;k<9;k++)
	{
	  fneq[k] = sites[x][y].f[k] - sites[x][y].fEq(k, rho, u0);
	}

      fneq[6] = fneq[8] + 0.5*(fneq[4]-fneq[2]);
      fneq[7] = fneq[5] + 0.5*(fneq[2]-fneq[4]);
      fneq[3] = fneq[1];

      sites[x][y].f[6] = sites[x][y].fEq(6, rho, u0) + fneq[6];
      sites[x][y].f[7] = sites[x][y].fEq(7, rho, u0) + fneq[7];
      sites[x][y].f[3] = sites[x][y].fEq(3, rho, u0) + fneq[3];

    }
}

      
void Outlet::BoundaryCondition(LatticeSite **sites)
{
}
