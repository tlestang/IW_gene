#include <iostream>
#include "Boundary.h"
#include "LatticeSite.h"

HotWall::HotWall(const int d[2])
{

  nbNodes = d[0];
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  int cc = 0;
  for (int x=0;x<d[0];x++)
    {
      nodes[cc][0] = x; nodes[cc][1] = 0;
      cc += 1;
    }

}
void HotWall::BoundaryCondition()
{
}

// void HotWall::BoundaryCondition(VelSite **velSites, ThermalSite **thermalSites,
// 				double **T, double ***u)
// {
//   int x,y; double u0[2] = {0.0, 0.0};
//   double Thot = 1.0;
//   double TT, pp, uu[2];
//   for (int i=0;i<nbNodes;i++)
//     {
//       x = nodes[i][0]; y = nodes[i][1];
      
//       thermalSites[x+1][y].computeRhoAndU(TT);
//       velSites[x+1][y].computeRhoAndU(pp, uu);
      
//       for (int k=0;k<4;k++)
//   	{
//   	  thermalSites[x][y].f[k] = thermalSites[x][y].fEq(k,Thot,u0) + thermalSites[x+1][y].f[k]
//   	    - thermalSites[x+1][y].fEq(k,TT,uu);
//   	}
//     }
//   }

void HotWall::BoundaryCondition(ThermalSite **sites)
{
  int x,y;
  double Thot = 1.0;

  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];
      sites[x][y].f[1] = Thot-sites[x][y].f[0]-sites[x][y].f[2]-sites[x][y].f[3];
    }
}

ColdWall::ColdWall(const int d[2])
{

  nbNodes = d[0];
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  int cc = 0;
  for (int x=0;x<d[0];x++)
    {
      nodes[cc][0] = x; nodes[cc][1] = d[1]-1;
      cc += 1;
    }

}
void ColdWall::BoundaryCondition()
{
}
// void ColdWall::BoundaryCondition(VelSite **velSites, ThermalSite **thermalSites, double **T, double ***u)
// {
//   int x,y; double u0[2] = {0.0, 0.0};
//   double Tcold = 1.0;
//   double TT, pp, uu[2];
//   for (int i=0;i<nbNodes;i++)
//     {
//       x = nodes[i][0]; y = nodes[i][1];

//       thermalSites[x-1][y].computeRhoAndU(TT);
//       velSites[x-1][y].computeRhoAndU(pp, uu);
      
//       for (int k=0;k<4;k++)
// 	{
// 	  thermalSites[x][y].f[k] = thermalSites[x][y].fEq(k,Tcold,u0) + thermalSites[x-1][y].f[k]
// 	    - thermalSites[x-1][y].fEq(k,TT,uu);
// 	}
//     }
  
// }
void ColdWall::BoundaryCondition(ThermalSite **sites)
{
  int x,y;
  double Tcold = 0.0;

  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];
      sites[x][y].f[3] = Tcold-sites[x][y].f[0]-sites[x][y].f[1]-sites[x][y].f[2];
    }
}

TopWall::TopWall(const int d[2])
{

  nbNodes = d[0]-2;
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  int cc = 0;
  for (int x=1;x<d[0]-1;x++)
    {
      nodes[cc][0] = x; nodes[cc][1] = d[1]-1;
      cc += 1;
    }

}
void TopWall::BoundaryCondition()
{
}
void TopWall::BoundaryCondition(VelSite **velSites, ThermalSite **thermalSites,
				double **T, double ***u)
{
  int x,y; double u0[2] = {0.0, 0.0};
  double Tcold = 0.0;
  double TT, pp, uu[2];
  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];

      thermalSites[x][y-1].computeRhoAndU(TT);
      velSites[x][y-1].computeRhoAndU(pp, uu);
      
      for (int k=0;k<4;k++)
	{
	  thermalSites[x][y].f[k] = TT/4. + thermalSites[x][y-1].f[k]
	    - thermalSites[x][y-1].fEq(k,TT,uu);
	}
    }
}

  BottomWall::BottomWall(const int d[2])
{

  nbNodes = d[0]-2;
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  int cc = 0;
  for (int x=1;x<d[0]-1;x++)
    {
      nodes[cc][0] = x; nodes[cc][1] = 0;
      cc += 1;
    }

}
void BottomWall::BoundaryCondition()
{
}
void BottomWall::BoundaryCondition(VelSite **velSites, ThermalSite **thermalSites,
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

