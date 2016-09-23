#include <cmath>
#include <iostream>

#include "Boundary.h"
#include "LatticeSite.h"



Topography::Topography(const int d[2], int h, double ***u,
		       VelSite **sites, VelSite **_sites)
{
  nbNodes = d[0];
  double u0[2] = {0.0, 0.0};
  int x,y;
  lbl = new int[nbNodes];
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  double a;
  int cc = 0;
  for (int xx=0;xx<d[0];xx++)
    {
      nodes[cc][0] = xx;
      a = 2.*M_PI / d[0];
      nodes[cc][1] = floor(h*sin(a*xx)) + h;
      for(int yy=0;yy<nodes[cc][1];yy++)
	{
	  sites[xx][yy].init(LatticeSite::Solid, 1.0, u0, 0);
	  _sites[xx][yy].init(LatticeSite::Solid, 1.0, u0, 0);
	}
            
      cc += 1;
    }

  for(int i=1;i<nbNodes-1;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];
      if(nodes[i-1][1] == (y-1) && nodes[i+1][1] == (y+1)){lbl[i]=2; } //3 pts aligned along y
      else if(nodes[i-1][1] == nodes[i+1][1]){lbl[i]=1;}
      else if((nodes[i-1][1] == (y-1) && nodes[i+1][1] == y) || (nodes[i-1][1] == y && nodes[i+1][1] == y-1) ){lbl[i] = 2;}
      else{lbl[i] = 0;}
    }

  
  x = nodes[0][0]; y = nodes[0][1];
  if(nodes[nbNodes-1][0] == nodes[1][0]){lbl[0]=2;}
  else if(nodes[nbNodes-1][1] == nodes[1][1]){lbl[0]=1;}
  else if((nodes[nbNodes-1][1] == (y-1) && nodes[1][1] == y) || (nodes[nbNodes-1][1] == y && nodes[1][1] == y-1) ){lbl[0] = 2;}
  else{lbl[0] = 0;}
  
  x = nodes[nbNodes-1][0]; y = nodes[nbNodes-1][1];
  if(nodes[nbNodes-2][0] == nodes[0][0]){lbl[nbNodes-1]=2;}
  else if(nodes[nbNodes-2][1] == nodes[0][1]){lbl[nbNodes-1]=1;}
  else if((nodes[nbNodes-2][1] == (y-1) && nodes[0][1] == y) || (nodes[nbNodes-2][1] == y && nodes[0][1] == y-1) ){lbl[nbNodes-1] = 2;}
  else{lbl[nbNodes-1] = 0;}


  
}

void Topography::BoundaryCondition()
{
}


void Topography::FreeSlipBC(VelSite **sites, VelSite **_sites)
{
  int x,y;
  int op[9] = {0, 1, 4, 3, 2, 6, 5, 8, 7};
  double dd; int cc= 0;

  // ofstream bb, lns, rns, ts;
  // bb.open("bb.dat");lns.open("lns.dat");rns.open("rns.dat");ts.open("ts.dat");
  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];

      if(lbl[i] == 1)
	{
	  _sites[x][y].f[6] = _sites[x][y].f[7];
	  _sites[x][y].f[2] = _sites[x][y].f[4];
	  _sites[x][y].f[5] = _sites[x][y].f[8];
	  // ts << x << " " << y << endl;
	  // cc += 1;
	}
      else if(lbl[i] == 2)
	{
	  if(nodes[i][0] == 0 || nodes[i][0] == nbNodes-1){dd = +1;}
	  else{dd = (nodes[i+1][1]-nodes[i-1][1])/2.;}

	  if(dd > 0.) // if mur gauche
	    {
	      _sites[x][y].f[3] = sites[x][y].f[4];
	      _sites[x][y].f[2] = sites[x][y].f[1];
	      _sites[x][y].f[6] = sites[x][y].f[8];
	      // _sites[x][y].f[6] = _sites[x][y].f[5];
	      // _sites[x][y].f[7] = _sites[x][y].f[8];
	      // _sites[x][y].f[3] = _sites[x][y].f[1];
	      // lns << x << " " << y << endl;
	      // cc += 1;
	    }
	  else
	    {
	      _sites[x][y].f[2] = sites[x][y].f[3];
	      _sites[x][y].f[1] = sites[x][y].f[4];
	      _sites[x][y].f[5] = sites[x][y].f[7];
	      // _sites[x][y].f[5] = _sites[x][y].f[6];
	      // _sites[x][y].f[8] = _sites[x][y].f[7];
	      // _sites[x][y].f[1] = _sites[x][y].f[3];
	      // rns << x << " " << y << endl;
	      // cc += 1;
	    }
	}
      else
	{
	  for (int k=0;k<9;k++)
	    {
	      _sites[x][y].f[k] = sites[x][y].f[op[k]];
	    }
	  // bb << x << " " << y << endl;
	  // cc += 1;
	}
    }

}

void Topography::HalfWayFreeSlipBC(VelSite **sites, VelSite **_sites)
{  int x,y, xm, xp;
  int op[9] = {0, 1, 4, 3, 2, 6, 5, 8, 7};
  double dd; int cc= 0;

  // ofstream bb, lns, rns, ts;
  // bb.open("bb.dat");lns.open("lns.dat");rns.open("rns.dat");ts.open("ts.dat");
  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];
      xp = (x + 1 + nbNodes)%nbNodes; 
      xm = (x - 1 + nbNodes)%nbNodes; 
      //flat wall
      if(lbl[i] == 1)
	{
	  _sites[x][y].f[6] = sites[xm][y].f[7];
	  _sites[x][y].f[2] = sites[x][y].f[4];
	  _sites[x][y].f[5] = sites[xp][y].f[8];
	}

      //External corner
       else if(lbl[i] == 2)
	{
	  if(nodes[i][0] == 0 || nodes[i][0] == nbNodes-1){dd = +1;}
	  else{dd = (nodes[i+1][1]-nodes[i-1][1])/2.;}

	  if(dd > 0.) // if mur gauche
	    {
	      _sites[x][y].f[6] = sites[x][y].f[8];
	    }
	  else
	    {
	      _sites[x][y].f[5] = sites[x][y].f[7];
	    }
	}
      
      // internal corner
       else{
	 if(nodes[i][0] == 0 || nodes[i][0] == nbNodes-1){dd = +1;}
	 else{dd = (nodes[i+1][1]-nodes[i-1][1])/2.;}
	 if(dd > 0.)
	   {
	     _sites[x][y].f[2] = sites[x][y].f[4];
	     _sites[x][y].f[3] = sites[x][y].f[1];
	     _sites[x][y].f[6] = sites[x][y].f[8];
	     _sites[x][y].f[5] = sites[xm][y].f[8];
	   }
	 else
	   {
	     _sites[x][y].f[2] = sites[x][y].f[4];
	     _sites[x][y].f[1] = sites[x][y].f[3];
	     _sites[x][y].f[5] = sites[x][y].f[7];
	     _sites[x][y].f[6] = sites[xp][y].f[7];
	   }
       }
    }
}


void Topography::TemperatureBC(VelSite **velSites, ThermalSite **thermalSites,
				   double **T, double ***u)
{
  int c[4][2] = {{1,0}, {0,1}, {-1,0}, {0, -1}};
  int x,y; 
  double Tcold = 0.0;
  double Tsup, rho, ub[2], usup[2];
  for (int i=0;i<nbNodes;i++)
    {
      double eu;
      x = nodes[i][0]; y = nodes[i][1];

      thermalSites[x][y+1].computeRhoAndU(Tsup);
      velSites[x][y+1].computeRhoAndU(rho, usup);
      velSites[x][y].computeRhoAndU(rho, ub);
      
      for (int k=0;k<4;k++)
	{
	  eu = c[k][0]*ub[0] + c[k][1]*ub[1];
	  thermalSites[x][y].f[k] = (Tsup/4.)*(1.+ 2.*eu) + thermalSites[x][y+1].f[k]
	    - thermalSites[x][y+1].fEq(k,Tsup,usup);
	}
    }
}
