#include <iostream>
#include <cmath>

#include "LatticeBoltzmann.h"
#include "LatticeSite.h"
#include "Boundary.h"

using namespace std;

LatticeBoltzmann::LatticeBoltzmann(const int d[2], const double omega_[2],
				   double coef_force_)
{

  dims[0] = d[0]; dims[1] = d[1];
  omega[0] = omega_[0];
  omega[1] = omega_[1];
  coef_force = coef_force_;

  velSites = new VelSite*[dims[0]];
  velSites_ = new VelSite*[dims[0]];

  for (int i=0;i<dims[0];i++)
    {
      velSites[i] = new VelSite[dims[1]];
      velSites_[i] = new VelSite[dims[1]];
    }

  thermalSites = new ThermalSite*[dims[0]];
  thermalSites_ = new ThermalSite*[dims[0]];

  for (int i=0;i<dims[0];i++)
    {
      thermalSites[i] = new ThermalSite[dims[1]];
      thermalSites_[i] = new ThermalSite[dims[1]];
    }
  
  // allocate memory for T and rho
  T = new double*[dims[0]];
  for (int x=0; x<dims[0]; x++)
    T[x] = new double[dims[1]];

    rho = new double*[dims[0]];
  for (int x=0; x<dims[0]; x++)
    rho[x] = new double[dims[1]];
	
  // allocate memory for velocity
  u = new double**[dims[0]];
  for (int x=0; x<dims[0]; x++)
    {
      u[x] = new double*[dims[1]];
		
      for (int y=0; y<dims[1]; y++)
	u[x][y] = new double[2];
    }

  generateGeometry();
  
  w = new TopWall(d);
  topo = new Topography(d, floor(dims[1]/10), dims[0], velSites, velSites_);
  

}


void LatticeBoltzmann::streamToNeighbors(int x, int y)
{
	for (int k=0; k<9; k++)
	{
	  int nx = (x + VelSite::e[k][0] + dims[0])%dims[0];
	  int ny = (y + VelSite::e[k][1] + dims[1])%dims[1];

	  velSites_[nx][ny].f[k] = velSites[x][y].f[k];
	}
	for (int k=0;k<4;k++)
	  {
	    int nx = (x + ThermalSite::c[k][0] + dims[0])%dims[0];
	    int ny = (y + ThermalSite::c[k][1] + dims[1])%dims[1];
	    
	    thermalSites_[nx][ny].f[k] = thermalSites[x][y].f[k];
	  }
}

void LatticeBoltzmann::update()
{
  ThermalSite **swapT;
  VelSite **swapVel;

	for (int x=0; x<dims[0]; x++)
	{
		for (int y=0; y<dims[1]; y++)
		{
		  if(velSites[x][y].isFluid())
		    {
		      
		  velSites[x][y].computeRhoAndU(rho[x][y], u[x][y]);
		  thermalSites[x][y].computeRhoAndU(T[x][y]);
		  
		  velSites[x][y].collide(rho[x][y], T[x][y], u[x][y]);
		  thermalSites[x][y].collide(T[x][y], u[x][y]);
		  
		  streamToNeighbors(x, y);
		    }
		}
	}

	w->BoundaryCondition(velSites, velSites_);
	topo->FreeSlipBC(velSites, velSites_);

	w->TemperatureBC(velSites_, thermalSites_, T, u);
	topo->TemperatureBC(velSites_, thermalSites_, T, u);


	swapT = thermalSites;
	thermalSites = thermalSites_;
	thermalSites_ = swapT;

	swapVel = velSites;
	velSites = velSites_;
	velSites_ = swapVel;
}

void LatticeBoltzmann::generateGeometry()
{
	double u[2] = {0, 0};

	for (int x=0; x<dims[0]; x++)
	{
		for (int y=0; y<dims[1]; y++)
		{
			if (x == 0 || x == dims[0]-1 || y == 0 || y == dims[1]-1)
			{
			  velSites[x][y].init(LatticeSite::Fluid, 1.0, u,
					      omega[0], coef_force);
			  velSites_[x][y].init(LatticeSite::Fluid, 1.0, u,
					       omega[0], coef_force);
			}
			else
			{
			  velSites[x][y].init(LatticeSite::Fluid, 1.0, u,
					      omega[0], coef_force);
			  velSites_[x][y].init(LatticeSite::Fluid, 1.0, u,
					       omega[0], coef_force);
			}
			if (y==0)
			  {
			    thermalSites[x][y].init(LatticeSite::Fluid, 1.0, u,
						    omega[1], coef_force);
			    thermalSites_[x][y].init(LatticeSite::Fluid, 1.0, u,
						     omega[1], coef_force);
			    T[x][y] = 1.0;
			  }
			else if(x==(dims[0]-1)/2 && y == 1)
			  {
			    T[x][y] = 1.1;
			    thermalSites[x][y].init(LatticeSite::Fluid, 1.1, u,
						    omega[1], coef_force);
			    thermalSites_[x][y].init(LatticeSite::Fluid, 1.1, u,
						     omega[1], coef_force);
			}
			else
			  {
			    T[x][y]=0.0;
			    thermalSites[x][y].init(LatticeSite::Fluid, 0.0, u,
						    omega[1], coef_force);
			    thermalSites_[x][y].init(LatticeSite::Fluid, 0.0, u,
						     omega[1], coef_force);
			  }
		}
	}
}

void LatticeBoltzmann::getDensityAndVelocityField(double **&tp,
						  double **&rp, double ***&up)
{
  tp = T;
  rp = rho;
  up = u;
}
	     
LatticeBoltzmann::~LatticeBoltzmann()
{
  delete thermalSites;
  delete thermalSites_;
  delete velSites;
  delete velSites_;
  
  delete w;
  delete topo;
}
