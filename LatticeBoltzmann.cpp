#include <iostream>
#include <fstream>
#include <cmath>

#include "LatticeBoltzmann.h"
#include "LatticeSite.h"
#include "Boundary.h"

using namespace std;

LatticeBoltzmann::LatticeBoltzmann(const int d[2], const double omega_[2],
				   double coef_force_, double u0_, int h_, double N2_)
{

  dims[0] = d[0]; dims[1] = d[1];
  omega[0] = omega_[0];
  omega[1] = omega_[1];
  coef_force = coef_force_;
  u0 = u0_; N2 = N2_; h = h_;

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
  topo = new Topography(d, h, u, velSites, velSites_);
  

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
  double om, a, dSpge;
  //ofstream lulu("spgegeometry.dat");
  	for (int x=0; x<dims[0]; x++)
	{
		for (int y=0; y<dims[1]; y++)
		{
		  if(velSites[x][y].isFluid())
		    {
		      
		  velSites[x][y].computeRhoAndU(rho[x][y], u[x][y]);
		  thermalSites[x][y].computeRhoAndU(T[x][y]);

		  
		  if(y>(spgeFirstNode-1))
		    {
		      if(y>ySpge){om = 0.001*omega[0];}
		      else{
			dSpge = ySpge-spgeFirstNode;
			a = (y-spgeFirstNode)/dSpge;
			om = (1.-0.999*a*a)*omega[0];
		      }
		    }
		  else{om = omega[0];}

		  // if(x==250)
		  //   {
		  //     lulu << y << " " << om << endl;
		  //   }
		  
		  velSites[x][y].collide(rho[x][y], T[x][y], u[x][y], om);
		  thermalSites[x][y].collide(rho[x][y], T[x][y], u[x][y], omega[1]);
		  
		  streamToNeighbors(x, y);
		    }
		}

	}

	w->FreeSlipBC(velSites, velSites_);
	topo->FreeSlipBC(velSites, velSites_);
	//wb->FreeSlipBC(velSites, velSites_);

	w->TemperatureBC(velSites_, thermalSites_, T, u);
	topo->TemperatureBC(velSites_, thermalSites_, T, u);
	//wb->TemperatureBC(velSites_, thermalSites_, T, u);


	swapT = thermalSites;
	thermalSites = thermalSites_;
	thermalSites_ = swapT;

	swapVel = velSites;
	velSites = velSites_;
	velSites_ = swapVel;
}

void LatticeBoltzmann::generateGeometry()
{
	double u[2] = {u0, 0};
	double a = N2/coef_force; double TT;
	
	for (int x=0; x<dims[0]; x++)
	{
		for (int y=0; y<dims[1]; y++)
		{
		  //u[0] = InitialCondition_X(x,y);
		  //u[1] = InitialCondition_Y(x,y);
		  
		  velSites[x][y].init(LatticeSite::Fluid, 1.0, u,
				      coef_force);
		  velSites_[x][y].init(LatticeSite::Fluid, 1.0, u,
				       coef_force);
		  
		  TT = a*y;
		  T[x][y]=TT;
		  thermalSites[x][y].init(LatticeSite::Fluid, TT, u,
					  coef_force);
		  thermalSites_[x][y].init(LatticeSite::Fluid, TT, u,
					   coef_force);
		}
	}
}

void LatticeBoltzmann::TagFluidNodes()
{
  int nx, ny;
  
  for(int x=0;x<dims[0];x++)
    {
      for (int y=0;y<dims[1];y++)
	{
	  
	  if(velSites.[x][y].isFluid())
	    {
	      for(int k=0;k<9;k++)
		{
		  nx = (x + VelSite::e[k][0] + dims[0])%dims[0];
		  ny = (y + VelSite::e[k][1] + dims[1])%dims[1];

		  if(velSites.[nx][ny].isSolid())
		    {
		      if(!velSites[x][y].isFluidSolid())
			{
			  velSites.setFluidTag(LatticeSite::FluidSolid);
			  velSites_.setFluidTag(LatticeSite::FluidSolid);
			}
		      // Matk the node as Boundary node (fluidSOlid)
		      brLinks[k] = 1;
		    }
		  else{
		    brLinks = 0;
		  }

		  if (!isSet) {
		    switch (k) {
		    case 0:
		      break;

		    case 1:
		      if (velSites[x + LatticeSite::e[5][0]][y+ LatticeSite::e[5][1]].isSolid()
			  && velSites[x + LatticeSite::e[8][0]][y + LatticeSite::e[8][1]].isSolid()
			  && velSites[x + LatticeSite::e[2][0]][y + LatticeSite::e[2][1]].isFluid()
			  && velSites[x + LatticeSite::e[4][0]][y + LatticeSite::e[4][1]].isFluid()) {
			isSet = true;
			velSites[x][y].setNormalLink(k);
			velSites_[x][y].setNormalLink(k);
		      }
		      break;
		    case 2:
		      if (velSites[x + LatticeSite::e[5][0]][y + LatticeSite::e[5][1]].isSolid()
			  && velSites[x + LatticeSite::e[6][0]][y + LatticeSite::e[6][1]].isSolid()
			  && velSites[x + LatticeSite::e[1][0]][y + LatticeSite::e[1][1]].isFluid()
			  && velSites[x + LatticeSite::e[3][0]][y + LatticeSite::e[3][1]].isFluid()) {
			isSet = true;
			velSites[x][y].setNormalLink(k);
			velSites_[x][y].setNormalLink(k);
		      }
		      break;
		    case 3:
		      if (velSites[x + LatticeSite::e[6][0]][y + LatticeSite::e[6][1]].isSolid()
			  && velSites[x + LatticeSite::e[7][0]][y + LatticeSite::e[7][1]].isSolid()
			  && velSites[x + LatticeSite::e[2][0]][y + LatticeSite::e[2][1]].isFluid()
			  && velSites[x + LatticeSite::e[4][0]][y + LatticeSite::e[4][1]].isFluid()) {
			isSet = true;
			velSites[x][y].setNormalLink(k);
			velSites_[x][y].setNormalLink(k);
		      }
		      break;
		    case 4:
		      if (velSites[x + LatticeSite::e[7][0]][y + LatticeSite::e[7][1]].isSolid()
			  && velSites[x + LatticeSite::e[8][0]][y + LatticeSite::e[8][1]].isSolid()
			  && velSites[x + LatticeSite::e[1][0]][y + LatticeSite::e[1][1]].isFluid()
			  && velSites[x + LatticeSite::e[3][0]][y + LatticeSite::e[3][1]].isFluid()) {
			isSet = true;
			velSites[x][y].setNormalLink(k);
			velSites_[x][y].setNormalLink(k);
		      }
		      break;
		    case 5:
		      if ((velSites[x + LatticeSite::e[1][0]][y + LatticeSite::e[1][1]].isFluid()
			   && velSites[x + LatticeSite::e[2][0]][y + LatticeSite::e[2][1]].isFluid())
			  || (velSites[x + LatticeSite::e[1][0]][y + LatticeSite::e[1][1]].isSolid()
			      && velSites[x + LatticeSite::e[2][0]][y + LatticeSite::e[2][1]].isSolid())) {
			isSet = true;
			velSites[x][y].setNormalLink(k);
			velSites_[x][y].setNormalLink(k);
		      }
		      break;
		    case 6:
		      if ((velSites[x + LatticeSite::e[3][0]][y + LatticeSite::e[3][1]].isFluid()
			   && velSites[x + LatticeSite::e[2][0]][y + LatticeSite::e[2][1]].isFluid())
			  || (velSites[x + LatticeSite::e[3][0]][y + LatticeSite::e[3][1]].isSolid()
			      && velSites[x + LatticeSite::e[2][0]][y + LatticeSite::e[2][1]].isSolid())) {
			isSet = true;
			velSites[x][y].setNormalLink(k);
			velSites_[x][y].setNormalLink(k);
		      }
		      break;
		    case 7:
		      if ((velSites[x + LatticeSite::e[3][0]][y + LatticeSite::e[3][1]].isFluid()
			   && velSites[x + LatticeSite::e[4][0]][y + LatticeSite::e[4][1]].isFluid())
			  || (velSites[x + LatticeSite::e[3][0]][y + LatticeSite::e[3][1]].isSolid()
			      && velSites[x + LatticeSite::e[4][0]][y + LatticeSite::e[4][1]].isSolid())) {
			isSet = true;
			velSites[x][y].setNormalLink(k);
			velSites_[x][y].setNormalLink(k);
		      }
		      break;
		    case 8:
		      if ((velSites[x + LatticeSite::e[1][0]][y + LatticeSite::e[1][1]].isFluid()
			   && velSites[x + LatticeSite::e[4][0]][y + LatticeSite::e[4][1]].isFluid())
			  || (velSites[x + LatticeSite::e[1][0]][y + LatticeSite::e[1][1]].isSolid()
			      && velSites[x+ LatticeSite::e[4][0]][y + LatticeSite::e[4][1]].isSolid())) {
			isSet = true;
			velSites[x][y].setNormalLink(k);
			velSites_[x][y].setNormalLink(k);
		      }
		      break;

		    default:
		      cout << "error in tagging fluid nodes"
			   << endl;
		    }
		  }
		}
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

double LatticeBoltzmann::InitialCondition_X(int x, int y)
{
  double a = 0.0;
  double delta = 2.*h;
  double arg; double k = (2.*M_PI)/(dims[0]-1);
  arg = (y-h*sin(k*x))/delta;

  a = (h/delta)*sin(k*x)*exp(-arg);

  return u0*(1+a);
}

double LatticeBoltzmann::InitialCondition_Y(int x, int y)
{
  double a, b, c;
  double delta = 2.*h;
  double arg; double k = (2.*M_PI)/(dims[0]-1);
  arg = (y-h*sin(k*x))/delta;

  a = u0*h*k*cos(k*x);
  b = exp(-arg);
  c = 1.+(h/delta)*sin(k*x);

  return a*b*c;

}

void LatticeBoltzmann::setSpgeLayer(int nbOfSpgeNodes)
{
  spgeFirstNode = (dims[1]-1) - (nbOfSpgeNodes-1);
  ySpge = dims[1] - floor((1./4)*nbOfSpgeNodes);
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
