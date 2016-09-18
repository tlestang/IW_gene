#include <iostream>
#include "LatticeSite.h"

using namespace std;

const int ThermalSite::c[4][2] = {{1,0}, {0,1}, {-1,0}, {0, -1}};
const int VelSite::e[9][2] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
const double VelSite::w[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};


LatticeSite::LatticeSite()
{

}
void VelSite::init(SiteType t, double rho, double u[2],
		   double coef_force_)
{
	setType(t);

	coef_force = coef_force_;
	
	for (int k=0; k<q; k++)
		f[k] = fEq(k, rho, u);
}

void ThermalSite::init(SiteType t, double T, double u[2],
		       double coef_force_)
{
	setType(t);
	
	for (int k=0; k<q; k++)
		f[k] = fEq(k, T, u);
}

void LatticeSite::setType(SiteType t)
{
	type = t;
}

bool LatticeSite::isSolid()
{
  return type == Solid;
}
bool LatticeSite::isFluid()
{
  return type == Fluid;
}

VelSite::VelSite()
{
  q = 9;
  T0 = 0.5;
}
//ThermalSite::ThermalSite(double T0_) : T0(T0_)
ThermalSite::ThermalSite()
{
  q = 4; 
}


void VelSite::computeRhoAndU(double& rho, double u[2])
{
  rho = 0;

  for (int k=0; k<q; k++)
    {rho += f[k];}

  u[0] = 0;
  u[1] = 0;

  for (int k=0; k<9; k++)
    {
      u[0] += f[k]*e[k][0];
      u[1] += f[k]*e[k][1];
    }
  u[0] /= rho;
  u[1] /= rho;
}

void ThermalSite::computeRhoAndU(double & rho, double u[2])
{
}
void ThermalSite::computeRhoAndU(double &T)
{
  T = 0.0;
  for (int k=0; k<q; k++)
    T += f[k];
}

double VelSite::fEq(int k, double rho, double u[2])
{
  double u2 = u[0]*u[0] + u[1]*u[1];
  double eu = e[k][0]*u[0] + e[k][1]*u[1];

  return rho * w[k] * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*u2);
}
double ThermalSite::fEq(int k, double T, double u[2])
{
  double eu = c[k][0]*u[0] + c[k][1]*u[1];
  
  return (1./4.)*T*(1.0 + 2.*eu);
}


void VelSite::collide(double& rho, double& T, double u[2], double omega)
{
  //computeRhoAndU(rho, u);
  double eu, a, F;
  a = (1.-0.5*omega);
  for (int k=0; k<q; k++)
    {
      //Compute the force from Boussinesq approx
      eu = e[k][0]*u[0] + e[k][1]*u[1];
      F = rho*coef_force*T;
      //force = a * w[k] * (3.*(e[k][1]*F-u[1]*F) + 9.*eu*F);
      force =  3.*w[k]*rho*coef_force*T*e[k][1];
      f[k]= f[k]*(1.-omega) +fEq(k,rho,u)*omega + force;	  
      // f[k] *= (1.0-omega);
      // f[k] += omega*fEq(k, rho, u);
    }
  
}
void ThermalSite::collide(double& rho, double& T, double u[2], double omega)
{
  //computeRhoAndU(T);

		for (int k=0; k<q; k++)
		{
			f[k] *= (1.0-omega);
			f[k] += omega*fEq(k, T, u);
		}
}




