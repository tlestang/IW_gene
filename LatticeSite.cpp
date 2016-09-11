#include <iostream>
#include "LatticeSite.h"

using namespace std;

const int ThermalSite::c[4][2] = {{1,0}, {0,1}, {-1,0}, {0, -1}};
const int VelSite::e[9][2] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
const double VelSite::w[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};


const double VelSite::lbda = 1./3.; const double VelSite::gma = 1./12.; 
const double VelSite::sigma = 5./12.;
const double VelSite::pterm[9] = {-4.*sigma, lbda, lbda, lbda, lbda, gma, gma, gma, gma};


LatticeSite::LatticeSite()
{

}
void VelSite::init(SiteType t, double p, double u[2],
		   double om, double coef_force_)
{
	setType(t);

	omega = om; omega1 = 1.0-om;
	coef_force = coef_force_;
	
	for (int k=0; k<q; k++)
		f[k] = fEq(k, p, u);
}

void ThermalSite::init(SiteType t, double T, double u[2],
		       double om, double coef_force_)
{
	setType(t);
	omega = om; omega1 = 1.0-om;
	
	for (int k=0; k<q; k++)
		f[k] = fEq(k, T, u);
}

void LatticeSite::setType(SiteType t)
{
	type = t;
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


void VelSite::computeRhoAndU(double& p, double u[2])
{
  double u2;
  p = 0;

  for (int k=0; k<q; k++)
    {p += f[k];}

  u[0] = 0;
  u[1] = 0;

  for (int k=0; k<9; k++)
    {
      u[0] += f[k]*e[k][0];
      u[1] += f[k]*e[k][1];
    }
  u2 = u[0]*u[0] + u[1]*u[1];
  
  p = (1./(4.*sigma))*(p - 1.5*w[0]*u2);
}

void ThermalSite::computeRhoAndU(double & p, double u[2])
{
}
void ThermalSite::computeRhoAndU(double &T)
{
  T = 0.0;
  for (int k=0; k<q; k++)
    T += f[k];
}

double VelSite::fEq(int k, double p, double u[2])
{
  double u2 = u[0]*u[0] + u[1]*u[1];
  double eu = e[k][0]*u[0] + e[k][1]*u[1];
  
  return p*pterm[k] + w[k] * (3.0*eu + 4.5*eu*eu - 1.5*u2);
}
double ThermalSite::fEq(int k, double T, double u[2])
{
  double eu = c[k][0]*u[0] + c[k][1]*u[1];
  
  return (1./4.)*T*(1.0 + 2.*eu);
}


void VelSite::collide(double& p, double u[2])
{
}
void VelSite::collide(double& p, double T, double u[2])
{
  computeRhoAndU(p, u);

  for (int k=0; k<q; k++)
    {
      //Compute the force from Boussinesq approx
      //force = - 3.*w[k]*coef_force*(T-T0)*e[k][1];
      f[k] = f[k]*omega1 + omega*fEq(k, p, u); //+ force;
    }
  
}
void ThermalSite::collide(double& T, double u[2])
{
		computeRhoAndU(T);

		for (int k=0; k<q; k++)
		{
			f[k] *= (1.0-omega);
			f[k] += omega*fEq(k, T, u);
		}
}




