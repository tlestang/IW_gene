#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <fstream>

#include "LatticeBoltzmann.h"
#include "LatticeSite.h"
#include "Boundary.h"

using namespace std;

void write_fluid_vtk(int, int, int, double**, double***, const char*);

int main(int argc, char *argv[])
{

  string folderName = "test_no_topo/";
  string instru = "mkdir " + folderName;
  system(instru.c_str());
  instru = "mkdir " + folderName + "vtk_fluid/";
  system(instru.c_str());
  
  
  double **rho = NULL;;
  double ***velocity = NULL;
  double **temp = NULL;

  // declaration of parameters in PHYSICAL units to
  // be init. from input_file.
  double U0, nuPhi, kappaPhi, D, asp, hb, Pr, N2Phi;

  // declaration of parameters in LATTICE units.

  int Nx;
  double gr = 0.001; double beta = 1.;
  double nu, kappa, tau_prim, N2;
  int h;
  int Dx, Dy;

  ifstream input_file(argv[1]);
  if(input_file.is_open())
    {
      input_file >> U0;
      input_file >> Pr;
      input_file >> D;
      input_file >> asp;
      input_file >> hb;
      input_file >> Dx;
      input_file >> N2Phi;
      input_file.close();
    }

  else
    {
      cout << "ERROR - Could not open input file " << endl; return(0);
    }

  // Now compute simulation parameters from physical values.
  //Lattice velocity set to 0.1 (Low Ma limit).
  double u0 = 0.1;

  double delta_x = D/Dx;
  double delta_t = (delta_x*u0)/U0;

  //Lattice viscosity is set so that tau = 0.51 which is expected to be stable enough.
  double tau = 0.51;
  nu = (2*tau -1)/6.; nuPhi = nu*(delta_x*delta_x/delta_t);
  kappa = nu/Pr; kappaPhi = kappa*(delta_x*delta_x/delta_t);
  tau_prim = 2.*kappa + 0.5;

  h = floor(hb/delta_x);
  cout << asp << " " << Dx << " " << endl;
  Dy = asp*(Dx-1) + 1; N2 = delta_t*delta_t*N2Phi;

  ofstream paramFile("parameters.datout");
  paramFile << "PHYSICAL SETUP IN LATTICE UNITS" << endl;
  paramFile << "--------------------" << endl;
  paramFile << "  Grid is " << Dx << "x" << Dy << endl;
  paramFile << "  tau = " << tau << " | tau_prim = " << tau_prim << endl;
  paramFile << "  Dy = " << Dy << endl;
  paramFile << "  h = " << h << endl;

  paramFile << "PHYSICAL SETUP IN PHYSICAL UNITS" << endl;
  paramFile << "  nu = " << nuPhi << endl;
  paramFile << "  kappa = " << kappaPhi << endl;
  paramFile << "  dx = " << delta_x << " m" << endl;
  paramFile << "  dt = " << delta_t << " s" << endl;
  paramFile.close();

  int dims[2] = {Dx, Dy};
  double omega[2];
  omega[0] = 1./tau;
  omega[1]  = 1./tau_prim;


  LatticeBoltzmann *lb;

  lb = new LatticeBoltzmann(dims, omega, gr*beta, u0, h, N2);

  lb->getDensityAndVelocityField(temp, rho, velocity);

  int N = 60000;
  int k = 0; int tt = 0;
  for (int i=0;i<N;i++)
    {
      lb->update();
      if(i%(N/100)==0)
      	{
      	  cout << k << "%" << endl;
      	  k++;
      	}
      if(i%100 == 0)
      	{
      	  write_fluid_vtk(tt, dims[0], dims[1], temp, velocity, folderName.c_str());
      	  tt++;
      	}
    }    
  
}
