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

int main()
{

  string folderName = "test/";
  string instru = "mkdir " + folderName;
  system(instru.c_str());
  instru = "mkdir " + folderName + "vtk_fluid/";
  system(instru.c_str());
  
  
  double **p = NULL;;
  double ***velocity = NULL;
  double **temp = NULL;
  
  int dims[2] = {129, 65};

  double dx = 1./(dims[0]-1); double gr = 0.001; double beta = 1.;
  double dt = sqrt(gr*dx);
  double Pr = 0.71; double Ra = 20000.;

  double nu = sqrt(Pr/Ra)*dt/(dx*dx);
  double kappa  = sqrt(1./(Pr*Ra))*dt/(dx*dx);

  double omega[2];
  omega[0] = 1./(3.*nu+0.5);
  omega[1]  = 1./(2.*kappa+0.5);

  cout << "tau = " << 1./omega[0] << " tau_prim = " << 1./omega[1] << endl;
  int tt = 0; double T0;

  double Tref = 1.0;
  T0 = sqrt((dims[0]-1)/(beta*gr*Tref));

  double T = 1000;
  int N  = floor(T*T0);
  
  LatticeBoltzmann *lb;

  lb = new LatticeBoltzmann(dims, omega, gr*beta);

  lb->getDensityAndVelocityField(temp, p, velocity);

  int k = 0;
  for (int i=0;i<N;i++)
    {
      lb->update();
            if(i%(N/100)==0)
      	{
      	  cout << k << "%" << endl;
      	  k++;
      	}
      if(i%10 == 0)
      	{
      	  write_fluid_vtk(tt, dims[0], dims[1], temp, velocity, folderName.c_str());
      	  tt++;
      	}
    }    
  
}
