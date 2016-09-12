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

  string folderName = "test/";
  string instru = "mkdir " + folderName;
  system(instru.c_str());
  instru = "mkdir " + folderName + "vtk_fluid/";
  system(instru.c_str());
  
  
  double **rho = NULL;;
  double ***velocity = NULL;
  double **temp = NULL;
  
  double gr = 0.001; double beta = 1.;
  double nu, kappa, h, period;
  int Dx, Dy;

  ifstream input_file(argv[1]);
  if(input_file.is_open())
    {
      input_file >> Dx;
      input_file >> Dy;
      input_file >> h;
      input_file >> period;
      input_file >> kappa;
      input_file >> nu;
      
      input_file.close();
    }
  else
    {
      cout << "ERROR - Could not open input file " << endl; return(0);
    }
  
  int dims[2] = {Dx, Dy};
  double omega[2];
  omega[0] = 1./(3.*nu+0.5);
  omega[1]  = 1./(2.*kappa+0.5);


  LatticeBoltzmann *lb;

  lb = new LatticeBoltzmann(dims, omega, gr*beta);

  lb->getDensityAndVelocityField(temp, rho, velocity);

  int N = 100;
  int k = 0;
  for (int i=0;i<N;i++)
    {
      lb->update();
      // if(i%1000 == 0)
      // 	{
      // 	  write_fluid_vtk(tt, dims[0], dims[1], temp, velocity, folderName.c_str());
      // 	  tt++;
      // 	}
    }    
  
}
