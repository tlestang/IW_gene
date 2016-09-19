#include <iostream>
#include <fstream>
#include <cmath>

#include "LatticeBoltzmann.h"
#include "LatticeSite.h"
#include "Boundary.h"

using namespace std;

void LatticeBoltzmann::generateTopography(int h)
{
  double a; int yy;
  double uu[2] = {0.0, 0.0};
  
    for (int xx=0;xx<d[0];xx++)
    {
      a = 2.*M_PI / d[0];
      yy = floor(h*sin(a*xx)) + h;
      
      for(int y=0;y<yy;y++)
	{
	  velSites[xx][y].init(LatticeSite::Solid, 1.0, uu, 0);
	  velSites_[xx][y].init(LatticeSite::Solid, 1.0, uu, 0);
	}
            
    }
}

