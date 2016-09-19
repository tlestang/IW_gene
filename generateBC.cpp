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
  double unull[2] = {0.0, 0.0};
  
    for (int xx=1;xx<dims[0]-1;xx++)
    {
      a = 2.*M_PI / (dims[0]-2);
      yy = floor(h*sin(a*(xx-1))) + h + 1; // +1 because there are ghost nodes underneath
      
      for(int y=1;y<yy;y++)
	{
	  velSites[xx][y].init(LatticeSite::Solid, 1.0, unull, 0.0);
	  velSites_[xx][y].init(LatticeSite::Solid, 1.0, unull, 0.0);
	}
            
    }
}

