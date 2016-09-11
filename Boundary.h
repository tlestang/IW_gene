#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "LatticeSite.h"


class Boundary{
protected:
  int nbNodes;
  int **nodes;
public:
  Boundary();
  virtual void BoundaryCondition() = 0;
  
};

class TopWall : public Boundary{
 private:
 public:
  Walls(const int d[2]);
  virtual void BoundaryCondition();
  virtual void BoundaryCondition(VelSite **sites, VelSite **_sites);

};

class Topography : public Walls{
 public:
  Topography(const int d[2], int h, int period, VelSite**, VelSite**);
  void FreeSlipBC(VelSite **sites, VelSite **_sites);
  void TemperatureBC(VelSite **velSites, ThermalSite **thermalSites,
				   double **T, double ***u)
};

#endif
