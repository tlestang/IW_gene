#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "LatticeSite.h"


class Boundary{
protected:
  int nbNodes;
  int **nodes;
  int dims[2];
public:
  Boundary();
  virtual void BoundaryCondition() = 0;
  
};

class TopWall : public Boundary{
 private:
 public:
  TopWall(const int d[2]);
  virtual void BoundaryCondition();
  virtual void BoundaryCondition(VelSite **sites, VelSite **_sites);
  void FreeSlipBC(VelSite **sites, VelSite **_sites);
  void TemperatureBC(VelSite **velSites, ThermalSite **thermalSites,
		     double **T, double ***u);

};

class BottomWall : public Boundary{
 private:
 public:
  BottomWall(const int d[2]);
  virtual void BoundaryCondition();
  virtual void BoundaryCondition(VelSite **sites, VelSite **_sites);
  void FreeSlipBC(VelSite **sites, VelSite **_sites);
  void TemperatureBC(VelSite **velSites, ThermalSite **thermalSites,
		     double **T, double ***u);

};

class Topography : public Boundary{
 private:
  int *lbl;
 public:
  Topography(const int d[2], int h, double*** u, VelSite**, VelSite**);
  virtual void BoundaryCondition();
  void FreeSlipBC(VelSite **sites, VelSite **_sites);
  void TemperatureBC(VelSite **velSites, ThermalSite **thermalSites,
		     double **T, double ***u);
};

#endif
