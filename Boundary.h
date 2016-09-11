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

class Walls : public Boundary{
 private:
 public:
  Walls(const int d[2]);
  virtual void BoundaryCondition();
  virtual void BoundaryCondition(VelSite **sites, VelSite **_sites);

};

class HotWall : public Boundary{
 private:
  
 public:
  HotWall(const int d[2]);
  virtual void BoundaryCondition();
  virtual void BoundaryCondition(ThermalSite**);
  //virtual void BoundaryCondition(VelSite**, ThermalSite**, double **T, double ***u);
};

class ColdWall : public Boundary{
 private:
  
 public:
  ColdWall(const int d[2]);
  virtual void BoundaryCondition();
  virtual void BoundaryCondition(ThermalSite**);
  //virtual void BoundaryCondition(VelSite**, ThermalSite**, double **T, double ***u);
  };
class TopWall : public Boundary{
 private:
  
 public:
  TopWall(const int d[2]);
  virtual void BoundaryCondition();
  virtual void BoundaryCondition(VelSite**, ThermalSite**, double **T, double ***u);
  };
class BottomWall : public Boundary{
 private:
  
 public:
  BottomWall(const int d[2]);
  virtual void BoundaryCondition();
  virtual void BoundaryCondition(VelSite**, ThermalSite**, double **T, double ***u);
  };




#endif
