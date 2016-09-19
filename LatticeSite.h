#ifndef LATTICE_SITE_H
#define LATTICE_SITE_H

class LatticeSite
{
	public:
		enum SiteType {Fluid = 'f', Solid = 'b', Null = 'n'};
		enum FluidTag {FluidSolid = 'b', InnerFluid = 'i', DefaultTag='y'};

		double f[9]; int brLinks[9];

	protected:
		int q;
		double force; 
		double coef_force;
		SiteType type; FluidTag tag;

	public:

		LatticeSite();
		virtual double fEq(int k, double rho, double u[2]) = 0;
		virtual void computeRhoAndU(double& rho, double u[2]) = 0;
		virtual void collide(double& rho, double& T, double u[2], double omega) = 0;
		virtual void init(SiteType, double, double*, double) = 0;
		bool isSolid(); bool isFluid();
		bool isFluidSolid(); bool isInnerFluid();
		void setType(SiteType t);
		void setTag(FluidTag tag);
		
};

class VelSite : public LatticeSite
{
 public:
  static const int e[9][2];
 private:
  static const double w[9];
  double T0, g, beta;
 public:
  VelSite();
  virtual void init(SiteType, double, double*, double);
  virtual double fEq(int k, double rho, double u[2]);
  virtual void computeRhoAndU(double& rho, double u[2]);
  virtual void collide(double& rho, double& T, double u[2], double omega);
};

class ThermalSite : public LatticeSite
{
 public:
  static const int c[4][2];

  ThermalSite();
  virtual void init(SiteType, double, double*, double);
  virtual double fEq(int k, double rho, double u[2]);
  virtual void computeRhoAndU(double&, double u[2]);
  virtual void computeRhoAndU(double& T);
  virtual void collide(double &rho, double &T, double u[2], double omega);
};
#endif
