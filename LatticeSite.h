#ifndef LATTICE_SITE_H
#define LATTICE_SITE_H

class LatticeSite
{
	public:
		enum SiteType {Fluid = 'f', Boundary = 'b', Null = 'n'};

		double f[9];

	protected:
		int q;
		double force; double omega1, omega;
		SiteType type;
		double coef_force;

	public:
		LatticeSite();
		virtual double fEq(int k, double rho, double u[2]) = 0;
		virtual void computeRhoAndU(double& rho, double u[2]) = 0;
		virtual void collide(double& rho, double u[2]) = 0;
		virtual void init(SiteType, double, double*, double, double) = 0;
		bool isBoundary();
		bool isFluid();
		void setType(SiteType t);
};

class VelSite : public LatticeSite
{
 public:
  static const int e[9][2];
 private:
  static const double w[9];
  static const double pterm[9];
  static const double lbda, gma, sigma;
  double T0, g, beta;
 public:
  VelSite();
  virtual void init(SiteType, double, double*, double, double);
  virtual double fEq(int k, double rho, double u[2]);
  virtual void computeRhoAndU(double& rho, double u[2]);
  virtual void collide(double&, double u[2]);
  virtual void collide(double& rho, double T, double u[2]);
};

class ThermalSite : public LatticeSite
{
 public:
  static const int c[4][2];

  ThermalSite();
  virtual void init(SiteType, double, double*, double, double);
  virtual double fEq(int k, double rho, double u[2]);
  virtual void computeRhoAndU(double&, double u[2]);
  virtual void computeRhoAndU(double& T);
  virtual void collide(double& T, double u[2]);
};
#endif
