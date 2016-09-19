#ifndef LATTICE_BOLTZMANN_H
#define LATTICE_BOLTZMANN_H

#include "LatticeSite.h"
#include "Boundary.h"

class LatticeBoltzmann
{
	private:
		int dims[2];
		double omega[2];
		double coef_force;
		double u0, N2, h;
		double spgeFirstNode, ySpge;
		
		double **rho;
		double **T;
		double ***u;

		VelSite **velSites;
		VelSite **velSites_;
		
		ThermalSite **thermalSites;
		ThermalSite **thermalSites_;

		TopWall *w;
		Topography *topo;

		/* TopWall *tw; */
		/* BottomWall *bw; */

		void streamToNeighbors(int x, int y);

	public:
		LatticeBoltzmann(const int*, const double*, double, double, int, double);
		~LatticeBoltzmann();

		void generateGeometry();
		void TagFluidNodes();
		double InitialCondition_X(int, int);
		double InitialCondition_Y(int, int);

		void setSite(int x, int y, LatticeSite::SiteType type, double u[2]);
		void setSpgeLayer(int);
		void getDensityAndVelocityField(double **&, double **&, double ***&);
		void update();
		
};

#endif
