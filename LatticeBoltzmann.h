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
		
		double **p;
		double **T;
		double ***u;

		VelSite **velSites;
		VelSite **velSites_;
		
		ThermalSite **thermalSites;
		ThermalSite **thermalSites_;

		Walls *w;
		HotWall *hw;
		ColdWall *cw;
		/* TopWall *tw; */
		/* BottomWall *bw; */

		void streamToNeighbors(int x, int y);

	public:
		LatticeBoltzmann(const int*, const double*, double);
		~LatticeBoltzmann();

		void generateGeometry();

		void setSite(int x, int y, LatticeSite::SiteType type, double u[2]);
		void getDensityAndVelocityField(double **&, double **&, double ***&);
		void update();
};

#endif