#ifndef HLIB_PLOTTER_H
#define HLIB_PLOTTER_H

#include<vector>

#include<TDirectory.h>
#include<TH1D.h>
#include<TH2D.h>

#include<event.h>
#include<plot.h>

namespace hlib {

	class Plotter {

	  public:
		Plotter();

		void fill(const hlib::Event& event, int cutmask);

		void save(TDirectory* dir);

		double XMass;
		double XMom;
		double CalcBeamE;
		double RPDMult;
		double PrimVX;
		double PrimVY;
		double PrimVZ;
		double ProtonMass;
		double TPrim;
		double RPDDeltaPhi;
		double TrigMask;

	  private:

		std::vector<hlib::Plot> _plots;

	};

}

#endif
