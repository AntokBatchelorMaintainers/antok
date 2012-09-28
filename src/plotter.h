#ifndef HLIB_PLOTTER_H
#define HLIB_PLOTTER_H

#include<vector>

#include<TDirectory.h>
#include<TH1D.h>
#include<TH2D.h>

#include<event.h>
#include<plot.h>

namespace antok {

	class Plotter {

	  public:

		static Plotter* instance();

		void fill(const antok::Event& event, int cutmask);

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
		double TrigMask;
		double beam_time;
		double cedarTheta_X;
		double cedarTheta_Y;

		double RPDDeltaPhi;
/*		double RPDPhiRes;
		double RPDDeltaPhi_fhaas;
		double RPDPhiRes_fhaas;
		double RPDDeltaPhiAbs;
		double RPDDeltaPhiAbs_fhaas;
*/
	  private:

		Plotter();

		static Plotter* _plotter;

		std::vector<antok::Plot> _plots;

	};

}

#endif
