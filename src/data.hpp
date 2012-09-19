
#ifndef HLIB_DATA_H
#define HLIB_DATA_H

#include<constants.hpp>

namespace hlib {

	struct Data {

		int Run;
		int TrigMask;
		Long64_t EvNbr;
		int SpillNbr;

		double X_primV;
		double Y_primV;
		double Z_primV;

		double gradx;
		double grady;

		std::vector<double> Mom_x;
		std::vector<double> Mom_y;
		std::vector<double> Mom_z;

		std::vector<double> z_max;

		std::vector<double> theta_RICH;
		std::vector<int> PID_RICH;

		double chi2PV;

		double RPD_Px;
		double RPD_Py;
		double RPD_Pz;
		double RPD_E;
		double RPD_Tz;
		double RPD_z;
		double RPD_beta;
		double RPD_Phi;
		double RPD_dEA;
		double RPD_dEB;
		int nbrRPDTracks;

		int isKaon;

		Data() {
			Mom_x.resize(hlib::N_PARTICLES);
			Mom_y.resize(hlib::N_PARTICLES);
			Mom_z.resize(hlib::N_PARTICLES);
			z_max.resize(hlib::N_PARTICLES);
			theta_RICH.resize(hlib::N_PARTICLES);
			PID_RICH.resize(hlib::N_PARTICLES);
		}

	};

}

#endif
