
#ifndef HLIB_DATA_H
#define HLIB_DATA_H

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

		double Mom_x1;
		double Mom_x2;
		double Mom_x3;
		double Mom_x4;
		double Mom_x5;

		double Mom_y1;
		double Mom_y2;
		double Mom_y3;
		double Mom_y4;
		double Mom_y5;

		double Mom_z1;
		double Mom_z2;
		double Mom_z3;
		double Mom_z4;
		double Mom_z5;

		double chi2PV;

		double theta_RICH_1;
		double theta_RICH_2;
		double theta_RICH_3;
		double theta_RICH_4;
		double theta_RICH_5;

		int PID_RICH_1;
		int PID_RICH_2;
		int PID_RICH_3;
		int PID_RICH_4;
		int PID_RICH_5;

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

		double zmax1;
		double zmax2;
		double zmax3;
		double zmax4;
		double zmax5;

	};

}

#endif
