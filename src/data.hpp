
namespace hlib {

	struct Data {

		static int Run;
		static int TrigMask;
		static Long64_t EvNbr;
		static int SpillNbr;

		static double X_primV;
		static double Y_primV;
		static double Z_primV;

		static double gradx;
		static double grady;

		static double Mom_x1;
		static double Mom_x2;
		static double Mom_x3;
		static double Mom_x4;
		static double Mom_x5;

		static double Mom_y1;
		static double Mom_y2;
		static double Mom_y3;
		static double Mom_y4;
		static double Mom_y5;

		static double Mom_z1;
		static double Mom_z2;
		static double Mom_z3;
		static double Mom_z4;
		static double Mom_z5;

		static double chi2PV;

		static double theta_RICH_1;
		static double theta_RICH_2;
		static double theta_RICH_3;
		static double theta_RICH_4;
		static double theta_RICH_5;

		static int PID_RICH_1;
		static int PID_RICH_2;
		static int PID_RICH_3;
		static int PID_RICH_4;
		static int PID_RICH_5;

		static double RPD_Px;
		static double RPD_Py;
		static double RPD_Pz;
		static double RPD_E;
		static double RPD_Tz;
		static double RPD_z;
		static double RPD_beta;
		static double RPD_Phi;
		static double RPD_dEA;
		static double RPD_dEB;
		static int nbrRPDTracks;

		static int isKaon;

		static double zmax1;
		static double zmax2;
		static double zmax3;
		static double zmax4;
		static double zmax5;

	};

int Data::Run = -1;
int Data::TrigMask = -1;
Long64_t Data::EvNbr = -1;
int Data::SpillNbr = -1;

double Data::X_primV = -1.;
double Data::Y_primV = -1.;
double Data::Z_primV = -1.;

double Data::gradx = -1.;
double Data::grady = -1.;

double Data::Mom_x1 = -1.;
double Data::Mom_x2 = -1.;
double Data::Mom_x3 = -1.;
double Data::Mom_x4 = -1.;
double Data::Mom_x5 = -1.;

double Data::Mom_y1 = -1.;
double Data::Mom_y2 = -1.;
double Data::Mom_y3 = -1.;
double Data::Mom_y4 = -1.;
double Data::Mom_y5 = -1.;

double Data::Mom_z1 = -1.;
double Data::Mom_z2 = -1.;
double Data::Mom_z3 = -1.;
double Data::Mom_z4 = -1.;
double Data::Mom_z5 = -1.;

double Data::chi2PV = -1.;

double Data::theta_RICH_1 = -1.;
double Data::theta_RICH_2 = -1.;
double Data::theta_RICH_3 = -1.;
double Data::theta_RICH_4 = -1.;
double Data::theta_RICH_5 = -1.;

int Data::PID_RICH_1 = -1;
int Data::PID_RICH_2 = -1;
int Data::PID_RICH_3 = -1;
int Data::PID_RICH_4 = -1;
int Data::PID_RICH_5 = -1;

double Data::RPD_Px = -1.;
double Data::RPD_Py = -1.;
double Data::RPD_Pz = -1.;
double Data::RPD_E = -1.;
double Data::RPD_Tz = -1.;
double Data::RPD_z = -1.;
double Data::RPD_beta = -1.;
double Data::RPD_Phi = -1.;
double Data::RPD_dEA = -1.;
double Data::RPD_dEB = -1.;
int Data::nbrRPDTracks = -1;

int Data::isKaon = -1;

double Data::zmax1 = -1.;
double Data::zmax2 = -1.;
double Data::zmax3 = -1.;
double Data::zmax4 = -1.;
double Data::zmax5 = -1.;

}
