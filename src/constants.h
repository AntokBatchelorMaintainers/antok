#ifndef ANTOK_CONSTANTS_H
#define ANTOK_CONSTANTS_H

#include<assert.h>

namespace antok {

	class Initializer;

	class Constants {

		friend class antok::Initializer;

	  public:

		static const double& chargedKaonMass() { assert(initialized()); return _charged_kaon_mass; };
		static const double& chargedPionMass() { assert(initialized()); return _charged_pion_mass; };
		static const double& protonMass() { assert(initialized()); return _proton_mass; };

		static double RPDPeriod();
		static double RPDRingARadius() { return 12.0; /* cm, from 2008 spectro paper */ }
		static double RPDRingBRadius() { return 75.0; /* cm, from 2008 spectro paper */ }
		static double TargetRadius() { return 1.75; /* cm, from 2008 spectro paper */ }
		static double TargetUpstreamEdge() { return -68.5; /* cm */}
		static double TargetDownstreamEdge() { return -28.5; /* cm */}

		static const unsigned int nParticles() { assert(initialized()); return _n_particles; };

		static const bool& initialized() { return _initialized; };

	  private:

		static double _charged_kaon_mass;
		static double _charged_pion_mass;
		static double _proton_mass;

		static unsigned int _n_particles;

		static bool _initialized;

		static bool set_charged_kaon_mass(const double& charged_kaon_mass);
		static bool set_charged_pion_mass(const double& charged_pion_mass);
		static bool set_proton_mass(const double& proton_mass);

		static bool set_n_particles(const unsigned int& n_particles);

	};

}

#endif

