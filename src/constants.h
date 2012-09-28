#ifndef ANTOK_CONSTANTS_H
#define ANTOK_CONSTANTS_H

#include<assert.h>

namespace antok {

	class Initializer;

	class Constants {

	  public:

		static const double& charged_kaon_mass() { assert(initialized()); return _charged_kaon_mass; };
		static const double& charged_pion_mass() { assert(initialized()); return _charged_pion_mass; };
		static const double& proton_mass() { assert(initialized()); return _proton_mass; };

		static const unsigned int& n_particles() { assert(initialized()); return _n_particles; };

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

		friend class antok::Initializer;

	};
/*
//	const double PION_MASS = 0.13957018;
//	const double PROTON_MASS = 0.938272046;
	const double PION_MASS = 0.139567;
	const double PROTON_MASS = 0.9382723;

	const double CHARGED_KAON_MASS = 0.493677;

	const unsigned int N_PARTICLES = 5;
*/
}

#endif
