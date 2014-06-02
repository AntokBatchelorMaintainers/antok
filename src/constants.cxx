#include<constants.h>

#include<cmath>
#include<iostream>

#include<TMath.h>

double antok::Constants::_charged_kaon_mass = -1.;
double antok::Constants::_charged_pion_mass = -1.;
double antok::Constants::_proton_mass = -1.;

unsigned int antok::Constants::_n_particles = -1;

bool antok::Constants::_initialized = false;

bool antok::Constants::set_charged_kaon_mass(const double& charged_kaon_mass) {
	if(charged_kaon_mass > 0.) {
		_charged_kaon_mass = charged_kaon_mass;
		if(std::fabs(_charged_kaon_mass - 0.493677) > 0.001) {
			std::cerr<<"Warning: charged kaon mass seems to be set strangely: "<<_charged_kaon_mass<<"."<<std::endl;
		}
		return true;
	}
	return false;
};

double antok::Constants::RPDPeriod() {
	return TMath::Pi() / 6.;
}

bool antok::Constants::set_charged_pion_mass(const double& charged_pion_mass) {
	if(charged_pion_mass > 0.) {
		_charged_pion_mass = charged_pion_mass;
		if(std::fabs(_charged_pion_mass - 0.13957018) > 0.001) {
			std::cerr<<"Warning: charged pion mass seems to be set strangely: "<<_charged_pion_mass<<"."<<std::endl;
		}
		return true;
	}
	return false;
};

bool antok::Constants::set_proton_mass(const double& proton_mass) {
	if(proton_mass > 0.) {
		_proton_mass = proton_mass;
		if(std::fabs(_proton_mass - 0.938272046) > 0.001) {
			std::cerr<<"Warning: charged pion mass seems to be set strangely: "<<_proton_mass<<"."<<std::endl;
		}
		return true;
	}
	return false;
};

bool antok::Constants::set_n_particles(const unsigned int& n_particles) {
	if(n_particles > 0) {
		_n_particles = n_particles;
		return true;
	};
	return false;
};
