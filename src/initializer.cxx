
#include<yaml-cpp/yaml.h>

#include<constants.h>
#include<cutter.h>
#include<event.h>
#include<initializer.h>
#include<plotter.h>

antok::Initializer* antok::Initializer::_initializer = 0;

antok::Initializer* antok::Initializer::instance() {

	if(_initializer == 0) {
		_initializer = new antok::Initializer();
	};
	return _initializer;

};

antok::Initializer::Initializer()
	: _cutter(0),
	  _event(0),
	  _plotter(0)
{

};

bool antok::Initializer::readConfigFile(const std::string& filename) {

	bool error = false;

	// Load the config file
	YAML::Node config;
	try {
		config = YAML::LoadFile(filename);
	} catch (YAML::ParserException e) {
		std::cerr<<"Error parsing config file: "<<e.what()<<"."<<std::endl;
		error = true;
	}

	// Number of particles
	try {
		if(not (config["NumberOfParticles"] and antok::Constants::set_n_particles(config["NumberOfParticles"].as<unsigned int>()))) {
			std::cerr<<"Could not set number of particles."<<std::endl;
			error = true;
		}
	} catch (YAML::TypedBadConversion<int> e) {
		std::cerr<<"Conversion error when reading \"NumberOfParticles\": "<<e.what()<<"."<<std::endl;
		error = true;
	}

	// Get the constants
	YAML::Node constants = config["Constants"];
	try {
		if(not (constants["ChargedKaonMass"] and antok::Constants::set_charged_kaon_mass(constants["ChargedKaonMass"].as<double>()))) {
			std::cerr<<"Could not set charged kaon mass."<<std::endl;
			error = true;
		}
		if(not (constants["ChargedPionMass"] and antok::Constants::set_charged_pion_mass(constants["ChargedPionMass"].as<double>()))) {
			std::cerr<<"Could not set pion mass."<<std::endl;
			error = true;
		}
		if(not (constants["ProtonMass"] and antok::Constants::set_proton_mass(constants["ProtonMass"].as<double>()))) {
			std::cerr<<"Coud not set proton mass."<<std::endl;
			error = true;
		}
	} catch (YAML::TypedBadConversion<double> e) {
		std::cerr<<"Conversion error when reading constants: "<<e.what()<<"."<<std::endl;
		error = true;
	}

	antok::Constants::_initialized = (not error);

	return (not error);

};

antok::Cutter& antok::Initializer::get_cutter() {

	if(_cutter == 0) {
		_cutter = antok::Cutter::instance();
	}
	return *_cutter;

};

antok::Event& antok::Initializer::get_event() {

	if(_event == 0) {
		_event = antok::Event::instance();
	}
	return *_event;

};

antok::Plotter& antok::Initializer::get_plotter() {

	if(_plotter == 0) {
		_plotter = antok::Plotter::instance();
	}
	return *_plotter;

};

