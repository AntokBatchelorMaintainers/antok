
#include<sstream>

#include<TFile.h>
#include<TTree.h>

#include<constants.h>
#include<cutter.h>
#include<data.hpp>
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
	: _config(0),
	  _cutter(0),
	  _data(0),
	  _event(0),
	  _plotter(0)
{

};

bool antok::Initializer::readConfigFile(const std::string& filename) {

	bool error = false;

	if(_config != 0) {
		std::cerr<<"Attempting to read config file twice."<<std::endl;
		return false;
	}

	_config = new YAML::Node();
	YAML::Node& config = *_config;

	// Load the config file
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

antok::Data& antok::Initializer::get_data(TFile* infile, TTree*& inTree) {

	if(_data == 0) {
		if(_config == 0) {
			std::cerr<<"Cannot create data object without reading a config file first."<<std::endl;
			throw 1;
		};
		YAML::Node& config = *_config;

		if(not config["TreeBranches"]) {
			std::cerr<<"TreeBranches not found in configuration file."<<std::endl;
			throw 1;
		};
		if(not config["TreeBranches"]["onePerEvent"] or not config["TreeBranches"]["onePerParticle"]) {
			std::cerr<<"TreeBranches[{\"onePerEvent\"|\"onePerParticle\"}] not found in configuration file."<<std::endl;
			throw 1;
		}

		// Get all the branches in the tree and fill the data maps
		YAML::Node perEventTreeBranches= config["TreeBranches"]["onePerEvent"];
		YAML::Node perParticleTreeBranches= config["TreeBranches"]["onePerParticle"];
		_data = new antok::Data();
		antok::Data& data = *_data;
		for(YAML::const_iterator typeIt = perEventTreeBranches.begin(); typeIt != perEventTreeBranches.end(); ++typeIt) {
			for(YAML::const_iterator valIt = typeIt->second.begin(); valIt != typeIt->second.end(); ++valIt) {
				std::string type = typeIt->first.as<std::string>();
				if(type == "double") {
					data.doubles[valIt->as<std::string>()] = -8888.8;
				} else if(type == "int") {
					data.ints[valIt->as<std::string>()] = -8888;
				} else if(type == "Long64_t") {
					data.long64_ts[valIt->as<std::string>()] = -8888;
				} else {
					std::cerr<<"Data type \""<<type<<"\" not supported."<<std::endl;
					throw 1;
				}
				std::cout<<typeIt->first.as<std::string>()<<": "<<valIt->as<std::string>()<<std::endl;
			}
		}
		const unsigned int& N_PARTICLES = antok::Constants::n_particles();
		for(YAML::const_iterator typeIt = perParticleTreeBranches.begin(); typeIt != perParticleTreeBranches.end(); ++typeIt) {
			for(YAML::const_iterator valIt = typeIt->second.begin(); valIt != typeIt->second.end(); ++valIt) {
				std::string type = typeIt->first.as<std::string>();
				if(type == "double") {
					std::vector<double> vec(N_PARTICLES, -8888.8);
					data.particle_doubles[valIt->as<std::string>()] = vec;
				} else if(type == "int") {
					std::vector<int> vec(N_PARTICLES, -8888);
					data.particle_ints[valIt->as<std::string>()] = vec;
				} else {
					std::cerr<<"Data type \""<<type<<"\" not supported."<<std::endl;
					throw 1;
				}
				std::cout<<typeIt->first.as<std::string>()<<": "<<valIt->as<std::string>()<<std::endl;
			}
		}

		// Set the branch addresses of the tree
		if(not config["TreeName"]) {
			std::cerr<<"\"TreeName\" not found in configuration file."<<std::endl;
			throw 1;
		}
		inTree = dynamic_cast<TTree*>(infile->Get(config["TreeName"].as<std::string>().c_str()));
		if(inTree == 0) {
			std::cerr<<"Could not open input TTree."<<std::endl;
			throw 1;
		}

		for(std::map<std::string, double>::iterator it = data.doubles.begin(); it != data.doubles.end(); ++it) {
			inTree->SetBranchAddress(it->first.c_str(), &(it->second));
			std::cout<<"double name: "<<it->first<<" value="<<it->second<<std::endl;
		}
		for(std::map<std::string, int>::iterator it = data.ints.begin(); it != data.ints.end(); ++it) {
			inTree->SetBranchAddress(it->first.c_str(), &(it->second));
			std::cout<<"int name: "<<it->first<<" value="<<it->second<<std::endl;
		}
		for(std::map<std::string, Long64_t>::iterator it = data.long64_ts.begin(); it != data.long64_ts.end(); ++it) {
			inTree->SetBranchAddress(it->first.c_str(), &(it->second));
			std::cout<<"Long64_t name: "<<it->first<<" value="<<it->second<<std::endl;
		}
		for(std::map<std::string, std::vector<double> >::iterator it = data.particle_doubles.begin(); it != data.particle_doubles.end(); ++it) {
			std::cout<<"Particle double name: "<<it->first<<" value=["<<it->second[0];
			for(unsigned int i = 0; i < it->second.size(); ++i) {
				std::stringstream strStr;
				strStr<<it->first<<(i+1);
				inTree->SetBranchAddress(strStr.str().c_str(), &(it->second[i]));
				std::cout<<", "<<it->second[i];
			}
			std::cout<<"]"<<std::endl;
		}
		for(std::map<std::string, std::vector<int> >::iterator it = data.particle_ints.begin(); it != data.particle_ints.end(); ++it) {
			std::cout<<"Particle double name: "<<it->first<<" value=["<<it->second[0];
			for(unsigned int i = 0; i < it->second.size(); ++i) {
				std::stringstream strStr;
				strStr<<it->first<<(i+1);
				inTree->SetBranchAddress(strStr.str().c_str(), &(it->second[i]));
				std::cout<<", "<<it->second[i];
			}
			std::cout<<"]"<<std::endl;
		}
	}

	return *_data;

}

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

