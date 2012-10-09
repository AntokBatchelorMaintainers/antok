
#include<sstream>

#include<TFile.h>
#include<TTree.h>

#include<constants.h>
#include<cutter.h>
#include<data.hpp>
#include<event.h>
#include<functions.hpp>
#include<initializer.h>
#include<object_manager.h>
#include<plotter.h>

antok::Initializer* antok::Initializer::_initializer = 0;

antok::Initializer* antok::Initializer::instance() {

	if(_initializer == 0) {
		_initializer = new antok::Initializer();
	};
	return _initializer;

};

antok::Initializer::Initializer()
	: _config(0)
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

bool antok::Initializer::initializeCutter() {

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();

	if(objectManager->_cutter == 0) {
		objectManager->_cutter = antok::Cutter::instance();
	}
	return true;

};

bool antok::Initializer::initializeData() {

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();

	if(objectManager->_data != 0) {
		std::cerr<<"Data seems to be initialized already."<<std::endl;
		return false;
	}
	objectManager->_data = new antok::Data();
	antok::Data& data = objectManager->getData();

	if(_config == 0) {
		std::cerr<<"Cannot create data object without reading a config file first."<<std::endl;
		return false;
	};
	YAML::Node& config = *_config;

	if(not config["TreeBranches"]) {
		std::cerr<<"TreeBranches not found in configuration file."<<std::endl;
		return false;
	};
	if(not config["TreeBranches"]["onePerEvent"] or not config["TreeBranches"]["onePerParticle"]) {
		std::cerr<<"TreeBranches[{\"onePerEvent\"|\"onePerParticle\"}] not found in configuration file."<<std::endl;
		return false;
	}

	// Get all the branches in the tree and fill the data maps
	YAML::Node perEventTreeBranches= config["TreeBranches"]["onePerEvent"];
	YAML::Node perParticleTreeBranches= config["TreeBranches"]["onePerParticle"];
	for(YAML::const_iterator typeIt = perEventTreeBranches.begin(); typeIt != perEventTreeBranches.end(); ++typeIt) {
		for(YAML::const_iterator valIt = typeIt->second.begin(); valIt != typeIt->second.end(); ++valIt) {
			std::string type = typeIt->first.as<std::string>();
			if(type == "double") {
				if(not data.insertDouble(valIt->as<std::string>())) {
					std::cerr<<"Could not insert double \""<<valIt->as<std::string>()<<"\" (double entry?)."<<std::endl;
					return false;
				}
			} else if(type == "int") {
				if(not data.insertInt(valIt->as<std::string>())) {
					std::cerr<<"Could not insert int \""<<valIt->as<std::string>()<<"\" (double entry?)."<<std::endl;
					return false;
				}
			} else if(type == "Long64_t") {
				if(not data.insertLong64_t(valIt->as<std::string>())) {
					std::cerr<<"Could not insert Long64_t \""<<valIt->as<std::string>()<<"\" (double entry?)."<<std::endl;
					return false;
				}
			} else {
				std::cerr<<"Data type \""<<type<<"\" not supported."<<std::endl;
				return false;
			}
		}
	}
	const unsigned int& N_PARTICLES = antok::Constants::n_particles();
	for(YAML::const_iterator typeIt = perParticleTreeBranches.begin(); typeIt != perParticleTreeBranches.end(); ++typeIt) {
		for(YAML::const_iterator valIt = typeIt->second.begin(); valIt != typeIt->second.end(); ++valIt) {
			std::string type = typeIt->first.as<std::string>();
			if(type == "double") {
				for(unsigned int i = 0; i < N_PARTICLES; ++i) {
					std::stringstream strStr;
					strStr<<(valIt->as<std::string>())<<(i+1);
					if(not data.insertDouble(strStr.str())) {
						std::cerr<<"Could not insert double \""<<strStr.str()<<"\" (double entry?)."<<std::endl;
						return false;
					}
				}
			} else if(type == "int") {
				for(unsigned int i = 0; i < N_PARTICLES; ++i) {
					std::stringstream strStr;
					strStr<<(valIt->as<std::string>())<<(i+1);
					if(not data.insertInt(strStr.str())) {
						std::cerr<<"Could not insert int \""<<strStr.str()<<"\" (double entry?)."<<std::endl;
						return false;
					}
				}
			} else {
				std::cerr<<"Data type \""<<type<<"\" not supported."<<std::endl;
				return false;
			}
		}
	}

	// Set the branch addresses of the tree
	if(not config["TreeName"]) {
		std::cerr<<"\"TreeName\" not found in configuration file."<<std::endl;
		return false;
	}
	TFile* inFile = objectManager->getInFile();
	TTree* inTree = dynamic_cast<TTree*>(inFile->Get(config["TreeName"].as<std::string>().c_str()));
	if(inTree == 0) {
		std::cerr<<"Could not open input TTree."<<std::endl;
		return false;
	}
	objectManager->_inTree = inTree;

	for(std::map<std::string, double>::iterator it = data.doubles.begin(); it != data.doubles.end(); ++it) {
		inTree->SetBranchAddress(it->first.c_str(), &(it->second));
	}
	for(std::map<std::string, int>::iterator it = data.ints.begin(); it != data.ints.end(); ++it) {
		inTree->SetBranchAddress(it->first.c_str(), &(it->second));
	}
	for(std::map<std::string, Long64_t>::iterator it = data.long64_ts.begin(); it != data.long64_ts.end(); ++it) {
		inTree->SetBranchAddress(it->first.c_str(), &(it->second));
	}

	return true;

}

bool antok::Initializer::initializeEvent() {

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();

	if(objectManager->_event != 0) {
		std::cerr<<"Event seems to be initialized already."<<std::endl;
		return false;
	}
	objectManager->_event = antok::Event::instance();
	YAML::Node& config = *_config;

	if(not config["CalculatedQuantities"]) {
		std::cerr<<"Warning: \"CalculatedQuantities\" not found in configuration file."<<std::endl;
	}
	for(YAML::const_iterator calcQuantity_it = config["CalculatedQuantities"].begin(); calcQuantity_it != config["CalculatedQuantities"].end(); calcQuantity_it++) {

		YAML::Node calcQuantity = (*calcQuantity_it);
		std::vector<std::string> quantityBaseNames;
		if(calcQuantity["Name"].IsSequence()) {
			quantityBaseNames = calcQuantity["Name"].as<std::vector<std::string> >();
		} else {
			quantityBaseNames.push_back(calcQuantity["Name"].as<std::string>());
		}
		if(quantityBaseNames.size() < 1) {
			std::cerr<<"Did not find a name to save calculated quantity to."<<std::endl;
			return false;
		}

		if(not (calcQuantity["Function"] and calcQuantity["Function"]["Name"])) {
			std::cerr<<"No Function or no function name for calculated quantity \""<<quantityBaseNames[0]<<"\"."<<std::endl;
			return false;
		}

		std::vector<int> indices;
		if(calcQuantity["Indices"]) {
			indices = calcQuantity["Indices"].as<std::vector<int> >();
		} else {
			indices.push_back(-1);
		}

		const YAML::Node& function = calcQuantity["Function"];
		std::string functionName = function["Name"].as<std::string>();

		for(unsigned int indices_i = 0; indices_i < indices.size(); ++indices_i) {

			std::vector<std::string> quantityNames;
			for(unsigned int baseName_i = 0; baseName_i < quantityBaseNames.size(); ++baseName_i) {
				if(indices[indices_i] > 0) {
					std::stringstream strStr;
					strStr<<quantityBaseNames[baseName_i]<<indices[indices_i];
					quantityNames.push_back(strStr.str());
				} else {
					quantityNames.push_back(quantityBaseNames[baseName_i]);
				}
			}

			if(functionName != "getTs" and quantityNames.size() > 1) {
				std::cerr<<"Too many names for function \""<<functionName<<"\"."<<std::endl;
				return false;
			}

			if(functionName == "getLorentzVec") {
				if(not registerGetLorentzVec(function, quantityNames[0], indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "mass") {
				if(not registerMass(function, quantityNames[0], indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "sum") {
				if(not registerSum(function, quantityNames[0], indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "getBeamLorentzVector") {
				if(not registerGetBeamLorentzVector(function, quantityNames[0], indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "getTs") {
				if(not registerGetTs(function, quantityNames, indices[indices_i])) {
					return false;
				}
				continue;
			} else {
				std::cerr<<"Function type \""<<functionName<<"\" not supported."<<std::endl;
				return false;
			}

		} // End of loop over all indices

	} // End of loop over all CalculatedQuantities

	return true;

};

bool antok::Initializer::initializePlotter() {

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();

	if(objectManager->_plotter == 0) {
		objectManager->_plotter = antok::Plotter::instance();
	}
	return true;

};

bool antok::Initializer::registerGetBeamLorentzVector(const YAML::Node& function, std::string& quantityName, int index)
{

	antok::Data& data = antok::ObjectManager::instance()->getData();
	antok::Event& event = antok::ObjectManager::instance()->getEvent();

	if(not (function["dX"] and function["dY"] and function["XLorentzVec"])) {
		std::cerr<<"Argument missing in configuration file for variable \""<<quantityName<<"\"."<<std::endl;
		return false;
	}

	std::string dXArg = function["dX"].as<std::string>();
	std::string dYArg = function["dY"].as<std::string>();
	std::string xLorentzVecArg = function["XLorentzVec"].as<std::string>();

	if(index > 0) {
		std::stringstream strStr;
		strStr<<dXArg<<index;
		dXArg = strStr.str();
		strStr.str("");
		strStr<<dYArg<<index;
		dYArg = strStr.str();
		strStr.str("");
		strStr<<xLorentzVecArg<<index;
		xLorentzVecArg = strStr.str();
	}

	double* dXaddr = data.getDoubleAddr(dXArg);
	double* dYaddr = data.getDoubleAddr(dYArg);
	TLorentzVector* xLorentzVecAddr = data.getLorentzVectorAddr(xLorentzVecArg);
	if((dXaddr == 0) or (dYaddr == 0) or (xLorentzVecAddr == 0)) {
		std::cerr<<"Arguments not found in Data's doubles map for variable entry  \""<<quantityName<<"\"."<<std::endl;
		return false;
	}

	if(not data.insertLorentzVector(quantityName)) {
		std::cerr<<"Could not insert TLorentzVector when creating variable entry \""<<quantityName<<"\"."<<std::endl;
		return false;
	}

	event._functions.push_back(new antok::GetBeamLorentzVec(dXaddr, dYaddr, xLorentzVecAddr, data.getLorentzVectorAddr(quantityName)));

	return true;

};

bool antok::Initializer::registerGetLorentzVec(const YAML::Node& function, std::string& quantityName, int index)
{

	antok::Data& data = antok::ObjectManager::instance()->getData();
	antok::Event& event = antok::ObjectManager::instance()->getEvent();

	bool pType;
	if(function["X"] and function["Y"] and function["Z"] and function["M"]) {
		pType = false;
	} else if (function["Px"] and function["Py"] and function["Pz"] and function["E"]) {
		pType = true;
	} else {
		std::cerr<<"Argument missing in configuration file for variable \""<<quantityName<<"\"."<<std::endl;
		return false;
	}
	std::string xArg;
    std::string yArg;
    std::string zArg;
	std::string mArg;
	double* xAddr;
	double* yAddr;
	double* zAddr;
	double* mAddr;
	if(pType) {
		xArg = function["Px"].as<std::string>();
		yArg = function["Py"].as<std::string>();
		zArg = function["Pz"].as<std::string>();
		mArg = function["E"].as<std::string>();
	} else {
		xArg = function["X"].as<std::string>();
		yArg = function["Y"].as<std::string>();
		zArg = function["Z"].as<std::string>();
		try {
			function["M"].as<double>();
		} catch(YAML::TypedBadConversion<double> e) {
			std::cerr<<"Mass in function \"mass\" should be of type double."<<std::endl;
			return false;
		}
		mAddr = new double();
		(*mAddr) = function["M"].as<double>();
	}
	if(index > 0) {
		std::stringstream strStr;
		strStr<<xArg<<index;
		xArg = strStr.str();
		strStr.str("");
		strStr<<yArg<<index;
		yArg = strStr.str();
		strStr.str("");
		strStr<<zArg<<index;
		zArg = strStr.str();
		if(pType) {
			strStr.str("");
			strStr<<mArg<<index;
			mArg = strStr.str();
		}
	}
	xAddr = data.getDoubleAddr(xArg);
	yAddr = data.getDoubleAddr(yArg);
	zAddr = data.getDoubleAddr(zArg);
	if(pType) {
		mAddr = data.getDoubleAddr(mArg);
	}
	if((xAddr == 0) or (yAddr == 0) or (zAddr == 0)) {
		std::cerr<<"Arguments not found in Data's doubles map for function \"mass\"."<<std::endl;
		return false;
	}

	if(not data.insertLorentzVector(quantityName)) {
		std::cerr<<"Could not insert TLorentzVector when creating variable entry \""<<quantityName<<"\"."<<std::endl;
		return false;
	}

	event._functions.push_back(new antok::GetLorentzVec(xAddr, yAddr, zAddr, mAddr, data.getLorentzVectorAddr(quantityName), pType));
	return true;

};

bool antok::Initializer::registerGetTs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	antok::Data& data = antok::ObjectManager::instance()->getData();
	antok::Event& event = antok::ObjectManager::instance()->getEvent();

	if(quantityNames.size() != 3) {
		std::cerr<<"Need 3 names for function \"getTs\""<<std::endl;
		return false;
	}

	if(not (function["BeamLorentzVec"] and function["XLorentzVec"])) {
		std::cerr<<"Argument missing in configuration file for variables \"[";
		for(unsigned int i = 0; i < quantityNames.size()-1; ++i) {
			std::cerr<<quantityNames[i]<<", ";
		}
		std::cerr<<quantityNames[quantityNames.size()-1]<<"]\"."<<std::endl;
		return false;
	}

	std::string beamLVArg = function["BeamLorentzVec"].as<std::string>();
	std::string xLVArg = function["XLorentzVec"].as<std::string>();
	if(index > 0) {
		std::stringstream strStr;
		strStr<<beamLVArg<<index;
		beamLVArg = strStr.str();
		strStr.str("");
		strStr<<xLVArg<<index;
		xLVArg = strStr.str();
	}

	TLorentzVector* beamLVAddr = data.getLorentzVectorAddr(beamLVArg);
	TLorentzVector* xLVAddr = data.getLorentzVectorAddr(xLVArg);
	if((beamLVAddr == 0) or (xLVAddr == 0)) {
		std::cerr<<"Arguments not found in Data's doubles map for function \"GetTs\"."<<std::endl;
		return false;
	}

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insertDouble(quantityNames[i])) {
			std::cerr<<"Could not insert double \""<<quantityNames[i]<<"\" in function \"getTs\""<<std::endl;
			return false;
		}
		quantityAddrs.push_back(data.getDoubleAddr(quantityNames[i]));
	}
	event._functions.push_back(new antok::GetTs(xLVAddr, beamLVAddr, quantityAddrs[0], quantityAddrs[1], quantityAddrs[2]));

	return true;

};

bool antok::Initializer::registerMass(const YAML::Node& function, std::string& quantityName, int index)
{

	antok::Data& data = antok::ObjectManager::instance()->getData();
	antok::Event& event = antok::ObjectManager::instance()->getEvent();

	if(not function["Vector"]) {
		std::cerr<<"\"Vector\" argument not found in \"mass\" function for variable \""<<quantityName<<"\"."<<std::endl;
		return false;
	}
	std::string vectorName = function["Vector"].as<std::string>();
	if(index > 0) {
		std::stringstream strStr;
		strStr<<vectorName<<index;
		vectorName = strStr.str();
	}
	TLorentzVector* vector = data.getLorentzVectorAddr(vectorName);
	if(vector == 0) {
		std::cerr<<"Did not find TLorentzVector \""<<vectorName<<"\" in Data for variable \""<<quantityName<<"\"."<<std::endl;
		return false;
	}
	if(not data.insertDouble(quantityName)) {
		std::cerr<<"Coud not insert double \""<<quantityName<<"\" (double entry?)."<<std::endl;
		return false;
	}
	event._functions.push_back(new antok::Mass(vector, data.getDoubleAddr(quantityName)));
	return true;

};

bool antok::Initializer::registerSum(const YAML::Node& function, std::string& quantityName, int index)
{

	antok::Data& data = antok::ObjectManager::instance()->getData();
	antok::Event& event = antok::ObjectManager::instance()->getEvent();

	std::vector<double*> doubleInputAddrs;
	std::vector<int*> intInputAddrs;
	std::vector<Long64_t*> long64_tInputAddrs;
	std::vector<TLorentzVector*> lorentzVectorInputAddrs;
	std::string typeName = "notInitialized";
	std::vector<std::string> summandNames;

	// Summing over one variable with indices
	if(not function["Summands"]) {
		std::cerr<<"Argument \"Summands\" is missing in configuration file for variable \""<<quantityName<<"\"."<<std::endl;
		return false;
	}

	if(function["Summands"]["Indices"] or function["Summands"]["Name"]) {
		if(not function["Summands"]["Indices"] and function["Summands"]["Name"]) {
			std::cerr<<"Either \"Summands\" or \"Name\" found in sum function, but not both (Variable: \""<<quantityName<<"\")."<<std::endl;
			return false;
		}
		if(index > 0) {
			std::cerr<<"Cannot have sum over indices for every particle (Variable: \""<<quantityName<<"\")."<<std::endl;
			return false;
		}
		std::vector<int> inner_indices = function["Summands"]["Indices"].as<std::vector<int> >();
		typeName = function["Summands"]["Name"].as<std::string>();
		std::string summandBaseName = typeName;
		std::stringstream strStr;
		strStr<<typeName<<inner_indices[0];
		typeName = data.global_map[strStr.str()];
		for(unsigned int inner_indices_i = 0; inner_indices_i < inner_indices.size(); ++inner_indices_i) {
			int inner_index = inner_indices[inner_indices_i];
			std::stringstream strStr;
			strStr<<summandBaseName<<inner_index;
			summandNames.push_back(strStr.str());
		}
	// Summing over list of variable names
	} else {
		typeName = function["Summands"].begin()->as<std::string>();
		if(index > 0) {
			std::stringstream strStr;
			strStr<<typeName<<index;
			typeName = strStr.str();
		}
		typeName = data.global_map[typeName];
		for(YAML::const_iterator summand_it = function["Summands"].begin(); summand_it != function["Summands"].end(); summand_it++) {
			std::string variableName = summand_it->as<std::string>();
			if(index > 0) {
				std::stringstream strStr;
				strStr<<variableName<<index;
				variableName = strStr.str();
			}
			summandNames.push_back(variableName);
		}
	}

	// Now do type checking and get all the addresses
	for(unsigned int summandNames_i = 0; summandNames_i < summandNames.size(); ++summandNames_i) {

			std::string variableName = summandNames[summandNames_i];

			if(data.global_map[variableName] != typeName) {
				std::cerr<<"Cannot sum terms of different type (\""<<typeName<<"\"<>\""<<data.global_map[variableName]<<"\") when calculating variable \""<<quantityName<<"\")."<<std::endl;
				return false;
			}

			if(typeName == "double") {
				double* addr = data.getDoubleAddr(variableName);
				if(addr == 0) {
					std::cerr<<"Did not find double \""<<variableName<<"\" in Data."<<std::endl;
					return false;
				}
				doubleInputAddrs.push_back(addr);
			} else if(typeName == "int") {
				int* addr = data.getIntAddr(variableName);
				if(addr == 0) {
					std::cerr<<"Did not find int \""<<variableName<<"\" in Data."<<std::endl;
					return false;
				}
				intInputAddrs.push_back(addr);
			} else if(typeName == "Long64_t") {
				Long64_t* addr = data.getLong64_tAddr(variableName);
				if(addr == 0) {
					std::cerr<<"Did not find Long64_t \""<<variableName<<"\" in Data."<<std::endl;
					return false;
				}
				long64_tInputAddrs.push_back(&data.long64_ts[variableName]);
			} else if(typeName == "TLorentzVector") {
				TLorentzVector* addr = data.getLorentzVectorAddr(variableName);
				if(addr == 0) {
					std::cerr<<"Did not find TLorentzVector \""<<variableName<<"\" in Data."<<std::endl;
					return false;
				}
				lorentzVectorInputAddrs.push_back(addr);
			} else {
				std::cerr<<"Type \""<<typeName<<"\" is not supported in function \"sum\" (Calculated quantity: \""<<quantityName<<"\")."<<std::endl;
				return false;
			}

		}

	// And produce the function
	if(typeName == "double") {
		data.insertDouble(quantityName);
		event._functions.push_back(new antok::Sum<double>(doubleInputAddrs, data.getDoubleAddr(quantityName)));
	} else if(typeName == "int") {
		data.insertInt(quantityName);
		event._functions.push_back(new antok::Sum<int>(intInputAddrs, data.getIntAddr(quantityName)));
	} else if(typeName == "Long64_t") {
		data.insertLong64_t(quantityName);
		event._functions.push_back(new antok::Sum<Long64_t>(long64_tInputAddrs, data.getLong64_tAddr(quantityName)));
	} else if(typeName == "TLorentzVector") {
		data.insertLorentzVector(quantityName);
		event._functions.push_back(new antok::Sum<TLorentzVector>(lorentzVectorInputAddrs, data.getLorentzVectorAddr(quantityName)));
	}

	return true;

};

bool antok::Initializer::registerSum2(const YAML::Node& function, std::string& quantityName, int index)
{



};

bool antok::Initializer::argumentHandler(std::vector<std::string>& argNames, int index)
{



};

