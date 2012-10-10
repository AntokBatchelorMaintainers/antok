
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
			std::string type = getYAMLStringSafe(typeIt->first);
			std::string name = getYAMLStringSafe(*valIt);
			if(name == "") {
				std::cerr<<"Conversion to std::string failed for one of the \"TreeBranches\"' \"onePerEvent\" "<<type<<"s."<<std::endl;
				return false;
			}
			if(type == "double") {
				if(not data.insertDouble(name)) {
					std::cerr<<getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if(type == "int") {
				if(not data.insertInt(name)) {
					std::cerr<<getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if(type == "Long64_t") {
				if(not data.insertLong64_t(name)) {
					std::cerr<<getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if(type == "") {
				std::cerr<<"Could not convert branch type to string when parsing the \"TreeBranches\"' \"onePerEvent\" part."<<std::endl;
				return false;
			} else {
				std::cerr<<"Data type \""<<type<<"\" not supported."<<std::endl;
				return false;
			}
		}
	}
	const unsigned int& N_PARTICLES = antok::Constants::n_particles();
	for(YAML::const_iterator typeIt = perParticleTreeBranches.begin(); typeIt != perParticleTreeBranches.end(); ++typeIt) {
		for(YAML::const_iterator valIt = typeIt->second.begin(); valIt != typeIt->second.end(); ++valIt) {
			std::string type = getYAMLStringSafe(typeIt->first);
			std::string baseName = getYAMLStringSafe(*valIt);
			if(baseName == "") {
				std::cerr<<"Conversion to std::string failed for one of the \"TreeBranches\"' \"onePerParticle\" "<<type<<"s."<<std::endl;
				return false;
			}

			if(type == "double") {
				for(unsigned int i = 0; i < N_PARTICLES; ++i) {
					std::stringstream strStr;
					strStr<<baseName<<(i+1);
					if(not data.insertDouble(strStr.str())) {
						std::cerr<<getVariableInsertionErrorMsg(strStr.str());
						return false;
					}
				}
			} else if(type == "int") {
				for(unsigned int i = 0; i < N_PARTICLES; ++i) {
					std::stringstream strStr;
					strStr<<baseName<<(i+1);
					if(not data.insertInt(strStr.str())) {
						std::cerr<<getVariableInsertionErrorMsg(strStr.str());
						return false;
					}
				}
			} else if(type == "") {
				std::cerr<<"Could not convert branch type to std::string when parsing the \"TreeBranches\"' \"onePerParticle\" part."<<std::endl;
				return false;
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
	std::string treeName = getYAMLStringSafe(config["TreeName"]);
	if(treeName == "") {
		std::cerr<<"Could not convert entry \"TreeName\" to std::string."<<std::endl;
		return false;
	}
	TTree* inTree = dynamic_cast<TTree*>(inFile->Get(treeName.c_str()));
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
			try {
				quantityBaseNames = calcQuantity["Name"].as<std::vector<std::string> >();
			} catch (YAML::TypedBadConversion<std::vector<std::string> > e) {
				std::cerr<<"Could not convert YAML sequence to std::vector<std::string> when parsing \"CalculatedQuantities\"' \"Name\""<<std::endl;
				return false;
			} catch (YAML::TypedBadConversion<std::string> e) {
				std::cerr<<"Could not entries in YAML sequence to std::string when parsing \"CalculatedQuantities\"' \"Name\""<<std::endl;
				return false;
			}
		} else {
			std::string baseName = getYAMLStringSafe(calcQuantity["Name"]);
			if(baseName == "") {
				std::cerr<<"Could not convert one of the \"CalculatedQuantities\"' \"Name\"s to std::string."<<std::endl;
				return false;
			}
			quantityBaseNames.push_back(baseName);
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
			try {
				indices = calcQuantity["Indices"].as<std::vector<int> >();
			} catch (YAML::TypedBadConversion<std::vector<int> > e) {
				std::cerr<<"Could not convert YAML sequence to std::vector<int> when parsing CalculatedQuantities' \""<<quantityBaseNames[0]<<"\" \"Indices\"."<<std::endl;
				return false;
			} catch (YAML::TypedBadConversion<int> e) {
				std::cerr<<"Could not convert entries in YAML sequence to int when parsing CalculatedQuantities' \""<<quantityBaseNames[0]<<"\" \"Indices\"."<<std::endl;
				return false;
			}
		} else {
			indices.push_back(-1);
		}

		const YAML::Node& function = calcQuantity["Function"];
		std::string functionName = getYAMLStringSafe(function["Name"]);

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

			if(functionName == "abs") {
				if(not registerAbs(function, quantityNames, indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "diff") {
				if(not registerDiff(function, quantityNames, indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "getBeamLorentzVector") {
				if(not registerGetBeamLorentzVector(function, quantityNames, indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "getLorentzVec") {
				if(not registerGetLorentzVec(function, quantityNames, indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "getRpdPhi") {
				if(not registerGetRpdPhi(function, quantityNames, indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "getTs") {
				if(not registerGetTs(function, quantityNames, indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "mass") {
				if(not registerMass(function, quantityNames, indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "sum") {
				if(not registerSum(function, quantityNames, indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "sum2") {
				if(not registerSum2(function, quantityNames, indices[indices_i])) {
					return false;
				}
				continue;
			} else if(functionName == "") {
				std::cerr<<"Could not convert function name to std::string for CalculatedQuantity \""<<quantityNames[0]<<"\"."<<std::endl;
				return false;
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

bool antok::Initializer::registerAbs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return false;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Arg", "double"));

	if(not argumentHandler(args, function, index)) {
		std::cerr<<getArgumentHandlerErrorMsg(quantityNames);
		return false;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* argAddr = data.getDoubleAddr(args[0].first);

	if(not data.insertDouble(quantityName)) {
		std::cerr<<getVariableInsertionErrorMsg(quantityNames);
		return false;
	}

	antok::Event& event = antok::ObjectManager::instance()->getEvent();
	event._functions.push_back(new Abs(argAddr, data.getDoubleAddr(quantityName)));

	return true;

};

bool antok::Initializer::registerDiff(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return false;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Minuend", "double"));
	args.push_back(std::pair<std::string, std::string>("Subtrahend", "double"));

	if(not argumentHandler(args, function, index)) {
		std::cerr<<getArgumentHandlerErrorMsg(quantityNames);
		return false;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* arg1Addr = data.getDoubleAddr(args[0].first);
	double* arg2Addr = data.getDoubleAddr(args[1].first);

	if(not data.insertDouble(quantityName)) {
		std::cerr<<getVariableInsertionErrorMsg(quantityNames);
		return false;
	}

	antok::Event& event = antok::ObjectManager::instance()->getEvent();
	event._functions.push_back(new Diff(arg1Addr, arg2Addr, data.getDoubleAddr(quantityName)));

	return true;

};

bool antok::Initializer::registerGetBeamLorentzVector(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return false;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;

	args.push_back(std::pair<std::string, std::string>("dX", "double"));
	args.push_back(std::pair<std::string, std::string>("dY", "double"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not argumentHandler(args, function, index)) {
		std::cerr<<getArgumentHandlerErrorMsg(quantityNames);
		return false;
	}

	double* dXaddr = data.getDoubleAddr(args[0].first);
	double* dYaddr = data.getDoubleAddr(args[1].first);
	TLorentzVector* xLorentzVecAddr = data.getLorentzVectorAddr(args[2].first);

	if(not data.insertLorentzVector(quantityName)) {
		std::cerr<<getVariableInsertionErrorMsg(quantityNames);
		return false;
	}

	antok::Event& event = antok::ObjectManager::instance()->getEvent();
	event._functions.push_back(new antok::GetBeamLorentzVec(dXaddr, dYaddr, xLorentzVecAddr, data.getLorentzVectorAddr(quantityName)));

	return true;

};

bool antok::Initializer::registerGetLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return false;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	bool pType;
	if(function["X"] and function["M"]) {
		pType = false;
	} else if (function["Px"]) {
		pType = true;
	} else {
		std::cerr<<"Function \"getLorentzVec\" needs either variables \"[X, Y, Z, M]\" or \"[Px, Py, Pz, E]\" (variable \""<<quantityName<<"\")."<<std::endl;
		return false;
	}

	std::vector<std::pair<std::string, std::string> > args;

	double* xAddr;
	double* yAddr;
	double* zAddr;
	double* mAddr;
	if(pType) {
		args.push_back(std::pair<std::string, std::string>("Px", "double"));
		args.push_back(std::pair<std::string, std::string>("Py", "double"));
		args.push_back(std::pair<std::string, std::string>("Pz", "double"));
		args.push_back(std::pair<std::string, std::string>("E", "double"));
	} else {
		args.push_back(std::pair<std::string, std::string>("X", "double"));
		args.push_back(std::pair<std::string, std::string>("Y", "double"));
		args.push_back(std::pair<std::string, std::string>("Z", "double"));
		try {
			function["M"].as<double>();
		} catch(YAML::TypedBadConversion<double> e) {
			std::cerr<<"Argument \"M\" in function \"mass\" should be of type double (variable \""<<quantityName<<"\")."<<std::endl;
			return false;
		}
		mAddr = new double();
		(*mAddr) = function["M"].as<double>();
	}

	if(not argumentHandler(args, function, index)) {
		std::cerr<<getArgumentHandlerErrorMsg(quantityNames);
		return false;
	}

	xAddr = data.getDoubleAddr(args[0].first);
	yAddr = data.getDoubleAddr(args[1].first);
	zAddr = data.getDoubleAddr(args[2].first);
	if(pType) {
		mAddr = data.getDoubleAddr(args[3].first);
	}

	if(not data.insertLorentzVector(quantityName)) {
		std::cerr<<getVariableInsertionErrorMsg(quantityNames);
		return false;
	}

	antok::Event& event = antok::ObjectManager::instance()->getEvent();
	event._functions.push_back(new antok::GetLorentzVec(xAddr, yAddr, zAddr, mAddr, data.getLorentzVectorAddr(quantityName), pType));
	return true;

};

bool antok::Initializer::registerGetRpdPhi(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 2) {
		std::cerr<<"Need 2 names for function \"getRpdPhi\""<<std::endl;
		return false;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::string method = getYAMLStringSafe(function["Method"]);
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("RPDProtonLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not argumentHandler(args, function, index)) {
		std::cerr<<getArgumentHandlerErrorMsg(quantityNames);
		return false;
	}

	TLorentzVector* beamLVAddr = data.getLorentzVectorAddr(args[0].first);
	TLorentzVector* RPDprotLVAddr = data.getLorentzVectorAddr(args[1].first);
	TLorentzVector* xLVAddr = data.getLorentzVectorAddr(args[2].first);

	int methodSwitch = -1;
	if(method == "Projection") {
		methodSwitch = 0;
	} else if (method == "Rotation") {
		methodSwitch = 1;
	} else if (method == "") {
	} else {
		if(method == "") {
			std::cerr<<"Could not convert \"Method\" to std::string in function \"getRpdPhi\" when calculating variables \"["<<std::endl;
		} else {
			std::cerr<<"Method \""<<method<<"\" is not supported for function \"getRpdPhi\" when calculating variables \"["<<std::endl;
		}
		for(unsigned int i = 0; i < quantityNames.size()-1; ++i) {
			std::cerr<<quantityNames[i]<<", ";
		}
		std::cerr<<quantityNames[quantityNames.size()-1]<<"]\"."<<std::endl;
		return false;
	}

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insertDouble(quantityNames[i])) {
			std::cerr<<getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return false;
		}
		quantityAddrs.push_back(data.getDoubleAddr(quantityNames[i]));
	}

	antok::Event& event = antok::ObjectManager::instance()->getEvent();
	event._functions.push_back(new GetRpdPhi(beamLVAddr, RPDprotLVAddr, xLVAddr, quantityAddrs[0], quantityAddrs[1], methodSwitch));

	return true;

};

bool antok::Initializer::registerGetTs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 3) {
		std::cerr<<"Need 3 names for function \"getTs\""<<std::endl;
		return false;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;

	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if( not argumentHandler(args, function, index)) {
		std::cerr<<getArgumentHandlerErrorMsg(quantityNames);
		return false;
	}

	TLorentzVector* beamLVAddr = data.getLorentzVectorAddr(args[0].first);
	TLorentzVector* xLVAddr = data.getLorentzVectorAddr(args[1].first);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insertDouble(quantityNames[i])) {
			std::cerr<<getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return false;
		}
		quantityAddrs.push_back(data.getDoubleAddr(quantityNames[i]));
	}

	antok::Event& event = antok::ObjectManager::instance()->getEvent();
	event._functions.push_back(new antok::GetTs(xLVAddr, beamLVAddr, quantityAddrs[0], quantityAddrs[1], quantityAddrs[2]));

	return true;

};

bool antok::Initializer::registerMass(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return false;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));

	if(not argumentHandler(args, function, index)) {
		std::cerr<<getArgumentHandlerErrorMsg(quantityNames);
		return false;
	}

	TLorentzVector* vector = data.getLorentzVectorAddr(args[0].first);
	if(not data.insertDouble(quantityName)) {
		std::cerr<<getVariableInsertionErrorMsg(quantityNames);
		return false;
	}

	antok::Event& event = antok::ObjectManager::instance()->getEvent();
	event._functions.push_back(new antok::Mass(vector, data.getDoubleAddr(quantityName)));
	return true;

};

bool antok::Initializer::registerSum(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return false;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double*> doubleInputAddrs;
	std::vector<int*> intInputAddrs;
	std::vector<Long64_t*> long64_tInputAddrs;
	std::vector<TLorentzVector*> lorentzVectorInputAddrs;
	std::string typeName = "notInitialized";
	std::vector<std::pair<std::string, std::string> > summandNames;

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
		std::vector<int> inner_indices;
		try {
			inner_indices = function["Summands"]["Indices"].as<std::vector<int> >();
		} catch (YAML::TypedBadConversion<std::vector<int> > e) {
			std::cerr<<"Could not convert YAML sequence to std::vector<int> when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
			return false;
		} catch (YAML::TypedBadConversion<int> e) {
			std::cerr<<"Could not convert entries in YAML sequence to int when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
			return false;
		}
		typeName = getYAMLStringSafe(function["Summands"]["Name"]);
		if(typeName == "") {
			std::cerr<<"Could not convert \"Summands\"' \"Name\" to std::string when registering calculation of \""<<quantityName<<"\"."<<std::endl;
		}
		std::string summandBaseName = typeName;
		std::stringstream strStr;
		strStr<<typeName<<inner_indices[0];
		typeName = data.global_map[strStr.str()];
		for(unsigned int inner_indices_i = 0; inner_indices_i < inner_indices.size(); ++inner_indices_i) {
			int inner_index = inner_indices[inner_indices_i];
			std::stringstream strStr;
			strStr<<summandBaseName<<inner_index;
			summandNames.push_back(std::pair<std::string, std::string>(strStr.str(), typeName));
		}
	// Summing over list of variable names
	} else {
		typeName = getYAMLStringSafe(*function["Summands"].begin());
		if(typeName == "") {
			std::cerr<<"Could not convert one of the \"Summands\" to std::string when registering calculation of \""<<quantityName<<"\"."<<std::endl;
			return false;
		}
		if(index > 0) {
			std::stringstream strStr;
			strStr<<typeName<<index;
			typeName = strStr.str();
		}
		typeName = data.global_map[typeName];
		for(YAML::const_iterator summand_it = function["Summands"].begin(); summand_it != function["Summands"].end(); summand_it++) {
			std::string variableName = getYAMLStringSafe(*summand_it);
			if(variableName == "") {
				std::cerr<<"Could not convert one of the \"Summands\" to std::string when registering calculation of \""<<quantityName<<"\"."<<std::endl;
				return false;
			}
			summandNames.push_back(std::pair<std::string, std::string>(variableName, typeName));
		}
	}
	if(not argumentHandler(summandNames, function, index, true)) {
		std::cerr<<getArgumentHandlerErrorMsg(quantityNames);
		return false;
	}

	// Now do type checking and get all the addresses
	for(unsigned int summandNames_i = 0; summandNames_i < summandNames.size(); ++summandNames_i) {

			std::string variableName = summandNames[summandNames_i].first;

			if(typeName == "double") {
				doubleInputAddrs.push_back(data.getDoubleAddr(variableName));
			} else if(typeName == "int") {
				intInputAddrs.push_back(data.getIntAddr(variableName));
			} else if(typeName == "Long64_t") {
				long64_tInputAddrs.push_back(&data.long64_ts[variableName]);
			} else if(typeName == "TLorentzVector") {
				lorentzVectorInputAddrs.push_back(data.getLorentzVectorAddr(variableName));
			} else {
				std::cerr<<"Type \""<<typeName<<"\" is not supported in function \"sum\" (Calculated quantity: \""<<quantityName<<"\")."<<std::endl;
				return false;
			}

		}

	// And produce the function
	antok::Event& event = antok::ObjectManager::instance()->getEvent();
	if(typeName == "double") {
		if(not data.insertDouble(quantityName)) {
			std::cerr<<getVariableInsertionErrorMsg(quantityNames);
			return false;
		}
		event._functions.push_back(new antok::Sum<double>(doubleInputAddrs, data.getDoubleAddr(quantityName)));
	} else if(typeName == "int") {
		if(not data.insertInt(quantityName)) {
			std::cerr<<getVariableInsertionErrorMsg(quantityNames);
			return false;
		}
		event._functions.push_back(new antok::Sum<int>(intInputAddrs, data.getIntAddr(quantityName)));
	} else if(typeName == "Long64_t") {
		if(not data.insertLong64_t(quantityName)) {
			std::cerr<<getVariableInsertionErrorMsg(quantityNames);
			return false;
		}
		event._functions.push_back(new antok::Sum<Long64_t>(long64_tInputAddrs, data.getLong64_tAddr(quantityName)));
	} else if(typeName == "TLorentzVector") {
		if(not data.insertLorentzVector(quantityName)) {
			std::cerr<<getVariableInsertionErrorMsg(quantityNames);
			return false;
		}
		event._functions.push_back(new antok::Sum<TLorentzVector>(lorentzVectorInputAddrs, data.getLorentzVectorAddr(quantityName)));
	}

	return true;

};

bool antok::Initializer::registerSum2(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return false;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double*> doubleInputAddrs;
	std::vector<std::pair<std::string, std::string> > summandNames;

	// Summing over one variable with indices
	if(not function["Summands"]) {
		std::cerr<<"Either trying to sum different types or input variable not found when calculating \""<<quantityName<<"\"."<<std::endl;
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
		std::vector<int> inner_indices;
		try {
			inner_indices = function["Summands"]["Indices"].as<std::vector<int> >();
		} catch (YAML::TypedBadConversion<std::vector<int> > e) {
			std::cerr<<"Could not convert YAML sequence to std::vector<int> when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
			return false;
		} catch (YAML::TypedBadConversion<int> e) {
			std::cerr<<"Could not convert entries in YAML sequence to int when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
			return false;
		}
		std::string summandBaseName = getYAMLStringSafe(function["Summands"]["Name"]);
		if(summandBaseName == "") {
			std::cerr<<"Could not convert \"Summands\"' \"Name\" to std::string when registering calculation of \""<<quantityName<<"\"."<<std::endl;
			return false;
		}
		for(unsigned int inner_indices_i = 0; inner_indices_i < inner_indices.size(); ++inner_indices_i) {
			int inner_index = inner_indices[inner_indices_i];
			std::stringstream strStr;
			strStr<<summandBaseName<<inner_index;
			summandNames.push_back(std::pair<std::string, std::string>(strStr.str(), "double"));
		}
	// Summing over list of variable names
	} else {
		for(YAML::const_iterator summand_it = function["Summands"].begin(); summand_it != function["Summands"].end(); summand_it++) {
			std::string summandName = getYAMLStringSafe(*summand_it);
			if(summandName == "") {
				std::cerr<<"Could not convert one of the \"Summands\" to std::string when registering calculation of \""<<quantityName<<"\"."<<std::endl;
				return false;
			}
			summandNames.push_back(std::pair<std::string, std::string>(summandName, "double"));
		}
	}

	if(not argumentHandler(summandNames, function, index, true)) {
		std::cerr<<getArgumentHandlerErrorMsg(quantityNames);
		return false;
	}

	// Now do type checking and get all the addresses
	for(unsigned int summandNames_i = 0; summandNames_i < summandNames.size(); ++summandNames_i) {

			std::string variableName = summandNames[summandNames_i].first;

			double* addr = data.getDoubleAddr(variableName);
			doubleInputAddrs.push_back(addr);

		}

	// And produce the function
	if(not data.insertDouble(quantityName)) {
		std::cerr<<getVariableInsertionErrorMsg(quantityNames);
		return false;
	}

	antok::Event& event = antok::ObjectManager::instance()->getEvent();
	event._functions.push_back(new antok::Sum2(doubleInputAddrs, data.getDoubleAddr(quantityName)));

	return true;

};

bool antok::Initializer::argumentHandler(std::vector<std::pair<std::string, std::string> >& args, const YAML::Node& function, int index, bool argStringsAlreadyValues)
{

	antok::Data& data = antok::ObjectManager::instance()->getData();
	for(unsigned int i = 0; i < args.size(); ++i) {
		std::string& argName = args[i].first;
		if(not argStringsAlreadyValues) {
			if(not function[argName]) {
				std::cerr<<"Argument \""<<argName<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
				return false;
			}
			argName = getYAMLStringSafe(function[argName]);
			if(argName == "") {
				std::cerr<<"Could not convert one of the arguments to std::string in function \""<<function["Name"]<<"\"."<<std::endl;
				return false;
			}
		}
		if(index > 0) {
			std::stringstream strStr;
			strStr<<argName<<index;
			argName = strStr.str();
		}
		if(data.global_map.count(argName) < 1) {
			std::cerr<<"Argument \""<<argName<<"\" not found in Data's global map."<<std::endl;
			return false;
		}
		if(data.global_map[argName] != args[i].second) {
			std::cerr<<"Argument \""<<argName<<"\" has type \""<<data.global_map[argName]<<"\", expected \""<<args[i].second<<"\"."<<std::endl;
			return false;
		}
	}

	return true;

};

std::string antok::Initializer::getYAMLStringSafe(const YAML::Node& node) {
	try{
		return node.as<std::string>();
	} catch(YAML::TypedBadConversion<std::string> e) {
		return "";
	}
}

std::string antok::Initializer::getArgumentHandlerErrorMsg(std::vector<std::string> quantityNames) {
	std::stringstream msgStream;
	if(quantityNames.size() > 1) {
		msgStream<<"Error when registering calculation for quantities \"[";
		for(unsigned int i = 0; i < quantityNames.size()-1; ++i) {
			msgStream<<quantityNames[i]<<", ";
		}
		msgStream<<quantityNames[quantityNames.size()-1]<<"]\"."<<std::endl;
	} else {
		msgStream<<"Error when registering calculation for quantity \""<<quantityNames[0]<<"\"."<<std::endl;
	}
	return msgStream.str();
};

std::string antok::Initializer::getVariableInsertionErrorMsg(std::vector<std::string> quantityNames, std::string quantityName) {
	std::stringstream msgStream;
	if(quantityNames.size() > 1) {
		msgStream<<"Could not insert variable \""<<quantityName<<"\" when registering calculation for quantities \"[";
		for(unsigned int i = 0; i < quantityNames.size()-1; ++i) {
			msgStream<<quantityNames[i]<<", ";
		}
		msgStream<<quantityNames[quantityNames.size()-1]<<"]\" (double entry?)."<<std::endl;
	} else {
		getVariableInsertionErrorMsg(quantityNames[0]);
	}
	return msgStream.str();
};

std::string antok::Initializer::getVariableInsertionErrorMsg(std::string variableName) {

	std::stringstream msgStream;
	msgStream<<"Could not insert variable \""<<variableName<<"\" (double entry?)."<<std::endl;
	return msgStream.str();

};
