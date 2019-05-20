#include<initializer.h>

#include<sstream>

#include<TFile.h>
#include<TH1D.h>
#include<TTree.h>

#include<constants.h>
#include<cut.hpp>
#include<cutter.h>
#include<entryPoint.hpp>
#include<event.h>
#include<data.h>
#include<functions.hpp>
#include<generators_cuts.h>
#include<generators_functions.h>
#include<generators_plots.h>
#include<object_manager.h>
#include<plotter.h>
#include<yaml_utils.hpp>


/**
	* If no OutputTree node is given in the configuration, the output tree is a clone of the input tree.
	* If OutputTree is given, it has to be a list of variables wich are then stored in all output trees
	*/
TTree* createOutTree( TTree* const inTree, const YAML::Node& config);

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

	using antok::YAMLUtils::hasNodeKey;

	if(_config != 0) {
		std::cerr<<"Attempting to read config file twice."<<std::endl;
		return false;
	}

	_config = new YAML::Node();
	YAML::Node& config = *_config;

	// Load the config file
	try {
		config = YAML::LoadFile(filename);
	} catch (const YAML::ParserException& e) {
		std::cerr<<"Error parsing config file: "<<e.what()<<"."<<std::endl;
		return false;
	}

	// Number of particles
	try {
		if (hasNodeKey(config, "NumberOfParticles")) {
			if (not antok::Constants::set_n_particles(config["NumberOfParticles"].as<unsigned int>())) {
				std::cerr << "Could not set number of particles." << std::endl;
				return false;
			}
		}
	} catch (const YAML::TypedBadConversion<int>& e) {
		std::cerr << "Conversion error when reading \"NumberOfParticles\": " << e.what() << "." << std::endl;
		return false;
	}

	// Get the constants
	YAML::Node constants = config["Constants"];
	try {
		if(not (hasNodeKey(constants, "ChargedKaonMass") and antok::Constants::set_charged_kaon_mass(constants["ChargedKaonMass"].as<double>()))) {
			std::cerr<<"Could not set charged kaon mass."<<std::endl;
			return false;
		}
		if(not (hasNodeKey(constants, "ChargedPionMass") and antok::Constants::set_charged_pion_mass(constants["ChargedPionMass"].as<double>()))) {
			std::cerr<<"Could not set pion mass."<<std::endl;
			return false;
		}
		if(not (hasNodeKey(constants, "ProtonMass") and antok::Constants::set_proton_mass(constants["ProtonMass"].as<double>()))) {
			std::cerr<<"Coud not set proton mass."<<std::endl;
			return false;
		}
	} catch (const YAML::TypedBadConversion<double>& e) {
		std::cerr<<"Conversion error when reading constants: "<<e.what()<<"."<<std::endl;
		return false;
	}

	antok::Constants::_initialized = true;

	return true;

};

bool antok::Initializer::initializeCutter() {

	using antok::YAMLUtils::hasNodeKey;

	if(_config == 0) {
		std::cerr<<"Trying to initialize Cutter without having read the config file first."<<std::endl;
		return false;
	}

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();

	if(objectManager->_cutter != 0) {
		std::cerr<<"Cutter seems to be initialized already."<<std::endl;
		return false;
	}
	objectManager->_cutter = antok::Cutter::instance();
	antok::Cutter& cutter = objectManager->getCutter();
	YAML::Node& config = (*_config);

	if(not hasNodeKey(config, "Cuts")) {
		std::cerr<<"Could not find section \"Cuts\" in configuration file."<<std::endl;
		return false;
	}

	TFile* outFile = objectManager->getOutFile();
	if(outFile == 0) {
		std::cerr<<"Output file not registered."<<std::endl;
		return false;
	}
	outFile->mkdir("tmptmptmp");
	TTree* inTree = objectManager->getInTree();
	if(inTree == 0) {
		std::cerr<<"Input TTree not registered."<<std::endl;
		return false;
	}

	for(YAML::const_iterator cutTrain_it = config["CutTrains"].begin(); cutTrain_it != config["CutTrains"].end(); ++cutTrain_it) {

		const YAML::Node& cutTrain = (*cutTrain_it);
		if(not hasNodeKey(cutTrain, "Name")) {
			std::cerr<<"\"Name\" not found in one of the \"CutTrains\"."<<std::endl;
			return false;
		}
		std::string cutTrainName = antok::YAMLUtils::getString(cutTrain["Name"]);
		if(cutTrainName == "") {
			std::cerr<<"Could not convert \"Name\" to std::string for one of the \"CutTrains\"."<<std::endl;
			return false;
		}
		if(not hasNodeKey(cutTrain, "Cuts")) {
			std::cerr<<"\"Cuts\" not found in cutTrain \""<<cutTrainName<<"\"."<<std::endl;
			return false;
		}
		if(not hasNodeKey(cutTrain, "Pertinent")) {
			std::cerr<<"\"Pertinent\" not found in cutTrain \""<<cutTrainName<<"\"."<<std::endl;
			return false;
		}

		bool pertinent = false;
		{
			std::string pertinence = antok::YAMLUtils::getString(cutTrain["Pertinent"]);
			if(pertinence == "Yes") {
				pertinent = true;
			} else if (pertinence == "No") {
				pertinent = false;
			} else {
				std::cerr<<"Entry \"Pertinent\" has to be either \"Yes\" or \"No\" (cutTrain \""<<cutTrainName<<"\")."<<std::endl;
				return false;
			}
		}

		outFile->cd();
		outFile->mkdir(cutTrainName.c_str());
		outFile->cd("tmptmptmp");
		TDirectory::CurrentDirectory()->mkdir(cutTrainName.c_str());
		if(pertinent) {
			outFile->cd(cutTrainName.c_str());
			TTree* outTree = createOutTree( inTree, config );
			if (outTree == nullptr) return false;
			cutter._outTreeMap[cutTrainName] = outTree;
			assert(objectManager->registerObjectToWrite(TDirectory::CurrentDirectory(), outTree));
		}
		outFile->cd();

		for(YAML::const_iterator cuts_it = cutTrain["Cuts"].begin(); cuts_it != cutTrain["Cuts"].end(); ++cuts_it) {

			const YAML::Node& cutEntry = (*cuts_it);
			std::string shortName = antok::YAMLUtils::getString(cutEntry["ShortName"]);
			if(shortName == "") {
				std::cerr<<"Did not find one of the cut's \"ShortName\"."<<std::endl;
				return false;
			}
			if(cutter._cutTrainsMap[cutTrainName].count(shortName) > 0) {
				std::cerr<<"Cannot have two cuts with the same \"ShortName\" \""<<shortName<<"\"."<<std::endl;
				return false;
			}

			antok::Cut* antokCut = 0;
			bool* result = 0;
			if(not antok::generators::generateCut(cutEntry, antokCut, result)) {
				std::cerr<<"Could not generate cut \""<<shortName<<"\" in cutTrain \""<<cutTrainName<<"\"."<<std::endl;
				return false;
			}

			if(cutter._cutsMap[shortName]) {
				if(not (*cutter._cutsMap[shortName] == *antokCut)) {
					std::cerr<<"Cannot have two different cuts with the same \"ShortName\" \""<<shortName<<"\"."<<std::endl;
					return false;
				}
				delete antokCut;
				antokCut = cutter._cutsMap[shortName];
			} else {
				cutter._cutsMap[shortName] = antokCut;
				cutter._cuts.push_back(std::pair<antok::Cut*, bool*>(antokCut, result));
			}
			cutter._cutTrainsMap[cutTrainName][shortName] = antokCut;
			cutter._cutTrainsCutOrderMap[cutTrainName].push_back(antokCut);

		}

	} // End loop over CutTrains

	for(std::map<std::string, TTree*>::const_iterator outTree_it = cutter._outTreeMap.begin(); outTree_it != cutter._outTreeMap.end(); ++outTree_it) {
		cutter._treesToFill.push_back(std::pair<TTree*, long>(outTree_it->second, cutter.getAllCutsCutmaskForCutTrain(outTree_it->first)));
	}

	return true;

};

bool antok::Initializer::initializeData() {

	using antok::YAMLUtils::hasNodeKey;

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();

	if(objectManager->_data != 0) {
		std::cerr<<"Data seems to be initialized already."<<std::endl;
		return false;
	}

	if(_config == 0) {
		std::cerr<<"Cannot create data object without reading a config file first."<<std::endl;
		return false;
	};
	YAML::Node& config = *_config;

	objectManager->_data = new antok::Data();
	antok::Data& data = objectManager->getData();

	if(not hasNodeKey(config, "TreeBranches")) {
		std::cerr<<"TreeBranches not found in configuration file."<<std::endl;
		return false;
	};
	if(not hasNodeKey(config["TreeBranches"], "onePerEvent")) {
		std::cerr<<"TreeBranch[\"onePerEvent\"] not found in configuration file."<<std::endl;
		return false;
	}

	// Get all the branches in the tree and fill the data maps
	YAML::Node perEventTreeBranches= config["TreeBranches"]["onePerEvent"];
	YAML::Node perParticleTreeBranches= config["TreeBranches"]["onePerParticle"];
	for(YAML::const_iterator typeIt = perEventTreeBranches.begin(); typeIt != perEventTreeBranches.end(); ++typeIt) {
		for(YAML::const_iterator valIt = typeIt->second.begin(); valIt != typeIt->second.end(); ++valIt) {
			std::string type = antok::YAMLUtils::getString(typeIt->first);
			std::string name = antok::YAMLUtils::getString(*valIt);
			if(name == "") {
				std::cerr<<"Conversion to std::string failed for one of the \"TreeBranches\"' \"onePerEvent\" "<<type<<"s."<<std::endl;
				return false;
			}
			if(type == "double") {
				if(not data.insertInputVariable<double>(name)) {
					std::cerr<<antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if(type == "int") {
				if(not data.insertInputVariable<int>(name)) {
					std::cerr<<antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if(type == "Long64_t") {
				if(not data.insertInputVariable<Long64_t>(name)) {
					std::cerr<<antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if(type == "std::vector<double>") {
				if(not data.insertInputVariable<std::vector<double> >(name)) {
					std::cerr<<antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if(type == "TLorentzVector") {
				if(not data.insertInputVariable<TLorentzVector>(name)) {
					std::cerr<<antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if(type == "TVector3") {
				if(not data.insertInputVariable<TVector3>(name)) {
					std::cerr<<antok::Data::getVariableInsertionErrorMsg(name);
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

	const unsigned int N_PARTICLES = antok::Constants::nParticles();
	if (hasNodeKey(config["TreeBranches"], "onePerParticle")) {
		if( N_PARTICLES == 0){
			std::cerr << "Input branches \"onePerParticle\" given, but \"NumberOfParticles\" was not defined!" << std::endl;
			return false;
		}

		for (YAML::const_iterator typeIt = perParticleTreeBranches.begin(); typeIt != perParticleTreeBranches.end(); ++typeIt) {
			for (YAML::const_iterator valIt = typeIt->second.begin(); valIt != typeIt->second.end(); ++valIt) {
				std::string type = antok::YAMLUtils::getString(typeIt->first);
				std::string baseName = antok::YAMLUtils::getString(*valIt);
				if (baseName == "") {
					std::cerr << "Conversion to std::string failed for one of the \"TreeBranches\"' \"onePerParticle\" " << type << "s." << std::endl;
					return false;
				}

				if (type == "double") {
					for (unsigned int i = 0; i < N_PARTICLES; ++i) {
						std::stringstream strStr;
						strStr << baseName << (i + 1);
						if (not data.insertInputVariable<double>(strStr.str())) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(strStr.str());
							return false;
						}
					}
				} else if (type == "int") {
					for (unsigned int i = 0; i < N_PARTICLES; ++i) {
						std::stringstream strStr;
						strStr << baseName << (i + 1);
						if (not data.insertInputVariable<int>(strStr.str())) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(strStr.str());
							return false;
						}
					}
				} else if (type == "Long64_t") {
					for (unsigned int i = 0; i < N_PARTICLES; ++i) {
						std::stringstream strStr;
						strStr << baseName << (i + 1);
						if (not data.insertInputVariable<Long64_t>(strStr.str())) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(strStr.str());
							return false;
						}
					}
				} else if (type == "std::vector<double>") {
					for (unsigned int i = 0; i < N_PARTICLES; ++i) {
						std::stringstream strStr;
						strStr << baseName << (i + 1);
						if (not data.insertInputVariable<std::vector<double> >(strStr.str())) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(strStr.str());
							return false;
						}
					}
				} else if (type == "TLorentzVector") {
					for (unsigned int i = 0; i < N_PARTICLES; ++i) {
						std::stringstream strStr;
						strStr << baseName << (i + 1);
						if (not data.insertInputVariable<TLorentzVector>(strStr.str())) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(strStr.str());
							return false;
						}
					}
				} else if (type == "TVector3") {
					for (unsigned int i = 0; i < N_PARTICLES; ++i) {
						std::stringstream strStr;
						strStr << baseName << (i + 1);
						if (not data.insertInputVariable<TVector3>(strStr.str())) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(strStr.str());
							return false;
						}
					}
				} else if (type == "") {
					std::cerr << "Could not convert branch type to std::string when parsing the \"TreeBranches\"' \"onePerParticle\" part." << std::endl;
					return false;
				} else {
					std::cerr << "Data type \"" << type << "\" not supported." << std::endl;
					return false;
				}
			}
		}
	} else if (N_PARTICLES > 0){
		std::cerr << "\"NumberOfParticles\" given but input branch definition \"onePerParticle\" is missing!" << std::endl;
		return false;
	}


	return true;

}

bool antok::Initializer::initializeInput(){
	using antok::YAMLUtils::hasNodeKey;


	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	antok::Data& data = objectManager->getData();
	YAML::Node& config = *_config;

	// Set the branch addresses of the tree
	if(not hasNodeKey(config, "TreeName")) {
		std::cerr<<"\"TreeName\" not found in configuration file."<<std::endl;
		return false;
	}
	TFile* inFile = objectManager->getInFile();
	std::string treeName = antok::YAMLUtils::getString(config["TreeName"]);
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
		if( data.isInputVariable(it->first))
			inTree->SetBranchAddress(it->first.c_str(), &(it->second));
	}
	for(std::map<std::string, int>::iterator it = data.ints.begin(); it != data.ints.end(); ++it) {
		if( data.isInputVariable(it->first))
			inTree->SetBranchAddress(it->first.c_str(), &(it->second));
	}
	for(std::map<std::string, Long64_t>::iterator it = data.long64_ts.begin(); it != data.long64_ts.end(); ++it) {
		if( data.isInputVariable(it->first))
			inTree->SetBranchAddress(it->first.c_str(), &(it->second));
	}
	for(std::map<std::string, std::vector<double>* >::iterator it = data.doubleVectors.begin(); it != data.doubleVectors.end(); ++it) {
		std::vector<double>* const oldPtr  = it->second; // SetBranchAddress is not allowed to change when opening a (new) file
		if( data.isInputVariable(it->first))
			inTree->SetBranchAddress(it->first.c_str(), &(it->second));
		if(it->second != oldPtr){
			std::cout << "Pointer address of vector<double> '" << it->first << "' has changed while opening a new file." << std::endl;
			return false;
		}
	}
	for(std::map<std::string, TLorentzVector*>::iterator it = data.lorentzVectors.begin(); it != data.lorentzVectors.end(); ++it) {
		TLorentzVector* const oldPtr  = it->second; // SetBranchAddress is not allowed to change when opening a (new) file
		if( data.isInputVariable(it->first))
			inTree->SetBranchAddress(it->first.c_str(), &(it->second));
		if(it->second != oldPtr){
			std::cout << "Pointer address of TLorentzVector '" << it->first << "' has changed while opening a new file." << std::endl;
			return false;
		}
	}
	for(std::map<std::string, TVector3>::iterator it = data.vectors.begin(); it != data.vectors.end(); ++it) {
		if( data.isInputVariable(it->first))
			inTree->SetBranchAddress(it->first.c_str(), &(it->second));
	}

	return true;
}

bool antok::Initializer::updateInput(){
	using antok::YAMLUtils::hasNodeKey;
	bool ok = true;
	ok &= initializeInput(); // again set all input branches

	{ // include waterfall plot of the new input file to the total waterfall plot
		YAML::Node& config = *_config;

		antok::ObjectManager* objectManager = antok::ObjectManager::instance();
		antok::Cutter& cutter = objectManager->getCutter();
		antok::plotUtils::GlobalPlotOptions plotOptions(config["GlobalPlotOptions"]);

		if (plotOptions.statisticsHistInName != "") {
			TH1D* statsHist = dynamic_cast<TH1D*>(objectManager->getInFile()->Get(plotOptions.statisticsHistInName.c_str()));
			if (statsHist != 0)
				objectManager->getPlotter().addInputfileToWaterfallHistograms(statsHist);
			else {
				std::cerr << "Can not find stats histogram in input file!" << std::endl;
				ok = false;
			}
		}

		// set the addresses of all branches of all ouput trees to those of the input tree if no we have no deticated output trees
		if(not hasNodeKey(config, "OutputTree")) { // no definition for the output tree given -> write full input tree
			for( std::map<std::string,TTree *>::iterator it = cutter._outTreeMap.begin(); it != cutter._outTreeMap.end(); ++it){
				objectManager->getInTree()->CopyAddresses(it->second);

			}


		}

	}

	return ok;
}

bool antok::Initializer::initializeEvent() {

	using antok::YAMLUtils::hasNodeKey;

	if(_config == 0) {
		std::cerr<<"Trying to initialize Cutter without having read the config file first."<<std::endl;
		return false;
	}
	YAML::Node& config = *_config;

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();

	if(objectManager->_event != 0) {
		std::cerr<<"Event seems to be initialized already."<<std::endl;
		return false;
	}
	objectManager->_event = antok::Event::instance();

	if(not hasNodeKey(config, "CalculatedQuantities")) {
		std::cerr<<"Warning: \"CalculatedQuantities\" not found in configuration file."<<std::endl;
	}
	for(YAML::const_iterator calcQuantity_it = config["CalculatedQuantities"].begin(); calcQuantity_it != config["CalculatedQuantities"].end(); ++calcQuantity_it) {

		YAML::Node calcQuantity = (*calcQuantity_it);
		std::vector<std::string> quantityBaseNames;
		if(not hasNodeKey(calcQuantity, "Name")) {
			std::cerr<<"Could not convert a \"Name\" of an entry in \"CalculatedQuantities\"."<<std::endl;
			return false;
		}
		if(calcQuantity["Name"].IsSequence()) {
			try {
				quantityBaseNames = calcQuantity["Name"].as<std::vector<std::string> >();
			} catch (const YAML::TypedBadConversion<std::vector<std::string> >& e) {
				std::cerr<<"Could not convert YAML sequence to std::vector<std::string> when parsing \"CalculatedQuantities\"' \"Name\"."<<std::endl;
				return false;
			} catch (const YAML::TypedBadConversion<std::string>& e) {
				std::cerr<<"Could not entries in YAML sequence to std::string when parsing \"CalculatedQuantities\"' \"Name\"."<<std::endl;
				return false;
			}
		} else {
			std::string baseName = antok::YAMLUtils::getString(calcQuantity["Name"]);
			if(baseName == "") {
				std::cerr<<"Could not convert one of the \"CalculatedQuantities\"' \"Name\"s to std::string."<<std::endl;
				return false;
			}
			quantityBaseNames.push_back(baseName);
		}
		if(quantityBaseNames.empty()) {
			std::cerr<<"Did not find a name to save calculated quantity to."<<std::endl;
			return false;
		}

		if(not (hasNodeKey(calcQuantity, "Function") and hasNodeKey(calcQuantity["Function"], "Name"))) {
			std::cerr<<"No Function or no function name for calculated quantity \""<<quantityBaseNames[0]<<"\"."<<std::endl;
			return false;
		}

		std::vector<int> indices;
		if(hasNodeKey(calcQuantity, "Indices")) {
			try {
				indices = calcQuantity["Indices"].as<std::vector<int> >();
			} catch (const YAML::TypedBadConversion<std::vector<int> >& e) {
				std::cerr<<"Could not convert YAML sequence to std::vector<int> when parsing CalculatedQuantities' \""<<quantityBaseNames[0]<<"\" \"Indices\"."<<std::endl;
				return false;
			} catch (const YAML::TypedBadConversion<int>& e) {
				std::cerr<<"Could not convert entries in YAML sequence to int when parsing CalculatedQuantities' \""<<quantityBaseNames[0]<<"\" \"Indices\"."<<std::endl;
				return false;
			}
		} else {
			indices.push_back(-1);
		}

		const YAML::Node& function = calcQuantity["Function"];
		std::string functionName = antok::YAMLUtils::getString(function["Name"]);

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

			antok::Function* antokFunctionPtr = 0;
			if(functionName == "abs") {
				antokFunctionPtr = antok::generators::generateAbs(function, quantityNames, indices[indices_i]);
			} else if(functionName == "log") {
				antokFunctionPtr = antok::generators::generateLog(function, quantityNames, indices[indices_i]);
			} else if(functionName == "convertIntToDouble") {
				antokFunctionPtr = antok::generators::generateConvertIntToDouble(function, quantityNames, indices[indices_i]);
			} else if(functionName == "diff") {
				antokFunctionPtr = antok::generators::generateDiff(function, quantityNames, indices[indices_i]);
			} else if(functionName == "quotient") {
				antokFunctionPtr = antok::generators::generateQuotient(function, quantityNames, indices[indices_i]);
			} else if(functionName == "mul") {
				antokFunctionPtr = antok::generators::generateMul(function, quantityNames, indices[indices_i]);
			} else if(functionName == "energy") {
				antokFunctionPtr = antok::generators::generateEnergy(function, quantityNames, indices[indices_i]);
			} else if(functionName == "getBeamLorentzVector") {
				antokFunctionPtr = antok::generators::generateGetBeamLorentzVector(function, quantityNames, indices[indices_i]);
			} else if(functionName == "getGradXGradY") {
				antokFunctionPtr = antok::generators::generateGetGradXGradY(function, quantityNames, indices[indices_i]);
			} else if(functionName == "getLorentzVectorAttributes") {
				antokFunctionPtr = antok::generators::generateGetLorentzVectorAttributes(function, quantityNames, indices[indices_i]);
			} else if(functionName == "getLorentzVec") {
				antokFunctionPtr = antok::generators::generateGetLorentzVec(function, quantityNames, indices[indices_i]);
			} else if(functionName == "getTs") {
				antokFunctionPtr = antok::generators::generateGetTs(function, quantityNames, indices[indices_i]);
			} else if(functionName == "getVector3") {
				antokFunctionPtr = antok::generators::generateGetVector3(function, quantityNames, indices[indices_i]);
			} else if(functionName == "mass") {
				antokFunctionPtr = antok::generators::generateMass(function, quantityNames, indices[indices_i]);
			} else if(functionName == "radToDegree") {
				antokFunctionPtr = antok::generators::generateRadToDegree(function, quantityNames, indices[indices_i]);
			} else if(functionName == "sum") {
				antokFunctionPtr = antok::generators::generateSum(function, quantityNames, indices[indices_i]);
			} else if(functionName == "sum2") {
				antokFunctionPtr= antok::generators::generateSum2(function, quantityNames, indices[indices_i]);
			} else if(functionName == "") {
				std::cerr<<"Could not convert function name to std::string for CalculatedQuantity \""<<quantityNames[0]<<"\"."<<std::endl;
				return false;
			} else {
				antokFunctionPtr = antok::user::getUserFunction(function, quantityNames, indices[indices_i]);
				if(antokFunctionPtr == 0) {
					std::cerr<<"Function type \""<<functionName<<"\" not supported."<<std::endl;
					return false;
				}
			}
			if(antokFunctionPtr == 0) {
				std::cerr<<"Error initializing function with name \"" + functionName + "\" which should calculate [";
				for(unsigned int i = 0; i < (quantityNames.size() - 1); ++i) {
					std::cerr<<quantityNames[i]<<", ";
				}
				std::cerr<<quantityNames[quantityNames.size() - 1]<<"]."<<std::endl;
				return false;
			}
			antok::Event& event = objectManager->getEvent();
			event._functions.push_back(antokFunctionPtr);

		}

	}

	return true;

};

bool antok::Initializer::initializePlotter() {

	using antok::YAMLUtils::hasNodeKey;

	if(_config == 0) {
		std::cerr<<"Trying to initialize Plotter without having read the config file first."<<std::endl;
		return false;
	}
	YAML::Node& config = *_config;

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();

	if(objectManager->_plotter != 0) {
		std::cerr<<"Event seems to be initialized already."<<std::endl;
		return false;
	}
	objectManager->_plotter = antok::Plotter::instance();

	antok::Plotter& plotter = antok::ObjectManager::instance()->getPlotter();

	antok::plotUtils::GlobalPlotOptions plotOptions(config["GlobalPlotOptions"]);

	if((plotOptions.statisticsHistInName == "" or plotOptions.statisticsHistOutName == "") and
	   not (plotOptions.statisticsHistInName == "" and plotOptions.statisticsHistOutName == ""))
	{
		std::cerr<<"Something went wrong when parsing the \"GlobalPlotOptions\"."<<std::endl;
		return false;
	}

	if(plotOptions.statisticsHistInName != "") {
		TFile* inFile = objectManager->getInFile();
		TH1D* statsHistTemplate = dynamic_cast<TH1D*>(inFile->Get(plotOptions.statisticsHistInName.c_str()));
		if(statsHistTemplate == 0) {
			std::cerr<<"Could not get the input \"StatisticsHistogram\" from the input file."<<std::endl;
			return false;
		}
		antok::Cutter& cutter = objectManager->getCutter();
		TFile* outFile = objectManager->getOutFile();
		std::vector<antok::plotUtils::waterfallHistogramContainer> waterfallHists;
		for(std::map<std::string, std::vector<antok::Cut*> >::const_iterator cutTrain_it = cutter._cutTrainsCutOrderMap.begin();
		    cutTrain_it != cutter._cutTrainsCutOrderMap.end();
		    ++cutTrain_it)
		{
			std::string cutTrainName = cutTrain_it->first;
			outFile->cd(cutTrainName.c_str());
			TH1D* statsHist = (TH1D*)statsHistTemplate->Clone(plotOptions.statisticsHistOutName.c_str());
			if(statsHist == 0) {
				std::cerr<<"Could not generate the output \"StatisticsHistogram\" for \"CutTrain\" \""<<cutTrainName<<"\"."<<std::endl;
				return false;
			}
			objectManager->registerObjectToWrite(TDirectory::CurrentDirectory(), statsHist);
			std::vector<antok::Cut*> cuts = cutTrain_it->second;
			std::vector<std::pair<const char*, const bool*> > cutsAndResults;
			for(unsigned int i = 0; i < cuts.size(); ++i) {
				cutsAndResults.push_back(std::pair<const char*, const bool*>(cuts[i]->getLongName().c_str(), cutter.getCutResult(cuts[i])));
			}
			waterfallHists.push_back(antok::plotUtils::waterfallHistogramContainer(statsHist, cutsAndResults));
		}
		plotter._waterfallHistograms = waterfallHists;
		outFile->cd();
	}

	if(not hasNodeKey(config, "Plots")) {
		std::cerr<<"Warning: \"Plots\" not found in configuration file."<<std::endl;
	}
	for(YAML::const_iterator plots_it = config["Plots"].begin(); plots_it != config["Plots"].end(); ++plots_it) {

		const YAML::Node& plot = *plots_it;

		if(not hasNodeKey(plot, "Name")) {
			std::cerr<<"\"Name\" not found for one of the \"Plots\"."<<std::endl;
			return false;
		}
		std::string plotName = antok::YAMLUtils::getString(plot["Name"]);

        if( plotName.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_ ") != std::string::npos )
        {
            std::cerr<<"Invalid character in name for \"Plot\" \""<<plotName<<"\"."<<std::endl;
            return false;
        }
		if((not((hasNodeKey(plot, "Variable") and hasNodeKey(plot, "LowerBound") and hasNodeKey(plot, "UpperBound")) or
		       (hasNodeKey(plot, "Variables") and hasNodeKey(plot, "LowerBounds") and hasNodeKey(plot, "UpperBounds")))) or
			not hasNodeKey(plot, "NBins"))
		{
			std::cerr<<"Required variables not found for \"Plot\" \""<<plotName;
			std::cerr<<"\" (\"Variable(s)\"|\"LowerBound(s)\"|\"UpperBound(s)\"|\"NBins\")."<<std::endl;
			return false;
		}

		antok::Plot* antokPlot = 0;
		if(hasNodeKey(plot, "Variables")) {
			if(plot["Variables"].IsSequence()) {
				if(not (plot["Variables"].IsSequence() and
				        plot["LowerBounds"].IsSequence() and
				        plot["UpperBounds"].IsSequence() and
						plot["NBins"].IsSequence()))
				{
					std::cerr<<"\"Variables\", \"LowerBounds\", \"UpperBounds\" and \"NBins\" all need to be sequences (in \"Plot\" \""<<plotName<<"\")."<<std::endl;
					return false;
				}
				if(plot["Variables"].size() > 3) {
					std::cerr<<"Cannot have \"Plot\" \""<<plotName<<"\" with more than 2 \"Variables\""<<std::endl;
					return false;
				}
				if((plot["Variables"].size() != plot["LowerBounds"].size()) or
				   (plot["Variables"].size() != plot["UpperBounds"].size()) or
				   (plot["Variables"].size() != plot["NBins"].size())) {
					std::cerr<<"\"Variables\", \"LowerBounds\", \"UpperBounds\" and \"NBins\" need to have the same number of entries (in \"Plot\" \""<<plotName<<"\")."<<std::endl;
					return false;
				}
				if(plot["Variables"].size() == 2) {
					antokPlot = antok::generators::generate2DPlot(plot, plotOptions);
				} else if (plot["Variables"].size() == 1) {
					antokPlot = antok::generators::generate1DPlot(plot, plotOptions);
				} else if(plot["Variables"].size() == 3) {
					antokPlot = antok::generators::generate3DPlot(plot, plotOptions);
				} else {
					std::cerr<<"Empty \"Variables\" sequence found in \"Plot\" \""<<plotName<<"\"."<<std::endl;
					return false;
				}
			} else {
				std::cerr<<"\"Variables\" implies 2D/3D plot and needs to be a list (in \"Plot\" \""<<plotName<<"\")."<<std::endl;
				return false;
			}
		} else if(hasNodeKey(plot, "Variable")) {
			antokPlot = antok::generators::generate1DPlot(plot, plotOptions);
		} else {
			assert(false);
		}

		if(antokPlot == 0) {
			std::cerr<<"Could not generate \"Plot\" \""<<plotName<<"\"."<<std::endl;
			return false;
		}
		plotter._plots.push_back(antokPlot);
	}

	return true;

};


template< typename T>
void addToOutputBranch(TTree* const outTree, antok::Data& data, const std::string& variable_name){
	outTree->Branch(variable_name.c_str(), data.getAddr<T>(variable_name) );
}

TTree* createOutTree( TTree* const inTree, const YAML::Node& config ){
	using antok::YAMLUtils::hasNodeKey;

	TTree* outTree = nullptr;
	if(not hasNodeKey(config, "OutputTree")) { // no definition for the output tree given -> write full input tree
		outTree = inTree->CloneTree(0);
	} else { // own definition of output tree is given -> build output tree
		antok::Data& data = antok::ObjectManager::instance()->getData();
		outTree = new TTree( inTree->GetName(), inTree->GetTitle() );

		for( const auto variable_node: config["OutputTree"] ){
			const std::string variable_name = antok::YAMLUtils::getString( variable_node );
			const std::string variable_type = data.getType( variable_name );
			if( variable_name != "" ){
				if(      variable_type == "double" )	        addToOutputBranch<double>(         outTree, data, variable_name );
				else if( variable_type == "int" )	        addToOutputBranch<int>(            outTree, data, variable_name );
				else if( variable_type == "TVector3" )       addToOutputBranch<TVector3>(       outTree, data, variable_name );
				else if( variable_type == "TLorentzVector" ) addToOutputBranch<TLorentzVector>( outTree, data, variable_name );
				else if( variable_type == "std::vector<double>" ) addToOutputBranch<std::vector<double>>( outTree, data, variable_name );
				else {
					std::cerr << "Variable type \"" << variable_type << "\" of variable \"" << variable_name << "\" not implemented for own output tree." << std::endl;
					return nullptr;
				}

			} else {
				std::cerr<<"Variable \""<<variable_name<<"\" not found in Data's global map for the output tree."<<std::endl;
				return nullptr;
			}
		}

	}
	return outTree;
}
