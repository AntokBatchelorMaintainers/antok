#include <sstream>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include "constants.h"
#include "cut.hpp"
#include "cutter.h"
#include "entryPoint.hpp"
#include "event.h"
#include "data.h"
#include "functions.hpp"
#include "generators_cuts.h"
#include "generators_plots.h"
#include "initializer.h"
#include "object_manager.h"
#include "plotter.h"
#include "yaml_utils.hpp"


/**
	* If no OutputTree node is given in the configuration, the output tree is a clone of the input tree.
	* If OutputTree is given, it has to be a list of variables wich are then stored in all output trees
	*/
TTree* createOutTree(TTree* const      inTree,
                     const YAML::Node& config);


antok::Initializer* antok::Initializer::_initializer = 0;


antok::Initializer*
antok::Initializer::instance()
{
	if (_initializer == nullptr) {
		_initializer = new antok::Initializer();
	}
	return _initializer;
}


antok::Initializer::Initializer()
	: _config(nullptr)
{ }


bool
antok::Initializer::readConfigFile(const std::string& filename)
{
	using antok::YAMLUtils::hasNodeKey;
	if (_config != nullptr) {
		std::cerr << "Attempting to read config file twice." << std::endl;
		return false;
	}

	// Load the config file
	_config = new YAML::Node();
	YAML::Node& config = *_config;
	try {
		config = YAML::LoadFile(filename);
	} catch (const YAML::ParserException& e) {
		std::cerr << "Error parsing config file: " << e.what() << "." << std::endl;
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
		std::cerr << "Conversion error when reading 'NumberOfParticles': " << e.what() << "." << std::endl;
		return false;
	}

	// Get the constants
	YAML::Node constants = config["Constants"];
	try {
		if (not (hasNodeKey(constants, "ChargedKaonMass") and antok::Constants::set_charged_kaon_mass(constants["ChargedKaonMass"].as<double>()))) {
			std::cerr << "Could not set charged kaon mass." << std::endl;
			return false;
		}
		if (not (hasNodeKey(constants, "ChargedPionMass") and antok::Constants::set_charged_pion_mass(constants["ChargedPionMass"].as<double>()))) {
			std::cerr << "Could not set pion mass." << std::endl;
			return false;
		}
		if (not (hasNodeKey(constants, "ProtonMass") and antok::Constants::set_proton_mass(constants["ProtonMass"].as<double>()))) {
			std::cerr << "Coud not set proton mass." << std::endl;
			return false;
		}
	} catch (const YAML::TypedBadConversion<double>& e) {
		std::cerr << "Conversion error when reading constants: " << e.what() << "." << std::endl;
		return false;
	}

	antok::Constants::_initialized = true;
	return true;
}


bool
antok::Initializer::initializeCutter()
{
	using antok::YAMLUtils::hasNodeKey;
	if (_config == nullptr) {
		std::cerr << "Trying to initialize Cutter without having read the config file first." << std::endl;
		return false;
	}

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	if (objectManager->_cutter != nullptr) {
		std::cerr << "Cutter seems to be initialized already." << std::endl;
		return false;
	}
	objectManager->_cutter = antok::Cutter::instance();
	antok::Cutter& cutter = objectManager->getCutter();
	YAML::Node&    config = (*_config);
	if (not hasNodeKey(config, "Cuts")) {
		std::cerr << "Could not find section 'Cuts' in configuration file." << std::endl;
		return false;
	}

	TFile* outFile = objectManager->getOutFile();
	if (outFile == nullptr) {
		std::cerr << "Output file not registered." << std::endl;
		return false;
	}
	outFile->mkdir("tmptmptmp");
	TTree* inTree = objectManager->getInTree();
	if (inTree == nullptr) {
		std::cerr << "Input TTree not registered." << std::endl;
		return false;
	}

	for (const YAML::Node& cutTrain : config["CutTrains"]) {
		if (not hasNodeKey(cutTrain, "Name")) {
			std::cerr << "'Name' not found in one of the 'CutTrains'." << std::endl;
			return false;
		}
		const std::string cutTrainName = antok::YAMLUtils::getString(cutTrain["Name"]);
		if (cutTrainName == "") {
			std::cerr << "Could not convert 'Name' to std::string for one of the 'CutTrains'." << std::endl;
			return false;
		}
		if (not hasNodeKey(cutTrain, "Cuts")) {
			std::cerr << "'Cuts' not found in cutTrain '" << cutTrainName << "'." << std::endl;
			return false;
		}
		if (not hasNodeKey(cutTrain, "Pertinent")) {
			std::cerr << "'Pertinent' not found in cutTrain '" << cutTrainName << "'." << std::endl;
			return false;
		}

		bool pertinent = false;
		{
			const std::string pertinence = antok::YAMLUtils::getString(cutTrain["Pertinent"]);
			if (pertinence == "Yes") {
				pertinent = true;
			} else if (pertinence == "No") {
				pertinent = false;
			} else {
				std::cerr << "Entry 'Pertinent' has to be either 'Yes' or 'No' (cutTrain '" << cutTrainName << "')." << std::endl;
				return false;
			}
		}

		outFile->cd();
		outFile->mkdir(cutTrainName.c_str());
		outFile->cd("tmptmptmp");
		TDirectory::CurrentDirectory()->mkdir(cutTrainName.c_str());
		if (pertinent) {
			outFile->cd(cutTrainName.c_str());
			TTree* outTree = createOutTree(inTree, config);
			if (outTree == nullptr) {
				return false;
			}
			cutter._outTreeMap[cutTrainName] = outTree;
			assert(objectManager->registerObjectToWrite(TDirectory::CurrentDirectory(), outTree));
		}
		outFile->cd();

		for (const YAML::Node& cutEntry : cutTrain["Cuts"]) {
			std::string shortName = antok::YAMLUtils::getString(cutEntry["ShortName"]);
			if (shortName == "") {
				std::cerr << "Did not find one of the cut's 'ShortName'." << std::endl;
				return false;
			}
			if (cutter._cutTrainsMap[cutTrainName].count(shortName) > 0) {
				std::cerr << "Cannot have two cuts with the same 'ShortName' '" << shortName << "'." << std::endl;
				return false;
			}

			antok::Cut* antokCut = nullptr;
			bool*       result   = nullptr;
			if (not antok::generators::generateCut(cutEntry, antokCut, result)) {
				std::cerr << "Could not generate cut '" << shortName << "' in cutTrain '" << cutTrainName << "'." << std::endl;
				return false;
			}
			if (cutter._cutsMap[shortName]) {
				if (not (*cutter._cutsMap[shortName] == *antokCut)) {
					std::cerr << "Cannot have two different cuts with the same 'ShortName' '" << shortName << "'." << std::endl;
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
		}  // End loop over Cuts

	}  // End loop over CutTrains

	for (const std::pair<std::string, TTree*>& outTree : cutter._outTreeMap) {
		cutter._treesToFill.push_back(std::pair<TTree*, long>(outTree.second, cutter.getAllCutsCutmaskForCutTrain(outTree.first)));
	}

	return true;
}


bool
antok::Initializer::initializeData()
{
	using antok::YAMLUtils::hasNodeKey;
	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	if (objectManager->_data != nullptr) {
		std::cerr << "Data seems to be initialized already." << std::endl;
		return false;
	}
	if (_config == nullptr) {
		std::cerr << "Cannot create data object without reading a config file first." << std::endl;
		return false;
	}
	YAML::Node& config = *_config;

	if (not hasNodeKey(config, "TreeBranches")) {
		std::cerr << "TreeBranches not found in configuration file." << std::endl;
		return false;
	}
	if (not hasNodeKey(config["TreeBranches"], "onePerEvent")) {
		std::cerr << "TreeBranches['onePerEvent'] not found in configuration file." << std::endl;
		return false;
	}

	// Get all the branches in the tree and fill the data maps
	objectManager->_data = new antok::Data();
	antok::Data& data = objectManager->getData();
	YAML::Node perEventTreeBranches    = config["TreeBranches"]["onePerEvent"];
	YAML::Node perParticleTreeBranches = config["TreeBranches"]["onePerParticle"];
	for (const YAML::detail::iterator_value& treeBranch : perEventTreeBranches) {
		for (const YAML::Node& typeNode : treeBranch.second) {
			std::string type = antok::YAMLUtils::getString(treeBranch.first);
			std::string name = antok::YAMLUtils::getString(typeNode);
			if (name == "") {
				std::cerr << "Conversion to std::string failed for one of the TreeBranches['onePerEvent'] " << type << "s." << std::endl;
				return false;
			}
			if (type == "double") {
				if (not data.insertInputVariable<double>(name)) {
					std::cerr << antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if (type == "int") {
				if (not data.insertInputVariable<int>(name)) {
					std::cerr << antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if (type == "Long64_t") {
				if (not data.insertInputVariable<Long64_t>(name)) {
					std::cerr << antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if (type == "std::vector<double>") {
				if (not data.insertInputVariable<std::vector<double>>(name)) {
					std::cerr << antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if (type == "std::vector<int>") {
				if (not data.insertInputVariable<std::vector<int>>(name)) {
					std::cerr << antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if (type == "TLorentzVector") {
				if (not data.insertInputVariable<TLorentzVector>(name)) {
					std::cerr << antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if (type == "TVector3") {
				if (not data.insertInputVariable<TVector3>(name)) {
					std::cerr << antok::Data::getVariableInsertionErrorMsg(name);
					return false;
				}
			} else if (type == "") {
				std::cerr << "Could not convert branch type to string when parsing the TreeBranches['onePerEvent'] part." << std::endl;
				return false;
			} else {
				std::cerr << "Data type '" << type << "' not supported." << std::endl;
				return false;
			}
		}
	}

	// lambda to construct variable names
	auto variableName = [] (const std::string& baseName, const int i) -> std::string
	{
		std::stringstream varName;
		varName << baseName << i;
		return varName.str();
	};

	const size_t nParticles = antok::Constants::nParticles();
	if (hasNodeKey(config["TreeBranches"], "onePerParticle")) {
		if (nParticles == 0) {
			std::cerr << "Input branches TreeBranches['onePerParticle'] given, but 'NumberOfParticles' was not defined!" << std::endl;
			return false;
		}

		for (const YAML::detail::iterator_value& treeBranch : perParticleTreeBranches) {
			for (const YAML::Node& typeNode : treeBranch.second) {
				std::string type     = antok::YAMLUtils::getString(treeBranch.first);
				std::string baseName = antok::YAMLUtils::getString(typeNode);
				if (baseName == "") {
					std::cerr << "Conversion to std::string failed for one of the TreeBranches['onePerParticle'] " << type << "s." << std::endl;
					return false;
				}
				if (type == "double") {
					for (size_t i = 0; i < nParticles; ++i) {
						const std::string varName = variableName(baseName, i + 1);
						if (not data.insertInputVariable<double>(varName)) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(varName);
							return false;
						}
					}
				} else if (type == "int") {
					for (size_t i = 0; i < nParticles; ++i) {
						const std::string varName = variableName(baseName, i + 1);
						if (not data.insertInputVariable<int>(varName)) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(varName);
							return false;
						}
					}
				} else if (type == "Long64_t") {
					for (size_t i = 0; i < nParticles; ++i) {
						const std::string varName = variableName(baseName, i + 1);
						if (not data.insertInputVariable<Long64_t>(varName)) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(varName);
							return false;
						}
					}
				} else if (type == "std::vector<double>") {
					for (size_t i = 0; i < nParticles; ++i) {
						const std::string varName = variableName(baseName, i + 1);
						if (not data.insertInputVariable<std::vector<double>>(varName)) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(varName);
							return false;
						}
					}
				} else if (type == "std::vector<int>") {
					for (size_t i = 0; i < nParticles; ++i) {
						const std::string varName = variableName(baseName, i + 1);
						if (not data.insertInputVariable<std::vector<int>>(varName)) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(varName);
							return false;
						}
					}
				} else if (type == "TLorentzVector") {
					for (size_t i = 0; i < nParticles; ++i) {
						const std::string varName = variableName(baseName, i + 1);
						if (not data.insertInputVariable<TLorentzVector>(varName)) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(varName);
							return false;
						}
					}
				} else if (type == "TVector3") {
					for (size_t i = 0; i < nParticles; ++i) {
						const std::string varName = variableName(baseName, i + 1);
						if (not data.insertInputVariable<TVector3>(varName)) {
							std::cerr << antok::Data::getVariableInsertionErrorMsg(varName);
							return false;
						}
					}
				} else if (type == "") {
					std::cerr << "Could not convert branch type to std::string when parsing the TreeBranches['onePerParticle'] part." << std::endl;
					return false;
				} else {
					std::cerr << "Data type '" << type << "' not supported." << std::endl;
					return false;
				}
			}
		}  // loop over perParticleTreeBranches
	} else if (nParticles > 0) {
		std::cerr << "'NumberOfParticles' given but input branch definition TreeBranches['onePerParticle'] is missing!" << std::endl;
		return false;
	}  // onePerParticle branch

	return true;
}


bool
antok::Initializer::initializeInput()
{
	using antok::YAMLUtils::hasNodeKey;
	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	antok::Data&          data          = objectManager->getData();
	YAML::Node&           config        = *_config;

	// Set the branch addresses of the tree
	if (not hasNodeKey(config, "TreeName")) {
		std::cerr << "'TreeName' not found in configuration file." << std::endl;
		return false;
	}
	TFile* inFile = objectManager->getInFile();
	std::string treeName = antok::YAMLUtils::getString(config["TreeName"]);
	if (treeName == "") {
		std::cerr << "Could not convert entry 'TreeName' to std::string." << std::endl;
		return false;
	}
	TTree* inTree = dynamic_cast<TTree*>(inFile->Get(treeName.c_str()));
	if (inTree == nullptr) {
		std::cerr << "Could not open input TTree." << std::endl;
		return false;
	}
	objectManager->_inTree = inTree;

	for (auto& it : data._doubles) {
		if (data.isInputVariable(it.first))
			inTree->SetBranchAddress(it.first.c_str(), &(it.second));
	}
	for (auto& it : data._ints) {
		if (data.isInputVariable(it.first))
			inTree->SetBranchAddress(it.first.c_str(), &(it.second));
	}
	for (auto& it : data._long64_ts) {
		if (data.isInputVariable(it.first))
			inTree->SetBranchAddress(it.first.c_str(), &(it.second));
	}
	for (auto& it : data._doubleVectors) {
		//std::vector<double>* const oldPtr = it.second;  // SetBranchAddress is not allowed to change when opening a (new) file
		if (data.isInputVariable(it.first))
			inTree->SetBranchAddress(it.first.c_str(), &(it.second));
		/*if (it.second != oldPtr) {  //TODO unclear why this obscure test is needed
			std::cout << "Pointer address of vector<double> '" << it.first << "' has changed while opening a new file." << std::endl;
			return false;
		*/
	}
	for (auto& it : data._intVectors) {
		//std::vector<int>* const oldPtr = it.second;  // SetBranchAddress is not allowed to change when opening a (new) file
		if (data.isInputVariable(it.first))
			inTree->SetBranchAddress(it.first.c_str(), &(it.second));
		/*if (it.second != oldPtr) {  //TODO unclear why this obscure test is needed
			std::cout << "Pointer address of vector<int> '" << it.first << "' has changed while opening a new file." << std::endl;
			return false;
		*/
	}
	for (auto& it : data._lorentzVectors) {
		//TLorentzVector* const oldPtr = it.second;  // SetBranchAddress is not allowed to change when opening a (new) file
		if (data.isInputVariable(it.first))
			inTree->SetBranchAddress(it.first.c_str(), &(it.second));
		/*if (it.second != oldPtr) {  //TODO unclear why this obscure test is needed
			std::cout << "Pointer address of TLorentzVector '" << it.first << "' has changed while opening a new file." << std::endl;
			return false;
		}*/
	}
	for (auto& it : data._vector3s) {
		if (data.isInputVariable(it.first))
			inTree->SetBranchAddress(it.first.c_str(), &(it.second));
	}

	return true;
}


bool
antok::Initializer::updateInput()
{
	using antok::YAMLUtils::hasNodeKey;
	bool ok = true;
	ok &= initializeInput();  // again set all input branches

	// include waterfall plot of the new input file to the existing waterfall plot
	const YAML::Node&     config        = *_config;
	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	antok::Cutter&        cutter        = objectManager->getCutter();
	antok::plotUtils::GlobalPlotOptions plotOptions(config["GlobalPlotOptions"]);
	if (plotOptions.statisticsHistInName != "") {
		TH1D* statsHist = dynamic_cast<TH1D*>(objectManager->getInFile()->Get(plotOptions.statisticsHistInName.c_str()));
		if (statsHist != nullptr) {
			objectManager->getPlotter().addInputfileToWaterfallHistograms(statsHist);
		} else {
			std::cerr << "Cannot find stats histogram in input file!" << std::endl;
			ok = false;
		}
	}

	// set the addresses of all branches of all ouput trees to those of
	// the input tree if no we have no dedicated output trees
	if (not hasNodeKey(config, "OutputTree")) {  // no definition for the output tree given -> write full input tree
		for (std::pair<std::string,TTree *> outTree : cutter._outTreeMap) {
			objectManager->getInTree()->CopyAddresses(outTree.second);
		}
	}

	return ok;
}


bool
antok::Initializer::initializeEvent()
{
	using antok::YAMLUtils::hasNodeKey;
	if (_config == nullptr) {
		std::cerr << "Trying to initialize Cutter without having read the config file first." << std::endl;
		return false;
	}
	YAML::Node&           config        = *_config;
	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	if (objectManager->_event != nullptr) {
		std::cerr << "Event seems to be initialized already." << std::endl;
		return false;
	}
	objectManager->_event = antok::Event::instance();

	if (not hasNodeKey(config, "CalculatedQuantities")) {
		std::cerr << "Warning: 'CalculatedQuantities' not found in configuration file." << std::endl;
	}
	for (YAML::Node calcQuantity : config["CalculatedQuantities"]) {
		std::vector<std::string> quantityBaseNames;
		if (not hasNodeKey(calcQuantity, "Name")) {
			std::cerr << "Could not convert a 'Name' of an entry in 'CalculatedQuantities'." << std::endl;
			return false;
		}
		if (calcQuantity["Name"].IsSequence()) {
			try {
				quantityBaseNames = calcQuantity["Name"].as<std::vector<std::string>>();
			} catch (const YAML::TypedBadConversion<std::vector<std::string>>& e) {
				std::cerr << "Could not convert YAML sequence to std::vector<std::string> when parsing 'CalculatedQuantities'' 'Name'." << std::endl;
				return false;
			} catch (const YAML::TypedBadConversion<std::string>& e) {
				std::cerr << "Could not entries in YAML sequence to std::string when parsing 'CalculatedQuantities'' 'Name'." << std::endl;
				return false;
			}
		} else {
			std::string baseName = antok::YAMLUtils::getString(calcQuantity["Name"]);
			if (baseName == "") {
				std::cerr << "Could not convert one of the 'CalculatedQuantities'' 'Name's to std::string." << std::endl;
				return false;
			}
			quantityBaseNames.push_back(baseName);
		}
		if (quantityBaseNames.empty()) {
			std::cerr << "Did not find a name to save calculated quantity to." << std::endl;
			return false;
		}

		if (not (hasNodeKey(calcQuantity, "Function") and hasNodeKey(calcQuantity["Function"], "Name"))) {
			std::cerr << "No Function or no function name for calculated quantity '" << quantityBaseNames[0] << "'." << std::endl;
			return false;
		}

		std::vector<int> indices;
		if (hasNodeKey(calcQuantity, "Indices")) {
			try {
				indices = calcQuantity["Indices"].as<std::vector<int>>();
			} catch (const YAML::TypedBadConversion<std::vector<int>>& e) {
				std::cerr << "Could not convert YAML sequence to std::vector<int> when parsing CalculatedQuantities' '" << quantityBaseNames[0] << "' 'Indices'." << std::endl;
				return false;
			} catch (const YAML::TypedBadConversion<int>& e) {
				std::cerr << "Could not convert entries in YAML sequence to int when parsing CalculatedQuantities' '" << quantityBaseNames[0] << "' 'Indices'." << std::endl;
				return false;
			}
		} else {
			indices.push_back(-1);
		}

		const YAML::Node&  function     = calcQuantity["Function"];
		const std::string& functionName = antok::YAMLUtils::getString(function["Name"]);
		for (const int index : indices) {
			std::vector<std::string> quantityNames;
			for (const std::string& baseName : quantityBaseNames) {
				if (index > 0) {
					std::stringstream qName;
					qName << baseName << index;
					quantityNames.push_back(qName.str());
				} else {
					quantityNames.push_back(baseName);
				}
			}

			antok::Function* antokFunctionPtr = nullptr;
			antokFunctionPtr = antok::Initializer::getFunction(function, quantityNames, index);
			if (antokFunctionPtr == nullptr) {
				std::cerr << "Error initializing function with name '" + functionName + "' which should calculate [";
				for (size_t i = 0; i < (quantityNames.size() - 1); ++i) {
					std::cerr << quantityNames[i] << ", ";
				}
				std::cerr << quantityNames[quantityNames.size() - 1] << "]." << std::endl;
				return false;
			}
			antok::Event& event = objectManager->getEvent();
			event._functions.push_back(antokFunctionPtr);

		}  // loop over indices
	}  // loop over CalculatedQuantities

	return true;
}


bool
antok::Initializer::initializePlotter()
{
	using antok::YAMLUtils::hasNodeKey;
	if (_config == nullptr) {
		std::cerr << "Trying to initialize Plotter without having read the config file first." << std::endl;
		return false;
	}
	const YAML::Node&     config        = *_config;
	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	if (objectManager->_plotter != nullptr) {
		std::cerr << "Event seems to be initialized already." << std::endl;
		return false;
	}
	objectManager->_plotter = antok::Plotter::instance();
	antok::Plotter& plotter = antok::ObjectManager::instance()->getPlotter();
	antok::plotUtils::GlobalPlotOptions plotOptions(config["GlobalPlotOptions"]);

	if (        (plotOptions.statisticsHistInName == "" or  plotOptions.statisticsHistOutName == "")
	    and not (plotOptions.statisticsHistInName == "" and plotOptions.statisticsHistOutName == ""))
	{
		std::cerr << "Something went wrong when parsing the 'GlobalPlotOptions'." << std::endl;
		return false;
	}

	if (plotOptions.statisticsHistInName != "") {
		TFile* inFile = objectManager->getInFile();
		TH1D* statsHistTemplate = dynamic_cast<TH1D*>(inFile->Get(plotOptions.statisticsHistInName.c_str()));
		if (statsHistTemplate == nullptr) {
			std::cerr << "Could not get the input 'StatisticsHistogram' from key '" << plotOptions.statisticsHistInName << "' "
			          << "in the input file '" << inFile->GetName() << "'." << std::endl;
			return false;
		}
		const antok::Cutter& cutter = objectManager->getCutter();
		TFile* outFile = objectManager->getOutFile();
		std::vector<antok::plotUtils::waterfallHistogramContainer> waterfallHists;
		for (const std::pair<std::string, std::vector<antok::Cut*>>& cutTrain : cutter._cutTrainsCutOrderMap) {
			const std::string& cutTrainName = cutTrain.first;
			outFile->cd(cutTrainName.c_str());
			TH1D* statsHist = (TH1D*)statsHistTemplate->Clone(plotOptions.statisticsHistOutName.c_str());
			if (statsHist == nullptr) {
				std::cerr << "Could not create the output 'StatisticsHistogram' for 'CutTrain' '" << cutTrainName << "' "
				          << "in key '" << cutTrainName + "/" + plotOptions.statisticsHistOutName << "' "
				          << "in the output file '" << outFile->GetName() << "'." << std::endl;
				return false;
			}
			objectManager->registerObjectToWrite(TDirectory::CurrentDirectory(), statsHist);
			const std::vector<antok::Cut*>& cuts = cutTrain.second;
			std::vector<std::pair<std::string, const bool*>> cutsAndResults;
			for (size_t i = 0; i < cuts.size(); ++i) {
				cutsAndResults.push_back(std::pair<std::string, const bool*>(cuts[i]->getLongName(), cutter.getCutResult(cuts[i])));
			}
			waterfallHists.push_back(antok::plotUtils::waterfallHistogramContainer(statsHist, cutsAndResults));
		}
		plotter._waterfallHistograms = waterfallHists;
		outFile->cd();
	}

	if (not hasNodeKey(config, "Plots")) {
		std::cerr << "Warning: 'Plots' not found in configuration file." << std::endl;
	}
	for (const YAML::Node& plot : config["Plots"]) {
		if (not hasNodeKey(plot, "Name")) {
			std::cerr << "'Name' not found for one of the 'Plots'." << std::endl;
			return false;
		}
		const std::string plotName = antok::YAMLUtils::getString(plot["Name"]);

		if (plotName.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_ ") != std::string::npos) {
			std::cerr << "Invalid character in name for 'Plot' '" << plotName << "'." << std::endl;
				return false;
		}
		if ((not(   (hasNodeKey(plot, "Variable")  and hasNodeKey(plot, "LowerBound")  and hasNodeKey(plot, "UpperBound"))
		         or (hasNodeKey(plot, "Variables") and hasNodeKey(plot, "LowerBounds") and hasNodeKey(plot, "UpperBounds"))))
		     or not  hasNodeKey(plot, "NBins")) {
			std::cerr << "Required variables not found for 'Plot' '" << plotName;
			std::cerr << "' ('Variable(s)'|'LowerBound(s)'|'UpperBound(s)'|'NBins')." << std::endl;
			return false;
		}

		antok::Plot* antokPlot = nullptr;
		if (hasNodeKey(plot, "Variables")) {
			if (plot["Variables"].IsSequence()) {
				if (not (    plot["LowerBounds"].IsSequence()
				         and plot["UpperBounds"].IsSequence()
				         and plot["NBins"].IsSequence())) {
					std::cerr << "'Variables', 'LowerBounds', 'UpperBounds' and 'NBins' all need to be sequences (in 'Plot' '" << plotName << "')." << std::endl;
					return false;
				}
				if (plot["Variables"].size() > 3) {
					std::cerr << "Cannot have 'Plot' '" << plotName << "' with more than 2 'Variables'" << std::endl;
					return false;
				}
				if (   (plot["Variables"].size() != plot["LowerBounds"].size())
				    or (plot["Variables"].size() != plot["UpperBounds"].size())
				    or (plot["Variables"].size() != plot["NBins"].size())) {
					std::cerr << "'Variables', 'LowerBounds', 'UpperBounds' and 'NBins' need to have the same number of entries (in 'Plot' '" << plotName << "')." << std::endl;
					return false;
				}
				if (plot["Variables"].size() == 2) {
					antokPlot = antok::generators::generate2DPlot(plot, plotOptions);
				} else if (plot["Variables"].size() == 1) {
					antokPlot = antok::generators::generate1DPlot(plot, plotOptions);
				} else if (plot["Variables"].size() == 3) {
					antokPlot = antok::generators::generate3DPlot(plot, plotOptions);
				} else {
					std::cerr << "Empty 'Variables' sequence found in 'Plot' '" << plotName << "'." << std::endl;
					return false;
				}
			} else {
				std::cerr << "'Variables' implies 2D/3D plot and needs to be a list (in 'Plot' '" << plotName << "')." << std::endl;
				return false;
			}
		} else if (hasNodeKey(plot, "Variable")) {
			antokPlot = antok::generators::generate1DPlot(plot, plotOptions);
		} else {
			assert(false);
		}

		if (antokPlot == nullptr) {
			std::cerr << "Could not generate 'Plot' '" << plotName << "'." << std::endl;
			return false;
		}
		plotter._plots.push_back(antokPlot);
	}

	return true;
}


template <typename T>
void
addToOutputBranch(TTree* const       outTree,
                  antok::Data&       data,
                  const std::string& variable_name)
{
	outTree->Branch(variable_name.c_str(), data.getAddr<T>(variable_name));
}


TTree*
createOutTree(TTree* const      inTree,
              const YAML::Node& config)
{
	using antok::YAMLUtils::hasNodeKey;
	TTree* outTree = nullptr;
	if (not hasNodeKey(config, "OutputTree")) {  // no definition for the output tree given -> write full input tree
		outTree = inTree->CloneTree(0);
	} else {  // own definition of output tree is given -> build output tree
		antok::Data& data = antok::ObjectManager::instance()->getData();
		outTree = new TTree(inTree->GetName(), inTree->GetTitle());

		for (const auto& variableNode : config["OutputTree"]) {
			const std::string variableName = antok::YAMLUtils::getString(variableNode);
			const std::string variableType = data.getType(variableName);
			if (variableName != "") {
				if      (variableType == "double")	            addToOutputBranch<double>             (outTree, data, variableName);
				else if (variableType == "int")	                addToOutputBranch<int>                (outTree, data, variableName);
				else if (variableType == "TVector3")            addToOutputBranch<TVector3>           (outTree, data, variableName);
				else if (variableType == "TLorentzVector")      addToOutputBranch<TLorentzVector>     (outTree, data, variableName);
				else if (variableType == "std::vector<double>") addToOutputBranch<std::vector<double>>(outTree, data, variableName);
				else if (variableType == "std::vector<int>")    addToOutputBranch<std::vector<int>>   (outTree, data, variableName);
				else {
					std::cerr << "Variable type '" << variableType << "' of variable '"
					          << variableName << "' not implemented for own output tree." << std::endl;
					return nullptr;
				}
			} else {
				std::cerr << "Variable '" << variableName << "' not found in Data's global map for the output tree." << std::endl;
				return nullptr;
			}
		}
	}

	return outTree;
}


antok::Function*
antok::Initializer::getFunction(const YAML::Node&               function,
                                const std::vector<std::string>& quantityNames,
                                const int                       index)
{
	const std::string& functionName = antok::YAMLUtils::getString(function["Name"]);
	if        (functionName == "abs") {
		return antok::generators::generateAbs                       (function, quantityNames, index);
	} else if (functionName == "log") {
	  return antok::generators::generateLog                         (function, quantityNames, index);
	} else if (functionName == "Sqrt") {
	  return antok::generators::generateSqrt                        (function, quantityNames, index);
	} else if (functionName == "convertIntToDouble") {
		return antok::generators::generateConvertIntToDouble        (function, quantityNames, index);
	} else if (functionName == "diff") {
		return antok::generators::generateDiff                      (function, quantityNames, index);
	} else if (functionName == "quotient") {
		return antok::generators::generateQuotient                  (function, quantityNames, index);
	} else if (functionName == "mul") {
		return antok::generators::generateMul                       (function, quantityNames, index);
	} else if (functionName == "energy") {
		return antok::generators::generateEnergy                    (function, quantityNames, index);
	} else if (functionName == "getBeamLorentzVector") {
		return antok::generators::generateGetBeamLorentzVector      (function, quantityNames, index);
	} else if (functionName == "getGradXGradY") {
		return antok::generators::generateGetGradXGradY             (function, quantityNames, index);
	} else if (functionName == "getLorentzVectorAttributes") {
		return antok::generators::generateGetLorentzVectorAttributes(function, quantityNames, index);
	} else if (functionName == "getLorentzVec") {
		return antok::generators::generateGetLorentzVec             (function, quantityNames, index);
	} else if (functionName == "getTs") {
		return antok::generators::generateGetTs                     (function, quantityNames, index);
	} else if (functionName == "getVector3") {
		return antok::generators::generateGetVector3                (function, quantityNames, index);
	} else if (functionName == "getVectorEntry") {
		return antok::generators::generateGetVectorEntry            (function, quantityNames, index);
	} else if (functionName == "mass") {
		return antok::generators::generateMass                      (function, quantityNames, index);
	} else if (functionName == "mass2") {
		return antok::generators::generateMass2                     (function, quantityNames, index);
	} else if (functionName == "radToDegree") {
		return antok::generators::generateRadToDegree               (function, quantityNames, index);
	} else if (functionName == "sum") {
		return antok::generators::generateSum                       (function, quantityNames, index);
	} else if (functionName == "sum2") {
		return antok::generators::generateSum2                      (function, quantityNames, index);
	} else if (functionName == "") {
		std::cerr << "Could not convert function name to std::string for calculated quantity '" << quantityNames[0] << "'." << std::endl;
		return nullptr;
	} else {
		antok::Function* usrFunctionPtr = antok::user::getUserFunction(function, quantityNames, index);
		if (usrFunctionPtr == nullptr) {
		    std::cerr << "Function type '" << functionName << "' did not initialize correctly or is not supported." << std::endl;
		}
		return usrFunctionPtr;
	}
}
