#include<function_generators.h>

#include<assert.h>

#include<TLorentzVector.h>

#include<data.hpp>
#include<functions.hpp>
#include<object_manager.h>

antok::Function* antok::generators::registerAbs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Arg", "double"));

	if(not functionArgumentHandler(args, function, index)) {
		std::cerr<<getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* argAddr = data.getDoubleAddr(args[0].first);

	if(not data.insertDouble(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::Abs(argAddr, data.getDoubleAddr(quantityName)));

};

antok::Function* antok::generators::registerDiff(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Minuend", "double"));
	args.push_back(std::pair<std::string, std::string>("Subtrahend", "double"));

	if(not functionArgumentHandler(args, function, index)) {
		std::cerr<<getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* arg1Addr = data.getDoubleAddr(args[0].first);
	double* arg2Addr = data.getDoubleAddr(args[1].first);

	if(not data.insertDouble(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::Diff(arg1Addr, arg2Addr, data.getDoubleAddr(quantityName)));

};

antok::Function* antok::generators::registerGetBeamLorentzVector(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;

	args.push_back(std::pair<std::string, std::string>("dX", "double"));
	args.push_back(std::pair<std::string, std::string>("dY", "double"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not functionArgumentHandler(args, function, index)) {
		std::cerr<<getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	double* dXaddr = data.getDoubleAddr(args[0].first);
	double* dYaddr = data.getDoubleAddr(args[1].first);
	TLorentzVector* xLorentzVecAddr = data.getLorentzVectorAddr(args[2].first);

	if(not data.insertLorentzVector(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::GetBeamLorentzVec(dXaddr, dYaddr, xLorentzVecAddr, data.getLorentzVectorAddr(quantityName)));

};

antok::Function* antok::generators::registerGetLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
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
		return 0;
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
			return 0;
		}
		mAddr = new double();
		(*mAddr) = function["M"].as<double>();
	}

	if(not functionArgumentHandler(args, function, index)) {
		std::cerr<<getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	xAddr = data.getDoubleAddr(args[0].first);
	yAddr = data.getDoubleAddr(args[1].first);
	zAddr = data.getDoubleAddr(args[2].first);
	if(pType) {
		mAddr = data.getDoubleAddr(args[3].first);
	}

	if(not data.insertLorentzVector(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::GetLorentzVec(xAddr, yAddr, zAddr, mAddr, data.getLorentzVectorAddr(quantityName), pType));

};

antok::Function* antok::generators::registerGetRpdPhi(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 2) {
		std::cerr<<"Need 2 names for function \"getRpdPhi\""<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::string method = getYAMLStringSafe(function["Method"]);
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("RPDProtonLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not functionArgumentHandler(args, function, index)) {
		std::cerr<<getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
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
		return 0;
	}

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insertDouble(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getDoubleAddr(quantityNames[i]));
	}

	return (new antok::functions::GetRpdPhi(beamLVAddr, RPDprotLVAddr, xLVAddr, quantityAddrs[0], quantityAddrs[1], methodSwitch));

};

antok::Function* antok::generators::registerGetTs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 3) {
		std::cerr<<"Need 3 names for function \"getTs\""<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;

	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if( not functionArgumentHandler(args, function, index)) {
		std::cerr<<getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* beamLVAddr = data.getLorentzVectorAddr(args[0].first);
	TLorentzVector* xLVAddr = data.getLorentzVectorAddr(args[1].first);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insertDouble(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getDoubleAddr(quantityNames[i]));
	}

	return (new antok::functions::GetTs(xLVAddr, beamLVAddr, quantityAddrs[0], quantityAddrs[1], quantityAddrs[2]));

};

antok::Function* antok::generators::registerMass(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));

	if(not functionArgumentHandler(args, function, index)) {
		std::cerr<<getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* vector = data.getLorentzVectorAddr(args[0].first);
	if(not data.insertDouble(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::Mass(vector, data.getDoubleAddr(quantityName)));

};

antok::Function* antok::generators::registerSum(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
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
		return 0;
	}

	if(function["Summands"]["Indices"] or function["Summands"]["Name"]) {
		if(not function["Summands"]["Indices"] and function["Summands"]["Name"]) {
			std::cerr<<"Either \"Summands\" or \"Name\" found in sum function, but not both (Variable: \""<<quantityName<<"\")."<<std::endl;
			return 0;
		}
		if(index > 0) {
			std::cerr<<"Cannot have sum over indices for every particle (Variable: \""<<quantityName<<"\")."<<std::endl;
			return 0;
		}
		std::vector<int> inner_indices;
		try {
			inner_indices = function["Summands"]["Indices"].as<std::vector<int> >();
		} catch (YAML::TypedBadConversion<std::vector<int> > e) {
			std::cerr<<"Could not convert YAML sequence to std::vector<int> when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
			return 0;
		} catch (YAML::TypedBadConversion<int> e) {
			std::cerr<<"Could not convert entries in YAML sequence to int when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
			return 0;
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
			return 0;
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
				return 0;
			}
			summandNames.push_back(std::pair<std::string, std::string>(variableName, typeName));
		}
	}
	if(not functionArgumentHandler(summandNames, function, index, true)) {
		std::cerr<<getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
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
				return 0;
			}

		}

	// And produce the function
	antok::Function* returnFunction = 0;
	if(typeName == "double") {
		if(not data.insertDouble(quantityName)) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return 0;
		}
		returnFunction = new antok::functions::Sum<double>(doubleInputAddrs, data.getDoubleAddr(quantityName));
	} else if(typeName == "int") {
		if(not data.insertInt(quantityName)) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return 0;
		}
		returnFunction = new antok::functions::Sum<int>(intInputAddrs, data.getIntAddr(quantityName));
	} else if(typeName == "Long64_t") {
		if(not data.insertLong64_t(quantityName)) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return 0;
		}
		returnFunction = new antok::functions::Sum<Long64_t>(long64_tInputAddrs, data.getLong64_tAddr(quantityName));
	} else if(typeName == "TLorentzVector") {
		if(not data.insertLorentzVector(quantityName)) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return 0;
		}
		returnFunction = new antok::functions::Sum<TLorentzVector>(lorentzVectorInputAddrs, data.getLorentzVectorAddr(quantityName));
	}
	assert(returnFunction != 0);

	return returnFunction;

};

antok::Function* antok::generators::registerSum2(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double*> doubleInputAddrs;
	std::vector<std::pair<std::string, std::string> > summandNames;

	// Summing over one variable with indices
	if(not function["Summands"]) {
		std::cerr<<"Either trying to sum different types or input variable not found when calculating \""<<quantityName<<"\"."<<std::endl;
		return 0;
	}

	if(function["Summands"]["Indices"] or function["Summands"]["Name"]) {
		if(not function["Summands"]["Indices"] and function["Summands"]["Name"]) {
			std::cerr<<"Either \"Summands\" or \"Name\" found in sum function, but not both (Variable: \""<<quantityName<<"\")."<<std::endl;
			return 0;
		}
		if(index > 0) {
			std::cerr<<"Cannot have sum over indices for every particle (Variable: \""<<quantityName<<"\")."<<std::endl;
			return 0;
		}
		std::vector<int> inner_indices;
		try {
			inner_indices = function["Summands"]["Indices"].as<std::vector<int> >();
		} catch (YAML::TypedBadConversion<std::vector<int> > e) {
			std::cerr<<"Could not convert YAML sequence to std::vector<int> when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
			return 0;
		} catch (YAML::TypedBadConversion<int> e) {
			std::cerr<<"Could not convert entries in YAML sequence to int when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
			return 0;
		}
		std::string summandBaseName = getYAMLStringSafe(function["Summands"]["Name"]);
		if(summandBaseName == "") {
			std::cerr<<"Could not convert \"Summands\"' \"Name\" to std::string when registering calculation of \""<<quantityName<<"\"."<<std::endl;
			return 0;
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
				return 0;
			}
			summandNames.push_back(std::pair<std::string, std::string>(summandName, "double"));
		}
	}

	if(not functionArgumentHandler(summandNames, function, index, true)) {
		std::cerr<<getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	// Now do type checking and get all the addresses
	for(unsigned int summandNames_i = 0; summandNames_i < summandNames.size(); ++summandNames_i) {

			std::string variableName = summandNames[summandNames_i].first;

			double* addr = data.getDoubleAddr(variableName);
			doubleInputAddrs.push_back(addr);

		}

	// And produce the function
	if(not data.insertDouble(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::Sum2(doubleInputAddrs, data.getDoubleAddr(quantityName)));

};

bool antok::generators::functionArgumentHandler(std::vector<std::pair<std::string, std::string> >& args, const YAML::Node& function, int index, bool argStringsAlreadyValues)
{

	antok::Data& data = antok::ObjectManager::instance()->getData();
	for(unsigned int i = 0; i < args.size(); ++i) {
		std::string& argName = args[i].first;
		if(not argStringsAlreadyValues) {
			if(not function[argName]) {
				std::cerr<<"Argument \""<<argName<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
				return 0;
			}
			argName = getYAMLStringSafe(function[argName]);
			if(argName == "") {
				std::cerr<<"Could not convert one of the arguments to std::string in function \""<<function["Name"]<<"\"."<<std::endl;
				return 0;
			}
		}
		if(index > 0) {
			std::stringstream strStr;
			strStr<<argName<<index;
			argName = strStr.str();
		}
		if(data.global_map.count(argName) < 1) {
			std::cerr<<"Argument \""<<argName<<"\" not found in Data's global map."<<std::endl;
			return 0;
		}
		if(data.global_map[argName] != args[i].second) {
			std::cerr<<"Argument \""<<argName<<"\" has type \""<<data.global_map[argName]<<"\", expected \""<<args[i].second<<"\"."<<std::endl;
			return 0;
		}
	}

	return true;

};

std::string antok::generators::getYAMLStringSafe(const YAML::Node& node) {
	try{
		return node.as<std::string>();
	} catch(YAML::TypedBadConversion<std::string> e) {
		return "";
	}
}

std::string antok::generators::getFunctionArgumentHandlerErrorMsg(std::vector<std::string> quantityNames) {
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

