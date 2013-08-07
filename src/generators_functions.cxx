#include<generators_functions.h>

#include<assert.h>

#include<TLorentzVector.h>

#include<data.h>
#include<functions.hpp>
#include<initializer.h>
#include<object_manager.h>
#include<yaml_utils.hpp>

namespace {

	bool __functionArgumentHandler(std::vector<std::pair<std::string, std::string> >& args, const YAML::Node& function, int index, bool argStringsAlreadyValues = false)
	{

		using antok::YAMLUtils::hasNodeKey;

		antok::Data& data = antok::ObjectManager::instance()->getData();
		for(unsigned int i = 0; i < args.size(); ++i) {
			std::string& argName = args[i].first;
			if(not argStringsAlreadyValues) {
				if(not hasNodeKey(function, argName)) {
					std::cerr<<"Argument \""<<argName<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
					return false;
				}
				argName = antok::YAMLUtils::getString(function[argName]);
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
			std::string type = data.getType(argName);
			if(type == "") {
				std::cerr<<"Argument \""<<argName<<"\" not found in Data's global map."<<std::endl;
				return false;
			}
			if(type != args[i].second) {
				std::cerr<<"Argument \""<<argName<<"\" has type \""<<type<<"\", expected \""<<args[i].second<<"\"."<<std::endl;
				return false;
			}
		}

		return true;

	};

	std::string __getFunctionArgumentHandlerErrorMsg(std::vector<std::string> quantityNames) {
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

	template<typename T>
	antok::Function* __getSumFunction(std::vector<std::pair<std::string, std::string> >& summandNames, std::string quantityName) {

		std::vector<T*> inputAddrs;
		antok::Data& data = antok::ObjectManager::instance()->getData();

		// Now do type checking and get all the addresses
		for(unsigned int summandNames_i = 0; summandNames_i < summandNames.size(); ++summandNames_i) {

			std::string variableName = summandNames[summandNames_i].first;
			inputAddrs.push_back(data.getAddr<T>(variableName));

		}

		// And produce the function
		if(not data.insert<T>(quantityName)) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityName);
			return 0;
		}
		return (new antok::functions::Sum<T>(inputAddrs, data.getAddr<T>(quantityName)));

	};

	std::vector<std::pair<std::string, std::string> >* __getSummandNames(const YAML::Node& function, std::string& quantityName, int index) {

		using antok::YAMLUtils::hasNodeKey;

		std::string typeName = "notInitialized";
		std::vector<std::pair<std::string, std::string> >* summandNames = new std::vector<std::pair<std::string, std::string> >();

		// Summing over one variable with indices
		if(not hasNodeKey(function, "Summands")) {
			std::cerr<<"Argument \"Summands\" is missing in configuration file for variable \""<<quantityName<<"\"."<<std::endl;
			return 0;
		}

		antok::Data& data = antok::ObjectManager::instance()->getData();

		if(hasNodeKey(function["Summands"], "Indices") or hasNodeKey(function["Summands"], "Name")) {
			if(not hasNodeKey(function["Summands"], "Indices") and hasNodeKey(function["Summands"], "Name")) {
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
			} catch (const YAML::TypedBadConversion<std::vector<int> >& e) {
				std::cerr<<"Could not convert YAML sequence to std::vector<int> when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
				return 0;
			} catch (const YAML::TypedBadConversion<int>& e) {
				std::cerr<<"Could not convert entries in YAML sequence to int when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
				return 0;
			}
			typeName = antok::YAMLUtils::getString(function["Summands"]["Name"]);
			if(typeName == "") {
				std::cerr<<"Could not convert \"Summands\"' \"Name\" to std::string when registering calculation of \""<<quantityName<<"\"."<<std::endl;
			}
			std::string summandBaseName = typeName;
			std::stringstream strStr;
			strStr<<typeName<<inner_indices[0];
			typeName = data.getType(strStr.str());
			for(unsigned int inner_indices_i = 0; inner_indices_i < inner_indices.size(); ++inner_indices_i) {
				int inner_index = inner_indices[inner_indices_i];
				std::stringstream strStr;
				strStr<<summandBaseName<<inner_index;
				summandNames->push_back(std::pair<std::string, std::string>(strStr.str(), typeName));
			}
			// Summing over list of variable names
		} else {
			typeName = antok::YAMLUtils::getString(*function["Summands"].begin());
			if(typeName == "") {
				std::cerr<<"Could not convert one of the \"Summands\" to std::string when registering calculation of \""<<quantityName<<"\"."<<std::endl;
				return 0;
			}
			if(index > 0) {
				std::stringstream strStr;
				strStr<<typeName<<index;
				typeName = strStr.str();
			}
			typeName = data.getType(typeName);
			for(YAML::const_iterator summand_it = function["Summands"].begin(); summand_it != function["Summands"].end(); ++summand_it) {
				std::string variableName = antok::YAMLUtils::getString(*summand_it);
				if(variableName == "") {
					std::cerr<<"Could not convert one of the \"Summands\" to std::string when registering calculation of \""<<quantityName<<"\"."<<std::endl;
					return 0;
				}
				summandNames->push_back(std::pair<std::string, std::string>(variableName, typeName));
			}
		}
		return summandNames;

	};

}

antok::Function* antok::generators::generateAbs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Arg", "double"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* argAddr = data.getAddr<double>(args[0].first);

	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::Abs(argAddr, data.getAddr<double>(quantityName)));

};

antok::Function* antok::generators::generateConvertIntToDouble(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Int", "int"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	int* argAddr = data.getAddr<int>(args[0].first);

	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::ConvertIntToDouble(argAddr, data.getAddr<double>(quantityName)));

};

antok::Function* antok::generators::generateDiff(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Minuend", "double"));
	args.push_back(std::pair<std::string, std::string>("Subtrahend", "double"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* arg1Addr = data.getAddr<double>(args[0].first);
	double* arg2Addr = data.getAddr<double>(args[1].first);

	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::Diff(arg1Addr, arg2Addr, data.getAddr<double>(quantityName)));

};

antok::Function* antok::generators::generateEnergy(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	TLorentzVector* arg1Addr = data.getAddr<TLorentzVector>(args[0].first);

	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::Energy(arg1Addr, data.getAddr<double>(quantityName)));

};

antok::Function* antok::generators::generateGetBeamLorentzVector(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
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

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	double* dXaddr = data.getAddr<double>(args[0].first);
	double* dYaddr = data.getAddr<double>(args[1].first);
	TLorentzVector* xLorentzVecAddr = data.getAddr<TLorentzVector>(args[2].first);

	if(not data.insert<TLorentzVector>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::GetBeamLorentzVec(dXaddr, dYaddr, xLorentzVecAddr, data.getAddr<TLorentzVector>(quantityName)));

};

antok::Function* antok::generators::generateGetGradXGradY(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 2) {
		std::cerr<<"Need 3 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* lorentzVectorAddr = data.getAddr<TLorentzVector>(args[0].first);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return (new antok::functions::GetGradXGradY(lorentzVectorAddr, quantityAddrs[0], quantityAddrs[1]));

};

antok::Function* antok::generators::generateGetLorentzVectorAttributes(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 5) {
		std::cerr<<"Need 5 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* lorentzVectorAddr = data.getAddr<TLorentzVector>(args[0].first);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return (new antok::functions::GetLorentzVectorAttributes(lorentzVectorAddr,
	                                                         quantityAddrs[0],
	                                                         quantityAddrs[1],
	                                                         quantityAddrs[2],
	                                                         quantityAddrs[3],
	                                                         quantityAddrs[4]));

};

antok::Function* antok::generators::generateGetLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
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
		} catch(const YAML::TypedBadConversion<double>& e) {
			std::cerr<<"Argument \"M\" in function \"mass\" should be of type double (variable \""<<quantityName<<"\")."<<std::endl;
			return 0;
		}
		mAddr = new double();
		(*mAddr) = function["M"].as<double>();
	}

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	xAddr = data.getAddr<double>(args[0].first);
	yAddr = data.getAddr<double>(args[1].first);
	zAddr = data.getAddr<double>(args[2].first);
	if(pType) {
		mAddr = data.getAddr<double>(args[3].first);
	}

	if(not data.insert<TLorentzVector>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::GetLorentzVec(xAddr, yAddr, zAddr, mAddr, data.getAddr<TLorentzVector>(quantityName), pType));

};

antok::Function* antok::generators::generateGetRpdExpectedHitsParameters(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {

	using antok::YAMLUtils::hasNodeKey;

	if(quantityNames.size() != 4) {
		std::cerr<<"Need 4 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("Vertex", "TVector3"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* beamLorentzVecAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* XLorentzVecAddr = data.getAddr<TLorentzVector>(args[1].first);
	TVector3* vertexAddr = data.getAddr<TVector3>(args[2].first);

	std::vector<std::string> possiblyConstArgs;
	possiblyConstArgs.push_back("XOffset");
	possiblyConstArgs.push_back("YOffset");
	possiblyConstArgs.push_back("XAngle");
	possiblyConstArgs.push_back("YAngle");

	for(unsigned int i = 0; i < possiblyConstArgs.size(); ++i) {
		if(not hasNodeKey(function, possiblyConstArgs[i])) {
			std::cerr<<"Argument \""<<possiblyConstArgs[i]<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
	}

	double* xOffsetAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[0]]);
	double* yOffsetAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[1]]);
	double* xAngleAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[2]]);
	double* yAngleAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[3]]);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return (new antok::functions::GetRPDExpectedHitsParameters(beamLorentzVecAddr,
	                                                           XLorentzVecAddr,
	                                                           vertexAddr,
	                                                           xOffsetAddr,
	                                                           yOffsetAddr,
	                                                           xAngleAddr,
	                                                           yAngleAddr,
	                                                           quantityAddrs[0],
	                                                           quantityAddrs[1],
	                                                           quantityAddrs[2],
	                                                           quantityAddrs[3]));

}

antok::Function* antok::generators::generateGetRpdPhi(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 2) {
		std::cerr<<"Need 2 names for function \"getRpdPhi\""<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::string method = antok::YAMLUtils::getString(function["Method"]);
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("RPDProtonLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* beamLVAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* RPDprotLVAddr = data.getAddr<TLorentzVector>(args[1].first);
	TLorentzVector* xLVAddr = data.getAddr<TLorentzVector>(args[2].first);

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
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return (new antok::functions::GetRpdPhi(beamLVAddr, RPDprotLVAddr, xLVAddr, quantityAddrs[0], quantityAddrs[1], methodSwitch));

};

antok::Function* antok::generators::generateGetTs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 3) {
		std::cerr<<"Need 3 names for function \"getTs\""<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;

	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* beamLVAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* xLVAddr = data.getAddr<TLorentzVector>(args[1].first);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return (new antok::functions::GetTs(xLVAddr, beamLVAddr, quantityAddrs[0], quantityAddrs[1], quantityAddrs[2]));

};

antok::Function* antok::generators::generateGetVector3(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	const std::string& quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("X", "double"));
	args.push_back(std::pair<std::string, std::string>("Y", "double"));
	args.push_back(std::pair<std::string, std::string>("Z", "double"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	double* xAddr = data.getAddr<double>(args[0].first);
	double* yAddr = data.getAddr<double>(args[1].first);
	double* zAddr = data.getAddr<double>(args[2].first);

	if(not data.insert<TVector3>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::GetTVector3(xAddr, yAddr, zAddr, data.getAddr<TVector3>(quantityName)));

}

antok::Function* antok::generators::generateMass(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* vector = data.getAddr<TLorentzVector>(args[0].first);
	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::Mass(vector, data.getAddr<double>(quantityName)));

};

antok::Function* antok::generators::generateRadToDegree(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance() ->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Angle", "double"));

	if(not __functionArgumentHandler(args, function, index)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	double* angle = data.getAddr<double>(args[0].first);
	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::RadToDegree(angle, data.getAddr<double>(quantityName)));

};

antok::Function* antok::generators::generateSum(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> >* summandNamesPtr = __getSummandNames(function, quantityName, index);
	if(summandNamesPtr == 0) {
		std::cerr<<"Could not generate summands for function \"Sum\" when trying to register calculation of \""<<quantityName<<"\"."<<std::endl;
		return 0;
	}
	std::vector<std::pair<std::string, std::string> >& summandNames = (*summandNamesPtr);

	if(not __functionArgumentHandler(summandNames, function, index, true)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	std::string typeName = summandNames[0].second;

	antok::Function* antokFunction = 0;
	if(typeName == "double") {
		antokFunction = __getSumFunction<double>(summandNames, quantityName);
	} else if (typeName == "int") {
		antokFunction = __getSumFunction<int>(summandNames, quantityName);
	} else if (typeName == "Long64_t") {
		antokFunction = __getSumFunction<Long64_t>(summandNames, quantityName);
	} else if (typeName == "TLorentzVector") {
		antokFunction = __getSumFunction<TLorentzVector>(summandNames, quantityName);
	} else {
		std::cerr<<"Type \""<<typeName<<"\" not supported by \"sum\" (registering calculation of \""<<quantityName<<"\")."<<std::endl;
		return 0;
	}

	delete summandNamesPtr;

	return (antokFunction);

};

antok::Function* antok::generators::generateSum2(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double*> doubleInputAddrs;

	std::vector<std::pair<std::string, std::string> >* summandNamesPtr = __getSummandNames(function, quantityName, index);
	if(summandNamesPtr == 0) {
		std::cerr<<"Could not generate summands for function \"Sum\" when trying to register calculation of \""<<quantityName<<"\"."<<std::endl;
		return 0;
	}
	std::vector<std::pair<std::string, std::string> >& summandNames = (*summandNamesPtr);

	if(not __functionArgumentHandler(summandNames, function, index, true)) {
		std::cerr<<__getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	// Now do type checking and get all the addresses
	for(unsigned int summandNames_i = 0; summandNames_i < summandNames.size(); ++summandNames_i) {

			std::string variableName = summandNames[summandNames_i].first;

			double* addr = data.getAddr<double>(variableName);
			doubleInputAddrs.push_back(addr);

		}

	// And produce the function
	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	delete summandNamesPtr;

	return (new antok::functions::Sum2(doubleInputAddrs, data.getAddr<double>(quantityName)));

};

