#include "TLorentzVector.h"
#include "TVector3.h"

#include "constants.h"
#include "data.h"
#include "functions.hpp"
#include "generators_functions.h"
#include "cdreis.h"
#include "jbeckers.h"
#include "jbeckers_functions.hpp"
#include "object_manager.h"
#include "yaml_utils.hpp"

using antok::generators::functionArgumentHandler;
using antok::generators::functionArgumentHandlerConst;
using antok::generators::getFunctionArgumentHandlerErrorMsg;
using antok::generators::nmbArgsIsExactly;
using antok::YAMLUtils::hasNodeKey;

// type aliases to save some typing
template <class T>
using vecPairString = std::vector<std::pair<std::string, T>>;
template <class T>
using mapStringT = std::map<std::string, T>;


antok::Function*
antok::user::jbeckers::getUserFunction(const YAML::Node&               function,
                                       const std::vector<std::string>& quantityNames,
                                       const int                       index)
{
	const std::string& functionName = antok::YAMLUtils::getString(function["Name"]);
	if (functionName == "scale") {
		return antok::user::jbeckers::generateScale                     (function, quantityNames, index);
	} else if (functionName == "increase") {
		return antok::user::jbeckers::generateIncrease                  (function, quantityNames, index);
	} else if (functionName == "CalculateCollinearityAngle") {
		return antok::user::jbeckers::generateCalculateCollinearityAngle(function, quantityNames, index);
	} else if (functionName == "CompareIndices") {
		return antok::user::jbeckers::generateCompareIndices            (function, quantityNames, index);
	} else if (functionName == "PrintEvent") {
		return antok::user::jbeckers::generatePrintEvent                (function, quantityNames, index);
	}
	return nullptr;
}


namespace {

	template <typename T>
	antok::Function*
	__generateScaleHelper(antok::Data&                      data,
	                      const vecPairString<std::string>& variableArgs,
	                      const mapStringT<double>&         constArgs,
	                      const std::vector<std::string>&   quantityNames,
	                      const std::string&                quantityName )
	{
		// Register output variables
		if (not data.insert<T>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}

		//std::cout << "Helper has const argument value of " << constArgs.at("Scale") << std::endl;

		return new antok::user::jbeckers::functions::Scale<T>(constArgs.at("Scale"),                    // scaleFactor
		                                                      *data.getAddr<T>(variableArgs[0].first),  // value
		                                                      *data.getAddr<T>(quantityName));          // result
	}

}  // anonymous namespace


antok::Function*
antok::user::jbeckers::generateScale(const YAML::Node&               function,
                                     const std::vector<std::string>& quantityNames,
                                     const int                       index)
{
	if (not antok::generators::nmbArgsIsExactly(function, quantityNames.size(), 1)) { // number of outputs args
		return nullptr;
	}

	// Get input variables
	//TODO improve naming of YAML nodes
	const std::string arg1Name = "Scale";  // constant scale factor, must be double
	mapStringT<double> constArgs = {{arg1Name, 0.0}};
	if (not antok::generators::functionArgumentHandlerConst<double>(constArgs, function, true)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	const std::string arg2Name     = "ScaleOn";  // variable to scale
	const std::string arg2TypeName = antok::generators::getTypeOfArg(function, index, arg2Name);
	vecPairString<std::string> variableArgs = {{arg2Name, arg2TypeName}};
	if (not antok::generators::functionArgumentHandler(variableArgs, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	//std::cout << "Found const argument value of " << constArgs[arg1Name] << std::endl;

	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if        (arg2TypeName == "double") {
		return __generateScaleHelper<double>             (data, variableArgs, constArgs, quantityNames, quantityName);
	// } else if (arg2TypeName == "int") {
	// 	return __generateScaleHelper<int>                (data, variableArgs, constArgs, quantityNames, quantityName);
	} else if (arg2TypeName == "TVector3") {
		return __generateScaleHelper<TVector3>           (data, variableArgs, constArgs, quantityNames, quantityName);
	} else if (arg2TypeName == "TLorentzVector") {
		return __generateScaleHelper<TLorentzVector>     (data, variableArgs, constArgs, quantityNames, quantityName);
	} else if (arg2TypeName == "std::vector<double>") {
		return __generateScaleHelper<std::vector<double>>(data, variableArgs, constArgs, quantityNames, quantityName);
	} else {
		std::cerr << "'" << function["Name"] << "' is not (yet) implemented for type '" << arg2TypeName << "' of "
		          << "variable '" << variableArgs[0].first << "'." << std::endl;
		return nullptr;
	}
}


namespace {

	template <typename T>
	antok::Function*
	__generateIncreaseHelper(antok::Data&                      data,
	                      const vecPairString<std::string>& variableArgs,
	                      const mapStringT<double>&         constArgs,
	                      const std::vector<std::string>&   quantityNames,
	                      const std::string&                quantityName )
	{
		// Register output variables
		if (not data.insert<T>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}

		//std::cout << "Helper has const argument value of " << constArgs.at("Increment") << std::endl;

		return new antok::user::jbeckers::functions::Increase<T>(constArgs.at("Increment"),                // increment
		                                                        *data.getAddr<T>(variableArgs[0].first),  // value
		                                                        *data.getAddr<T>(quantityName));          // result
	}

}  // anonymous namespace


antok::Function*
antok::user::jbeckers::generateIncrease(const YAML::Node&               function,
                                        const std::vector<std::string>& quantityNames,
                                        const int                       index)
{
	if (not antok::generators::nmbArgsIsExactly(function, quantityNames.size(), 1)) { // number of outputs args
		return nullptr;
	}

	// Get input variables
	//TODO improve naming of YAML nodes
	const std::string arg1Name = "Increment";  // constant increment, must be double
	mapStringT<double> constArgs = {{arg1Name, 0.0}};
	if (not antok::generators::functionArgumentHandlerConst<double>(constArgs, function, true)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	const std::string arg2Name     = "IncreaseOn";  // variable to increase
	const std::string arg2TypeName = antok::generators::getTypeOfArg(function, index, arg2Name);
	vecPairString<std::string> variableArgs = {{arg2Name, arg2TypeName}};
	if (not antok::generators::functionArgumentHandler(variableArgs, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	//std::cout << "Found const argument value of " << constArgs[arg1Name] << std::endl;

	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if        (arg2TypeName == "double") {
		return __generateIncreaseHelper<double>             (data, variableArgs, constArgs, quantityNames, quantityName);
	// } else if (arg2TypeName == "int") {
	// 	return __generateScaleHelper<int>                (data, variableArgs, constArgs, quantityNames, quantityName);
	} else if (arg2TypeName == "std::vector<double>") {
		return __generateIncreaseHelper<std::vector<double>>(data, variableArgs, constArgs, quantityNames, quantityName);
	} else {
		std::cerr << "'" << function["Name"] << "' is not (yet) implemented for type '" << arg2TypeName << "' of "
		          << "variable '" << variableArgs[0].first << "'." << std::endl;
		return nullptr;
	}
}

antok::Function*
antok::user::jbeckers::generateCalculateCollinearityAngle(const YAML::Node&               function,
                                                          const std::vector<std::string>& quantityNames,
                                                          const int                       index)
{
	if (not antok::generators::nmbArgsIsExactly(function, quantityNames.size(), 1)) { // number of outputs args
		return nullptr;
	}

	const std::string arg1Name     = "PrimVector";
	const std::string arg2Name     = "SecVector";
	const std::string arg3Name     = "V0LV";
	const std::string arg1TypeName = antok::generators::getTypeOfArg(function, index, arg1Name);
	const std::string arg2TypeName = antok::generators::getTypeOfArg(function, index, arg2Name);
	const std::string arg3TypeName = antok::generators::getTypeOfArg(function, index, arg3Name);
	vecPairString<std::string> variableArgs = {{arg1Name, arg1TypeName},{arg2Name, arg2TypeName},{arg3Name, arg3TypeName}};
	if (not functionArgumentHandler(variableArgs, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();
	// Register output variables
	const std::string& quantityName = quantityNames[0];
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	
	return new antok::user::jbeckers::functions::CalculateCollinearityAngle(*data.getAddr<TVector3>       (variableArgs[0].first),
																            *data.getAddr<TVector3>       (variableArgs[1].first),
																            *data.getAddr<TLorentzVector> (variableArgs[2].first),
																            *data.getAddr<double>         (quantityNames[0]));

}


antok::Function*
antok::user::jbeckers::generateCompareIndices(const YAML::Node&               function,
                                              const std::vector<std::string>& quantityNames,
                                              const int                       index)
{
	if (not antok::generators::nmbArgsIsExactly(function, quantityNames.size(), 1)) { // number of outputs args
		return nullptr;
	}

	const std::string arg1Name     = "IndicesVector";
	const std::string arg2Name     = "IndicesVectorVector";
	const std::string arg1TypeName = antok::generators::getTypeOfArg(function, index, arg1Name);
	const std::string arg2TypeName = antok::generators::getTypeOfArg(function, index, arg2Name);
	vecPairString<std::string> variableArgs = {{arg1Name, arg1TypeName},{arg2Name, arg2TypeName}};
	if (not functionArgumentHandler(variableArgs, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();
	// Register output variables
	const std::string& quantityName = quantityNames[0];
	if (not data.insert<int>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	
	return new antok::user::jbeckers::functions::CompareIndices(*data.getAddr<std::vector<int>>             (variableArgs[0].first),
																*data.getAddr<std::vector<std::vector<int>>>(variableArgs[1].first),
																*data.getAddr<int>                          (quantityNames[0]));

}


antok::Function*
antok::user::jbeckers::generatePrintEvent(const YAML::Node&               function,
                                          const std::vector<std::string>& quantityNames,
                                          const int                       index)
{
	if (not antok::generators::nmbArgsIsExactly(function, quantityNames.size(), 1)) { // number of outputs args
		return nullptr;
	}

	const std::string arg1Name     = "RunNumber";
	const std::string arg2Name     = "SpillNumber";
	const std::string arg3Name     = "EventNumber";
	const std::string arg1TypeName = antok::generators::getTypeOfArg(function, index, arg1Name);
	const std::string arg2TypeName = antok::generators::getTypeOfArg(function, index, arg2Name);
	const std::string arg3TypeName = antok::generators::getTypeOfArg(function, index, arg3Name);
	vecPairString<std::string> variableArgs = {{arg1Name, arg1TypeName},{arg2Name, arg2TypeName},{arg3Name, arg3TypeName}};
	if (not functionArgumentHandler(variableArgs, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();
	
	return new antok::user::jbeckers::functions::PrintEvent(*data.getAddr<int>(variableArgs[0].first),*data.getAddr<int>(variableArgs[1].first),*data.getAddr<int>(variableArgs[2].first));

}