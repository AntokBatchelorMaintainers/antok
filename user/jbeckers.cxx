#include "TLorentzVector.h"
#include "TVector3.h"

#include "constants.h"
#include "data.h"
#include "functions.hpp"
#include "generators_functions.h"
#include "jbeckers.h"
#include "jbeckers_functions.hpp"
#include "object_manager.h"
#include "yaml_utils.hpp"


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
		return antok::user::jbeckers::generateScale(function, quantityNames, index);
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
