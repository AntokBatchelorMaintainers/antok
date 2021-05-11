#include <constants.h>
#include <data.h>
#include <object_manager.h>
#include <functions.hpp>
#include <generators_functions.h>
#include <yaml_utils.hpp>
#include "TVector3.h"
#include "TLorentzVector.h"

#include <jbeckers.h>
#include <jbeckers_functions.hpp>

// type aliases to save some typing
template <class T>
using vecPairString = std::vector<std::pair<std::string, T>>;
template <class T>
using mapStringTo = std::map<std::string, T>;


antok::Function* antok::user::jbeckers::getUserFunction(const YAML::Node& function,
                                                          const std::vector<std::string>& quantityNames,
                                                          int index)
{
	std::string functionName = antok::YAMLUtils::getString(function["Name"]);
	antok::Function* antokFunctionPtr = 0;
	if(functionName == "scale") {
		antokFunctionPtr = jbeckers::generateScale(function, quantityNames, index);
	}
	return antokFunctionPtr;
}


namespace {

	template <typename T>
	antok::Function*
	__generateScaleHelper(antok::Data&                      data,
	                      const vecPairString<std::string>& arg_non_const,
                          const mapStringTo<double>&        arg_const,
	                      const std::vector<std::string>&   quantityNames,
	                      const std::string&                quantityName )
	{
		// Register output variables
		if (not data.insert<T>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}

		//std::cout << "Helper has const argument value of " << arg_const.at("Scale") << std::endl;

		return new antok::user::jbeckers::functions::Scale<T>(arg_const.at("Scale"), // scale
                                                              *data.getAddr<T>(arg_non_const[0].first), // scaleOn
                                                              *data.getAddr<T>(quantityName)); // output
	}

}  // anonymous namespace


antok::Function*
antok::user::jbeckers::generateScale(const YAML::Node&             function,
                                     const std::vector<std::string>& quantityNames,
                                     const int                       index)
{
	if (not antok::generators::nmbArgsIsExactly(function, quantityNames.size(), 1)) { // number of outputs args
		return nullptr;
	}

	// Get input variables
	const YAML::Node& functionName = function["Name"];
	const std::string arg1Name     = "Scale";    // constant scale factor, must be double
	const std::string arg2Name     = "ScaleOn";  // variable to scale
	const std::string arg2TypeName = antok::generators::getTypeOfArg(function, index, arg2Name);

	vecPairString<std::string> arg_non_const
		= {{arg2Name, arg2TypeName}};

	if (not antok::generators::functionArgumentHandler(arg_non_const, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

    mapStringTo<double> arg_const
		= {{arg1Name, 0.0}};

    if (not antok::generators::functionArgumentHandlerConst<double>(arg_const, function, true)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	//std::cout << "Found const argument value of " << arg_const[arg1Name] << std::endl;

	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if        (arg2TypeName == "double") {
		return __generateScaleHelper<double>             (data, arg_non_const, arg_const, quantityNames, quantityName);
    } else if (arg2TypeName == "TVector3") {
		return __generateScaleHelper<TVector3>           (data, arg_non_const, arg_const, quantityNames, quantityName);
    } else if (arg2TypeName == "TLorentzVector") {
		return __generateScaleHelper<TLorentzVector>     (data, arg_non_const, arg_const, quantityNames, quantityName);
	//} else if (arg2TypeName == "int") {
	//	return __generateScaleHelper<int>                (data, arg_non_const, arg_const, quantityNames, quantityName);
	} else if (arg2TypeName == "std::vector<double>") {
		return __generateScaleHelper<std::vector<double>>(data, arg_non_const, arg_const, quantityNames, quantityName);
	} else {
		std::cerr << "'" << functionName << "' is not (yet) implemented for type '" << arg2TypeName << "' of "
		          << "variable '" << arg_non_const[0].first << "'." << std::endl;
		return nullptr;
	}
}
