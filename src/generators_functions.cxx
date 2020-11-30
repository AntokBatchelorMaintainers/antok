#include <cassert>

#include "constants.h"
#include "data.h"
#include "functions.hpp"
#include "generators_functions.h"
#include "initializer.h"
#include "object_manager.h"
#include "yaml_utils.hpp"

using antok::YAMLUtils::hasNodeKey;

// type aliases to save some typing
template <class T>
using vecPairString = std::vector<std::pair<std::string, T>>;
template <class T>
using mapStringTo = std::map<std::string, T>;


std::string
antok::generators::mergeNameIndex(const std::string& name,
	                                const int          index)
{
	if (index > 0) {
		std::stringstream strStr;
		strStr << name << index;
		return strStr.str();
	}
	return name;
}


bool
antok::generators::nmbArgsIsExactly(const YAML::Node& function,
                                    const size_t&     actualNmb,
                                    const size_t&     requiredNmb)
{
	if (actualNmb == requiredNmb) {
		return true;
	}
	std::cerr << "Need " << requiredNmb << " argument name(s) instead of " << actualNmb << " name(s) "
	          << "for function '" << function["Name"] << "'." << std::endl;
	return false;
}


std::string
antok::generators::getTypeOfArg(const YAML::Node&  function,
                                const int          index,
                                const std::string& argName)
{
	const YAML::Node& functionName = function["Name"];
	antok::Data&      data         = antok::ObjectManager::instance()->getData();
	if (hasNodeKey(function, argName)) {
		return data.getType(antok::generators::mergeNameIndex(antok::YAMLUtils::getString(function[argName]), index));
	} else {
		std::cerr << "Argument '" << argName << "' not found (required for function '" << functionName << "')." << std::endl;
		return "";
	}
}


bool
antok::generators::functionArgumentHandler(vecPairString<std::string>& args,
                                           const YAML::Node&           function,
                                           const int                   index,
                                           const bool                  argStringsAlreadyValues)
{
	const YAML::Node& functionName = function["Name"];
	for (size_t i = 0; i < args.size(); ++i) {
		std::string& argName = args[i].first;
		if (not argStringsAlreadyValues) {
			if (not hasNodeKey(function, argName)) {
				std::cerr << "Argument '" << argName << "' not found (required for function '" << functionName << "')." << std::endl;
				return false;
			}
			argName = antok::YAMLUtils::getString(function[argName]);
			if (argName == "") {
				std::cerr << "Could not convert one of the arguments to std::string in function '" << functionName << "'." << std::endl;
				return false;
			}
		}
		argName = antok::generators::mergeNameIndex(argName, index);
		const antok::Data& data = antok::ObjectManager::instance()->getData();
		const std::string  type = data.getType(argName);
		if (type == "") {
			std::cerr << "Argument '" << argName << "' not found in Data's global map." << std::endl;
			return false;
		}
		if (type != args[i].second) {
			std::cerr << "Argument '" << argName << "' has type '" << type << "', expected '" << args[i].second << "'." << std::endl;
			return false;
		}
	}

	return true;
}


/**
 * Sets values in given args vector to the constants found in YAML node
 * @param args Map node name to node value (values will be set in this function)
 * @param function: Node of the function
 * @return true if everything was ok
 */
template <typename T>
bool
antok::generators::functionArgumentHandlerConst(mapStringTo<T>&   args,
                                                const YAML::Node& function)
{
	const YAML::Node& functionName = function["Name"];
	for (auto& arg : args) {
		const std::string& argName = arg.first;
		T&                 argVal  = arg.second;
		argVal = T();
		if (not hasNodeKey(function, argName)) {
			std::cerr << "Argument '" << argName << "' not found (required for function '" << functionName << "')." << std::endl;
			return false;
		}
		const YAML::Node& argNode = function[argName];
		try {
			argVal = argNode.as<T>();
		} catch (const YAML::TypedBadConversion<T>& e) {
			std::cerr << "Argument '" << argName << "' has wrong type (required for function '" << functionName << "')." << std::endl;
			return false;
		}
	}
	return true;
}

// specializations for common data types, add more as needed
template bool antok::generators::functionArgumentHandlerConst<double>     (mapStringTo<double>&      args, const YAML::Node& function);
template bool antok::generators::functionArgumentHandlerConst<int>        (mapStringTo<int>&         args, const YAML::Node& function);
template bool antok::generators::functionArgumentHandlerConst<std::string>(mapStringTo<std::string>& args, const YAML::Node& function);


/**
 * Sets the data pointers in the args vector to the address of the variable or to a constant if no variable name, but a number is given
 * @param args Vector of pairs where first: node/variable name, second: data pointer (will be set in this function)
 * @param function: Node of the function
 * @param index: Index of the function call (0 if this arguments have no index)
 * @return true if everything was ok
 */
template <typename T>
bool
antok::generators::functionArgumentHandlerPossibleConst(vecPairString<T*>& args,
                                                        const YAML::Node&  function,
                                                        const int          index)
{
	const YAML::Node& functionName = function["Name"];
	antok::Data&      data         = antok::ObjectManager::instance()->getData();

	// find all arguments which are given in the function node
	for (auto& arg : args) {
		if (hasNodeKey(function, arg.first)) {
			const YAML::Node& node = function[arg.first];
			std::string variable_name = antok::YAMLUtils::getString(node);
			try {
				// find unused name in antok::Data memory allocation
				const T val = node.as<T>();
				size_t nmb = 0;
				while (!data.insert<T>("__generatorsFunctionsConst" + std::to_string(nmb))) {
					if (nmb > 10000) { // TODO arbitrary limit
						std::cerr << "Could not allocate memory for antok::generators::functionArgumentHandlerPossibleConst in Function '" << functionName << "'.";
					}
					nmb++;
				}
				T* retval = data.getAddr<T>("__generatorsFunctionsConst" + std::to_string(nmb));
				*retval = val;
				arg.second = retval;
			} catch (const YAML::TypedBadConversion<T>& e) {  // test if variable is a variable name
				if (variable_name == "") {
					std::cerr << "Entry has to be either a variable name or a convertible type." << std::endl;
					return false;
				}
				variable_name = antok::generators::mergeNameIndex(variable_name, index);
				arg.second = data.getAddr<T>(variable_name);
				if (arg.second == nullptr) {
					std::cerr << "Cannot find variable << '" << variable_name << "' (required for function '" << functionName << "')." << std::endl;
					return false;
				}
			}
		} else {
			std::cerr << "Argument '" << arg.first << "' not found (required for function '" << functionName << "')." << std::endl;
			return false;
		}
	}

	return true;
}

// specializations for common data types, add more if needed
template bool antok::generators::functionArgumentHandlerPossibleConst<double>(vecPairString<double*>& args, const YAML::Node& function, int index);
template bool antok::generators::functionArgumentHandlerPossibleConst<int>   (vecPairString<int*>&    args, const YAML::Node& function, int index);


std::string
antok::generators::getFunctionArgumentHandlerErrorMsg(const std::vector<std::string>& quantityNames)
{
	std::stringstream msgStream;
	if (quantityNames.size() > 1) {
		msgStream << "Error when registering calculation for quantities '[";
		for (size_t i = 0; i < quantityNames.size() - 1; ++i) {
			msgStream << quantityNames[i] << ", ";
		}
		msgStream << quantityNames[quantityNames.size() - 1] << "]'." << std::endl;
	} else {
		msgStream << "Error when registering calculation for quantity '" << quantityNames[0] << "'." << std::endl;
	}
	return msgStream.str();
}


antok::Function*
antok::generators::generateAbs(const YAML::Node&               function,
                               const std::vector<std::string>& quantityNames,
                               const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	const std::string argName     = "Arg";
	const std::string argTypeName = antok::generators::getTypeOfArg(function, index, argName);
	vecPairString<std::string> args = {{argName, argTypeName}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	//TODO reduce boiler-plate code using helper function similar to __generateDiffHelper()
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (argTypeName == "std::vector<double>") {
		if (not data.insert<std::vector<double>>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}
	}	else {
		if (not data.insert<double>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}
	}

	using antok::functions::Abs;
	if        (argTypeName == "double") {
		return new Abs<double>             (*data.getAddr<double>(args[0].first),
		                                    *data.getAddr<double>(quantityName));
	} else if (argTypeName == "int") {
		return new Abs<int>                (*data.getAddr<int>(args[0].first),
		                                    *data.getAddr<double>(quantityName));
	} else if (argTypeName == "TVector3") {
		return new Abs<TVector3>           (*data.getAddr<TVector3>(args[0].first),
		                                    *data.getAddr<double>(quantityName));
	} else if (argTypeName == "std::vector<double>") {
		return new Abs<std::vector<double>>(*data.getAddr<std::vector<double>>(args[0].first),
		                                    *data.getAddr<std::vector<double>>(quantityName));
	} else {
		std::cerr << "'" << function["Name"] << "' is not (yet) implemented for input type '" << argTypeName << "'." << std::endl;
	}
	return nullptr;

}


antok::Function*
antok::generators::generateLog(const YAML::Node&               function,
                               const std::vector<std::string>& quantityNames,
                               const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	const std::string argName     = "Arg";
	const std::string argTypeName = antok::generators::getTypeOfArg(function, index, argName);
	vecPairString<std::string> args = {{argName, argTypeName}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get optional constant argument
	std::map<std::string, double> constArgs = {{"Base", 0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		// set default value
		constArgs["Base"] = std::numeric_limits<double>::quiet_NaN();
	}

	// Register output variables
	//TODO reduce boiler-plate code using helper function similar to __generateDiffHelper()
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (argTypeName == "std::vector<double>") {
		if (not data.insert<std::vector<double>>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}
	}
	else {
		if (not data.insert<double>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}
	}

	using antok::functions::Log;
	if        (argTypeName == "double") {
		return new Log<double>             (*data.getAddr<double>(args[0].first),
		                                    constArgs["Base"],
		                                    *data.getAddr<double>(quantityName));
	} else if (argTypeName == "int") {
		return new Log<int>                (*data.getAddr<int>(args[0].first),
		                                    constArgs["Base"],
		                                    *data.getAddr<double>(quantityName));
	} else if (argTypeName == "std::vector<double>") {
		return new Log<std::vector<double>>(*data.getAddr<std::vector<double>>(args[0].first),
		                                    constArgs["Base"],
		                                    *data.getAddr<std::vector<double>>(quantityName));
	} else {
		std::cerr << "'" << function["Name"] << "' is not (yet) implemented for input type '" << argTypeName << "'." << std::endl;
	}
	return nullptr;

}


antok::Function*
antok::generators::generateConvertIntToDouble(const YAML::Node&               function,
                                              const std::vector<std::string>& quantityNames,
                                              const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args = {{"Int", "int"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::functions::ConvertIntToDouble(*data.getAddr<int>(args[0].first), *data.getAddr<double>(quantityName));
}


namespace {

	template <typename T>
	antok::Function*
	__generateDiffHelper(antok::Data&                      data,
	                     const vecPairString<std::string>& args,
	                     const std::vector<std::string>&   quantityNames,
	                     const std::string&                quantityName)
	{
		// Register output variables
		if (not data.insert<T>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}

		return new antok::functions::Diff<T>(*data.getAddr<T>(args[0].first),
		                                     *data.getAddr<T>(args[1].first),
		                                     *data.getAddr<T>(quantityName));
	}

}  // anonymous namespace


antok::Function*
antok::generators::generateDiff(const YAML::Node&               function,
                                const std::vector<std::string>& quantityNames,
                                const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	const YAML::Node& functionName = function["Name"];
	const std::string arg1Name     = "Minuend";
	const std::string arg2Name     = "Subtrahend";
	const std::string arg1TypeName = antok::generators::getTypeOfArg(function, index, arg1Name);
	const std::string arg2TypeName = antok::generators::getTypeOfArg(function, index, arg2Name);
	if (arg1TypeName != arg2TypeName) {
		std::cerr << "Argument '" << arg1Name << "' (" << arg1TypeName << ") and '" << arg2Name << "' (" << arg2TypeName << ") "
		          << "have different types (must be the same for function '" << functionName << "')." << std::endl;
	}
	vecPairString<std::string> args
		= {{arg1Name, arg1TypeName},
		   {arg2Name, arg2TypeName}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if        (arg1TypeName == "double") {
		return __generateDiffHelper<double>             (data, args, quantityNames, quantityName);
	} else if (arg1TypeName == "int") {
		return __generateDiffHelper<int>                (data, args, quantityNames, quantityName);
	} else if (arg1TypeName == "std::vector<double>") {
		return __generateDiffHelper<std::vector<double>>(data, args, quantityNames, quantityName);
	} else {
		std::cerr << "'" << functionName << "' is not (yet) implemented for type '" << arg1TypeName << "' of "
		          << "variables '" << args[0].first << "' and '" << args[1].first << "'." << std::endl;
		return nullptr;
	}
}

namespace {

	template <typename T>
	antok::Function*
	__generateQuotientHelper(antok::Data&                      data,
	                         const vecPairString<std::string>& args,
	                         const std::vector<std::string>&   quantityNames,
	                         const std::string&                quantityName)
	{
		// Register output variables
		if (not data.insert<T>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}

		return new antok::functions::Quotient<T>(*data.getAddr<T>(args[0].first),
		                                         *data.getAddr<T>(args[1].first),
		                                         *data.getAddr<T>(quantityName));
	}

}  // anonymous namespace


antok::Function*
antok::generators::generateQuotient(const YAML::Node&               function,
                                    const std::vector<std::string>& quantityNames,
                                    const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	const YAML::Node& functionName = function["Name"];
	const std::string arg1Name     = "Dividend";
	const std::string arg2Name     = "Divisor";
	const std::string arg1TypeName = antok::generators::getTypeOfArg(function, index, arg1Name);
	const std::string arg2TypeName = antok::generators::getTypeOfArg(function, index, arg2Name);
	if (arg1TypeName != arg2TypeName) {
		std::cerr << "Argument '" << arg1Name << "' (" << arg1TypeName << ") and '" << arg2Name << "' (" << arg2TypeName << ") "
		          << "have different types (must be the same for function '" << functionName << "')." << std::endl;
	}
	vecPairString<std::string> args
		= {{arg1Name, arg1TypeName},
		   {arg2Name, arg2TypeName}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if        (arg1TypeName == "double") {
		return __generateQuotientHelper<double>             (data, args, quantityNames, quantityName);
	} else if (arg1TypeName == "int") {
		return __generateQuotientHelper<int>                (data, args, quantityNames, quantityName);
	} else if (arg1TypeName == "std::vector<double>") {
		return __generateQuotientHelper<std::vector<double>>(data, args, quantityNames, quantityName);
	} else {
		std::cerr << "'" << functionName << "' is not (yet) implemented for type '" << arg1TypeName << "' of "
		          << "variables '" << args[0].first << "' and '" << args[1].first << "'." << std::endl;
		return nullptr;
	}
}


namespace {

	template <typename T>
	antok::Function*
	__generateMulHelper(antok::Data&                      data,
	                    const vecPairString<std::string>& args,
	                    const std::vector<std::string>&   quantityNames,
	                    const std::string&                quantityName )
	{
		// Register output variables
		if (not data.insert<T>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}

		return new antok::functions::Mul<T>(*data.getAddr<T>(args[0].first),
		                                    *data.getAddr<T>(args[1].first),
		                                    *data.getAddr<T>(quantityName));
	}

}  // anonymous namespace


antok::Function*
antok::generators::generateMul(const YAML::Node&               function,
                               const std::vector<std::string>& quantityNames,
                               const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	const YAML::Node& functionName = function["Name"];
	const std::string arg1Name     = "Factor1";
	const std::string arg2Name     = "Factor2";
	const std::string arg1TypeName = antok::generators::getTypeOfArg(function, index, arg1Name);
	const std::string arg2TypeName = antok::generators::getTypeOfArg(function, index, arg2Name);
	if (arg1TypeName != arg2TypeName) {
		std::cerr << "Argument '" << arg1Name << "' (" << arg1TypeName << ") and '" << arg2Name << "' (" << arg2TypeName << ") "
		          << "have different types (must be the same for function '" << functionName << "')." << std::endl;
	}
	vecPairString<std::string> args
		= {{arg1Name, arg1TypeName},
		   {arg2Name, arg2TypeName}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if        (arg1TypeName == "double") {
		return __generateMulHelper<double>             (data, args, quantityNames, quantityName);
	} else if (arg1TypeName == "int") {
		return __generateMulHelper<int>                (data, args, quantityNames, quantityName);
	} else if (arg1TypeName == "std::vector<double>") {
		return __generateMulHelper<std::vector<double>>(data, args, quantityNames, quantityName);
	} else {
		std::cerr << "'" << functionName << "' is not (yet) implemented for type '" << arg1TypeName << "' of "
		          << "variables '" << args[0].first << "' and '" << args[1].first << "'." << std::endl;
		return nullptr;
	}
}


antok::Function*
antok::generators::generateEnergy(const YAML::Node&               function,
                                  const std::vector<std::string>& quantityNames,
                                  const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args = {{"Vector", "TLorentzVector"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::functions::Energy(*data.getAddr<TLorentzVector>(args[0].first), *data.getAddr<double>(quantityName));
}


antok::Function*
antok::generators::generateGetBeamLorentzVector(const YAML::Node&               function,
                                                const std::vector<std::string>& quantityNames,
                                                const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args
		= {{"dX",          "double"},
		   {"dY",          "double"},
		   {"XLorentzVec", "TLorentzVector"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const bool hasNodeKeyBeamMass   = antok::YAMLUtils::hasNodeKey(function, "BeamMass");
	const bool hasNodeKeyTargetMass = antok::YAMLUtils::hasNodeKey(function, "TargetMass");
	vecPairString<double*> possibleConstArgs;
	if (hasNodeKeyBeamMass)
		possibleConstArgs.push_back(std::pair<std::string, double*>("BeamMass", 0));
	if (hasNodeKeyTargetMass)
		possibleConstArgs.push_back(std::pair<std::string, double*>("TargetMass", 0));
	if (not antok::generators::functionArgumentHandlerPossibleConst(possibleConstArgs, function, 0)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	const double* massBeamAddr   = nullptr;
	const double* massTargetAddr = nullptr;
	// keep this order!!!
	if (hasNodeKeyTargetMass) {
		massTargetAddr = possibleConstArgs.back().second;
		possibleConstArgs.pop_back();
	}
	if (hasNodeKeyBeamMass) {
		massBeamAddr = possibleConstArgs.back().second;
		possibleConstArgs.pop_back();
	}

	// Register output variables
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (not data.insert<TLorentzVector>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::functions::GetBeamLorentzVec(*data.getAddr<double>(args[0].first),          // gradX
	                                               *data.getAddr<double>(args[1].first),          // gradY
	                                               *data.getAddr<TLorentzVector>(args[2].first),  // xLorentzVec
	                                               *data.getAddr<TLorentzVector>(quantityName),   // out
	                                               massBeamAddr,
	                                               massTargetAddr);
}


antok::Function*
antok::generators::generateGetGradXGradY(const YAML::Node&               function,
                                         const std::vector<std::string>& quantityNames,
                                         const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 2)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args = {{"Vector", "TLorentzVector"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data&         data = antok::ObjectManager::instance()->getData();
	std::vector<double*> quantityAddrs;
	for (size_t i = 0; i < quantityNames.size(); ++i) {
		if (not data.insert<double>(quantityNames[i])) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return nullptr;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return new antok::functions::GetGradXGradY(*data.getAddr<TLorentzVector>(args[0].first), *quantityAddrs[0], *quantityAddrs[1]);
}


antok::Function*
antok::generators::generateGetLorentzVectorAttributes(const YAML::Node&               function,
                                                      const std::vector<std::string>& quantityNames,
                                                      const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 5)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args = {{"Vector", "TLorentzVector"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data&         data = antok::ObjectManager::instance()->getData();
	std::vector<double*> quantityAddrs;
	for (size_t i = 0; i < quantityNames.size(); ++i) {
		if (not data.insert<double>(quantityNames[i])) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return nullptr;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return new antok::functions::GetLorentzVectorAttributes(*data.getAddr<TLorentzVector>(args[0].first),
	                                                        *quantityAddrs[0],
	                                                        *quantityAddrs[1],
	                                                        *quantityAddrs[2],
	                                                        *quantityAddrs[3],
	                                                        *quantityAddrs[4]);
}


antok::Function*
antok::generators::generateGetLorentzVec(const YAML::Node&               function,
                                         const std::vector<std::string>& quantityNames,
                                         const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}
	// Get input variables
	const YAML::Node&  functionName = function["Name"];
	const std::string& quantityName = quantityNames[0];
	using antok::functions::GetLorentzVec;
	GetLorentzVec::lorentzVecDefType defType;
	if        (function["X"] and function["M"]) {
		defType = GetLorentzVec::XYZM;
	} else if (function["Px"] and function["E"]) {
		defType = GetLorentzVec::PxPyPzE;
	} else if (function["Vec3"] and function["M"]) {
		defType = GetLorentzVec::Vec3M;
	} else if (function["Vec3"] and function["E"]) {
		defType = GetLorentzVec::Vec3E;
	} else {
		std::cerr << "Function '" << functionName << "' needs either input variables "
		          << "'[X, Y, Z, M]', '[Px, Py, Pz, E]', '[Vec3, M]', or '[Vec3, E]' "
		          << "to calculate variable '" << quantityName << "'." << std::endl;
		return nullptr;
	}
	antok::Data& data = antok::ObjectManager::instance()->getData();
	vecPairString<std::string> args;
	double* mAddr = nullptr;
	if (defType == GetLorentzVec::XYZM or defType == GetLorentzVec::Vec3M) {
		size_t nmb = 0;
		// search for unused memory space
		while (!data.insert<double>("__getLorentzVecMassVec" + std::to_string(nmb))) {
			if (nmb > 10000) { // TODO arbitrary limit
				std::cerr << "Could not allocate memory for __getLorentzVecMassVec in Function '" << functionName << "'.";
			}
			nmb++;
		}
		mAddr = data.getAddr<double>("__getLorentzVecMassVec" + std::to_string(nmb));
		std::map<std::string, double> constArgs = {{"M", 0}};
		if (not antok::generators::functionArgumentHandlerConst<double>(constArgs, function)) {
			std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
			return nullptr;
		}
		(*mAddr) = constArgs["M"];
	}
	switch (defType) {
		case GetLorentzVec::XYZM: {
			args = {{"X", "double"},
			        {"Y", "double"},
			        {"Z", "double"}};
			break;
		}
		case GetLorentzVec::PxPyPzE:
			args = {{"Px", "double"},
			        {"Py", "double"},
			        {"Pz", "double"},
			        {"E",  "double"}};
			break;
		case GetLorentzVec::Vec3M: {
			args = {{"Vec3", "TVector3"}};
			break;
		}
		case GetLorentzVec::Vec3E:
			args = {{"Vec3", "TVector3"},
			        {"E",    "double"}};
			break;
	}
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	// Register output variables
	if (not data.insert<TLorentzVector>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}
	switch (defType) {
		case GetLorentzVec::XYZM:
			return new GetLorentzVec(*data.getAddr<double>(args[0].first),         // x
			                         *data.getAddr<double>(args[1].first),         // y
			                         *data.getAddr<double>(args[2].first),         // z
			                         *mAddr,                                       // m
			                         *data.getAddr<TLorentzVector>(quantityName),  // out
			                         defType);
		case GetLorentzVec::PxPyPzE:
			return new GetLorentzVec(*data.getAddr<double>(args[0].first),         // x
			                         *data.getAddr<double>(args[1].first),         // y
			                         *data.getAddr<double>(args[2].first),         // z
			                         *data.getAddr<double>(args[3].first),         // E
			                         *data.getAddr<TLorentzVector>(quantityName),  // out
			                         defType);
		case GetLorentzVec::Vec3M:
			return new GetLorentzVec(*data.getAddr<TVector3>(args[0].first),       // (x, y, z)^T
			                         *mAddr,                                       // m
			                         *data.getAddr<TLorentzVector>(quantityName),  // out
			                         defType);
		case GetLorentzVec::Vec3E:
			return new GetLorentzVec(*data.getAddr<TVector3>(args[0].first),       // (x, y, z)^T
			                         *data.getAddr<double>(args[1].first),         // E
			                         *data.getAddr<TLorentzVector>(quantityName),  // out
			                         defType);
		default: {
			std::cerr << "Unclear what to do in function '" << functionName << "'." << std::endl;
			return nullptr;
		}
	}
	return nullptr;
}


antok::Function*
antok::generators::generateGetTs(const YAML::Node&               function,
                                 const std::vector<std::string>& quantityNames,
                                 const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 3)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args
		= {{"BeamLorentzVec", "TLorentzVector"},
		   {"XLorentzVec",    "TLorentzVector"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	vecPairString<double*> possibleConstArgs = {{"TargetMass", nullptr}};
	const double* targetMassAddr = nullptr;
	if (hasNodeKey(function, "TargetMass")) {
		if (not antok::generators::functionArgumentHandlerPossibleConst(possibleConstArgs, function, 0)) {
			std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
			return nullptr;
		}
		targetMassAddr = possibleConstArgs[0].second;
	}

	// Register output variables
	antok::Data&         data = antok::ObjectManager::instance()->getData();
	std::vector<double*> quantityAddrs;
	for (size_t i = 0; i < quantityNames.size(); ++i) {
		if (not data.insert<double>(quantityNames[i])) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return nullptr;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return new antok::functions::GetTs(*data.getAddr<TLorentzVector>(args[1].first),             // XLorentzVec
	                                   *data.getAddr<TLorentzVector>(args[0].first),             // BeamLorentzVec
	                                   *quantityAddrs[0], *quantityAddrs[1], *quantityAddrs[2],  // out
	                                   targetMassAddr);
}


antok::Function*
antok::generators::generateGetVector3(const YAML::Node&               function,
                                      const std::vector<std::string>& quantityNames,
                                      const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	const YAML::Node&  functionName = function["Name"];
	const std::string& quantityName = quantityNames[0];
	vecPairString<std::string> args;
	enum vecDefType {fromCoords = 0, fromVectors = 1, fromTLorentzVector = 2};
	vecDefType defType;
	if (hasNodeKey(function, "X")) {
		defType = fromCoords;
		args = {{"X", "double"},
		        {"Y", "double"},
		        {"Z", "double"}};
	} else if (hasNodeKey(function,"VectorX")) {
		defType = fromVectors;
		args = {{"VectorX", "std::vector<double>"},
		        {"VectorY", "std::vector<double>"},
		        {"VectorZ", "std::vector<double>"}};
	} else if (hasNodeKey(function,"LVector")) {
		defType = fromTLorentzVector;
		args = {{"LVector", "TLorentzVector"}};
	} else {
		std::cerr << "Function '" << functionName << "' needs either input variables "
		          << "'[X, Y, Z]', '[VectorX, VectorY, VectorZ]', or 'LVector' "
		          << "to calculate variable '" << quantityName << "'." << std::endl;
		return nullptr;
	}
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables and create functor
	antok::Data& data = antok::ObjectManager::instance()->getData();
	switch (defType) {
		case fromCoords: {
			if (not data.insert<TVector3>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
				return nullptr;
			}
			return new antok::functions::GetTVector3(*data.getAddr<double>(args[0].first),    // x
			                                         *data.getAddr<double>(args[1].first),    // y
			                                         *data.getAddr<double>(args[2].first),    // z
			                                         *data.getAddr<TVector3>(quantityName));  // out
		}
		case fromVectors: {
			if (not data.insert<std::vector<TVector3>>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
				return nullptr;
			}
			return new antok::functions::GetVectorTVector3(*data.getAddr<std::vector<double>>(args[0].first),    // x
			                                               *data.getAddr<std::vector<double>>(args[1].first),    // y
			                                               *data.getAddr<std::vector<double>>(args[2].first),    // z
			                                               *data.getAddr<std::vector<TVector3>>(quantityName));  // out
		}
		case fromTLorentzVector: {
			if (not data.insert<TVector3>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
				return nullptr;
			}
			return new antok::functions::GetTVector3FromTLorenzVector(*data.getAddr<TLorentzVector>(args[0].first),
			                                                          *data.getAddr<TVector3>(quantityName));
		}
		default: {
			std::cerr << "Unclear what to do in function '" << functionName << "'." << std::endl;
			return nullptr;
		}
	}

	return nullptr;
}


namespace {

	// registers output variables and creates functor
	template <typename T>
	antok::Function*
	__getVectorEntryFunction(const int&                        entry,
	                         const vecPairString<std::string>& args,
	                         const std::vector<std::string>&   quantityNames)
	{
		const std::string& quantityName = quantityNames[0];
		antok::Data&       data         = antok::ObjectManager::instance()->getData();
		if (not data.insert<T>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}

		return new antok::functions::GetVectorEntry<T>(*data.getAddr<std::vector<T>>(args[0].first),  // vector
		                                               entry,
		                                               *data.getAddr<T>(quantityName));               // result
	}

}


antok::Function*
antok::generators::generateGetVectorEntry(const YAML::Node&               function,
                                          const std::vector<std::string>& quantityNames,
                                          const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	const YAML::Node& functionName = function["Name"];
	vecPairString<std::string> args;
	enum vectorElementType {Int = 0, Double = 1, Vector3 = 2, LorentzVector = 3};
	vectorElementType elemType;
	if        (hasNodeKey(function, "VectorInt")) {
		elemType = Int;
		args     = {{"VectorInt", "std::vector<int>"}};
	} else if (hasNodeKey(function, "VectorDouble")) {
		elemType = Double;
		args     = {{"VectorDouble", "std::vector<double>"}};
	} else if (hasNodeKey(function, "VectorTVector3")) {
		elemType = Vector3;
		args     = {{"VectorTVector3", "std::vector<TVector3>"}};
	} else if (hasNodeKey(function, "VectorTLorentzVector")) {
		elemType = LorentzVector;
		args     = {{"VectorTLorentzVector", "std::vector<TLorentzVector>"}};
	} else {
		std::cerr << "Function '" << functionName << "' needs either input variables "
		          << "'VectorInt', 'VectorDouble', 'VectorTVector3', or 'VectorTLorentzVector' "
		          << "to calculate variable '" << quantityNames[0] << "'." << std::endl;
		return nullptr;
	}
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, int> constArgs = {{"Entry", 0}};
	if (not antok::generators::functionArgumentHandlerConst<int>(constArgs, function)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// return functor
	switch (elemType) {
		case Int:
		  return __getVectorEntryFunction<int>           (constArgs["Entry"], args, quantityNames);
		case Double:
		  return __getVectorEntryFunction<double>        (constArgs["Entry"], args, quantityNames);
		case Vector3:
		  return __getVectorEntryFunction<TVector3>      (constArgs["Entry"], args, quantityNames);
		case LorentzVector:
		  return __getVectorEntryFunction<TLorentzVector>(constArgs["Entry"], args, quantityNames);
		default: {
			std::cerr << "Unclear what to do in function '" << functionName << "'." << std::endl;
			return nullptr;
		}
	}
	return nullptr;
}


antok::Function*
antok::generators::generateMass(const YAML::Node&               function,
                                const std::vector<std::string>& quantityNames,
                                const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args = {{"Vector", "TLorentzVector"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::functions::Mass(*data.getAddr<TLorentzVector>(args[0].first),
	                                  *data.getAddr<double>(quantityName));
}


antok::Function*
antok::generators::generateMass2(const YAML::Node&               function,
                                 const std::vector<std::string>& quantityNames,
                                 const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args = {{"Vector", "TLorentzVector"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::functions::Mass2(*data.getAddr<TLorentzVector>(args[0].first),
	                                   *data.getAddr<double>(quantityName));
}


antok::Function*
antok::generators::generateRadToDegree(const YAML::Node&               function,
                                       const std::vector<std::string>& quantityNames,
                                       const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	antok::Data& data = antok::ObjectManager::instance() ->getData();
	vecPairString<std::string> args = {{"Angle", "double"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& quantityName = quantityNames[0];
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::functions::RadToDegree(*data.getAddr<double>(args[0].first), *data.getAddr<double>(quantityName));
}


namespace {

	template <typename T>
	antok::Function*
	__getSumFunction(const vecPairString<std::string>& summandNames,
	                 const vecPairString<std::string>& subtrahendsNames,
	                 const std::string&                quantityName)
	{
		antok::Data& data = antok::ObjectManager::instance()->getData();

		// Do type checking and get all the addresses
		// for summands
		std::vector<T*> inputAddrsSummands;
		for (size_t summandNames_i = 0; summandNames_i < summandNames.size(); ++summandNames_i) {
			const std::string& variableName = summandNames[summandNames_i].first;
			inputAddrsSummands.push_back(data.getAddr<T>(variableName));
		}
		// for subtrahends
		std::vector<T*> inputAddrsSubtrahends;
		for (size_t subtrahendsNames_i = 0; subtrahendsNames_i < subtrahendsNames.size(); ++subtrahendsNames_i) {
			const std::string& variableName = subtrahendsNames[subtrahendsNames_i].first;
			inputAddrsSubtrahends.push_back(data.getAddr<T>(variableName));
		}

		// Register output variables
		if (not data.insert<T>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
			return nullptr;
		}

		return new antok::functions::Sum<T>(inputAddrsSummands, inputAddrsSubtrahends, *data.getAddr<T>(quantityName));
	}


	vecPairString<std::string>*
	__getSummandNames(const YAML::Node&  function,
	                  const std::string& quantityName,
	                  const int          index,
	                  const std::string& type)
	{
		if (not hasNodeKey(function, type.c_str())) {
			return nullptr;
		}

		antok::Data& data     = antok::ObjectManager::instance()->getData();
		std::string  typeName = "notInitialized";
		vecPairString<std::string>* summandNames = new vecPairString<std::string>();
		const bool hasNodeKeyIndices = hasNodeKey(function[type.c_str()], "Indices");
		const bool hasNodeKeyName    = hasNodeKey(function[type.c_str()], "Name");
		if (hasNodeKeyIndices or hasNodeKeyName) {

			// Summing over one variable with indices
			if (not (hasNodeKeyIndices and hasNodeKeyName)) {
				std::cerr << "Either 'Indices' or 'Name' found in sum function, but not both (Variable: '" << quantityName << "')." << std::endl;
				return nullptr;
			}
			if (index > 0) {
				std::cerr << "Cannot have sum over indices for every particle (Variable: '" << quantityName << "')." << std::endl;
				return nullptr;
			}
			std::vector<int> inner_indices;
			try {
				inner_indices = function[type.c_str()]["Indices"].as<std::vector<int>>();
			} catch (const YAML::TypedBadConversion<std::vector<int>>& e) {
				std::cerr << "Could not convert YAML sequence to std::vector<int> when parsing 'Indices' for sum function "
				          << "(for variable '" << quantityName << "')." << std::endl;
				return nullptr;
			} catch (const YAML::TypedBadConversion<int>& e) {
				std::cerr << "Could not convert entries in YAML sequence to int when parsing 'Indices' for sum function "
				          << "(for variable '" << quantityName << "')." << std::endl;
				return nullptr;
			}
			typeName = antok::YAMLUtils::getString(function[type.c_str()]["Name"]);
			if (typeName == "") {
				std::cerr << "Could not convert 'Name' of summand to std::string when registering calculation of '" << quantityName << "'." << std::endl;
			}
			const std::string summandBaseName = typeName;
			{
				std::stringstream strStr;
				strStr << typeName << inner_indices[0];
				typeName = data.getType(strStr.str());
			}
			for (size_t inner_indices_i = 0; inner_indices_i < inner_indices.size(); ++inner_indices_i) {
				int inner_index = inner_indices[inner_indices_i];
				std::stringstream strStr;
				strStr << summandBaseName << inner_index;
				summandNames->push_back(std::pair<std::string, std::string>(strStr.str(), typeName));
			}

		} else {

			// Summing over list of variable names
			typeName = antok::YAMLUtils::getString(*function[type.c_str()].begin());
			if (typeName == "") {
				std::cerr << "Could not convert one of the 'Summands' to std::string when registering calculation of '" << quantityName << "'." << std::endl;
				return nullptr;
			}
			if (index > 0) {
				std::stringstream strStr;
				strStr << typeName << index;
				typeName = strStr.str();
			}
			typeName = data.getType(typeName);
			for (YAML::const_iterator summand_it = function[type.c_str()].begin(); summand_it != function[type.c_str()].end(); ++summand_it) {
				const std::string& variableName = antok::YAMLUtils::getString(*summand_it);
				if (variableName == "") {
					std::cerr << "Could not convert one of the 'Summands' to std::string when registering calculation of '" << quantityName << "'." << std::endl;
					return nullptr;
				}
				summandNames->push_back(std::pair<std::string, std::string>(variableName, typeName));
			}

		}

		return summandNames;
	}

}  // anonymous namespace


antok::Function*
antok::generators::generateSum(const YAML::Node&               function,
                               const std::vector<std::string>& quantityNames,
                               const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	const YAML::Node&  functionName = function["Name"];
	const std::string& quantityName = quantityNames[0];
	vecPairString<std::string>* summandNamesPtr    = __getSummandNames(function, quantityName, index, "Summands");
	vecPairString<std::string>* subtrahendNamesPtr = __getSummandNames(function, quantityName, index, "Subtrahends");
	if ((summandNamesPtr == nullptr) and (subtrahendNamesPtr == nullptr)) {
		std::cerr << "Could not generate summands for function '" << functionName << "' when trying to register calculation of '" << quantityName << "'." << std::endl;
		return nullptr;
	}
	vecPairString<std::string> summandNames;
	if (summandNamesPtr != nullptr) {
		summandNames = (*summandNamesPtr);
	}
	if (not antok::generators::functionArgumentHandler(summandNames, function, index, true)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	vecPairString<std::string> subtrahendNames;
	if (subtrahendNamesPtr != nullptr) {
		subtrahendNames = (*subtrahendNamesPtr);
	}
	if (not antok::generators::functionArgumentHandler(subtrahendNames, function, index, true)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	std::string typeName;
	if (summandNamesPtr != nullptr) {
		typeName = summandNames[0].second;
	} else {
		typeName = subtrahendNames[0].second;
	}

	antok::Function* antokFunction = nullptr;
	if        (typeName == "double") {
		antokFunction = __getSumFunction<double>        (summandNames, subtrahendNames, quantityName);
	} else if (typeName == "int") {
		antokFunction = __getSumFunction<int>           (summandNames, subtrahendNames, quantityName);
	} else if (typeName == "Long64_t") {
		antokFunction = __getSumFunction<Long64_t>      (summandNames, subtrahendNames, quantityName);
	} else if (typeName == "TLorentzVector") {
		antokFunction = __getSumFunction<TLorentzVector>(summandNames, subtrahendNames, quantityName);
	} else if (typeName == "TVector3") {
		antokFunction = __getSumFunction<TVector3>      (summandNames, subtrahendNames, quantityName);
	} else {
		std::cerr << "'" << functionName << "' is not (yet) implemented for type '" << typeName << "' "
		          << "(calculation of '" << quantityName << "')." << std::endl;
		return nullptr;
	}
	delete summandNamesPtr;
	delete subtrahendNamesPtr;
	return antokFunction;
}


antok::Function*
antok::generators::generateSum2(const YAML::Node&               function,
                                const std::vector<std::string>& quantityNames,
                                const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Define arguments
	const YAML::Node&  functionName = function["Name"];
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	vecPairString<std::string>* summandNamesPtr = __getSummandNames(function, quantityName, index, "Summands");
	if (summandNamesPtr == nullptr) {
		std::cerr << "Could not generate summands for function '" << functionName << "' "
		          << "when trying to register calculation of '" << quantityName << "'." << std::endl;
		return nullptr;
	}
	vecPairString<std::string>& summandNames = (*summandNamesPtr);
	if (not antok::generators::functionArgumentHandler(summandNames, function, index, true)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Now do type checking and get all the addresses
	std::vector<double*> doubleInputAddrs;
	for (size_t summandNames_i = 0; summandNames_i < summandNames.size(); ++summandNames_i) {
			const std::string& variableName = summandNames[summandNames_i].first;
			doubleInputAddrs.push_back(data.getAddr<double>(variableName));
		}

	// Register output variables
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}
	delete summandNamesPtr;

	return new antok::functions::Sum2(doubleInputAddrs, *data.getAddr<double>(quantityName));
}
