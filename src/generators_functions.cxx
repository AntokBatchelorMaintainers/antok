#include <cassert>

#include "TVector3.h"
#include "TLorentzVector.h"

#include "constants.h"
#include "data.h"
#include "functions.hpp"
#include "generators_functions.h"
#include "initializer.h"
#include "object_manager.h"
#include "yaml_utils.hpp"


using antok::YAMLUtils::hasNodeKey;


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
antok::generators::functionArgumentHandler(std::vector<std::pair<std::string, std::string> >& args,
                                           const YAML::Node&                                  function,
                                           const int                                          index,
                                           const bool                                         argStringsAlreadyValues)
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
 * Sets the data pointers in the args vector to the address of the variable or to a constant if no variable name, but a number is given
 * @param args Vector of pairs where first: node/variable name, second: data pointer (will be set in this function)
 * @param function: Node of the function
 * @param index: Index of the function call (0 if this arguments have no index)
 * @return true if everything was ok
 */
template <typename T>
bool
antok::generators::functionArgumentHandlerPossibleConst(std::vector< std::pair< std::string, T* > >& args,
                                                        const YAML::Node&                            function,
                                                        const int                                    index)
{
	const YAML::Node& functionName = function["Name"];
	antok::Data&      data         = antok::ObjectManager::instance()->getData();

	// find all arguments which are given in the function node
	for (size_t i = 0; i < args.size(); ++i) {
		auto& arg = args[i];
		if (hasNodeKey(function, arg.first)) {
			const YAML::Node& node = function[arg.first];
			try {
				const T val = node.as<T>();
				arg.second = new T(val);
			} catch (const YAML::TypedBadConversion<T>& e) {  // test if variable is a variable name
				std::string variable_name = antok::YAMLUtils::getString(node);
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


// initialize for a bunch of data types, add more if needed
template bool antok::generators::functionArgumentHandlerPossibleConst<double>(std::vector< std::pair< std::string, double* > >& args, const YAML::Node& function, int index);
template bool antok::generators::functionArgumentHandlerPossibleConst<int>   (std::vector< std::pair< std::string, int* > >&    args, const YAML::Node& function, int index);


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


//TODO there is a lot of boiler-plate code that is repeated in the functions below
//     -> move this code into separate functions
antok::Function*
antok::generators::generateAbs(const YAML::Node&               function,
                               const std::vector<std::string>& quantityNames,
                               const int                       index)
{
	const YAML::Node& functionName = function["Name"];
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << functionName << "'." << std::endl;
		return nullptr;
	}

	// Get type of argument
	antok::Data& data = antok::ObjectManager::instance()->getData();
	std::string  typeNameArg1;
	if (hasNodeKey(function, "Arg")) {
		typeNameArg1 = data.getType(antok::generators::mergeNameIndex(antok::YAMLUtils::getString(function["Arg"]), index));
	} else {
		std::cerr << "Argument 'Arg' not found (required for function '" << functionName << "')." << std::endl;
		return nullptr;
	}
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Arg", typeNameArg1));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	if        (typeNameArg1 == "double") {
		return new antok::functions::Abs<double>  (*data.getAddr<double>(args[0].first),   *data.getAddr<double>(quantityName));
	} else if (typeNameArg1 == "int") {
		return new antok::functions::Abs<int>     (*data.getAddr<int>(args[0].first),      *data.getAddr<double>(quantityName));
	} else if (typeNameArg1 == "TVector3") {
		return new antok::functions::Abs<TVector3>(*data.getAddr<TVector3>(args[0].first), *data.getAddr<double>(quantityName));
	} else {
		std::cerr << "'" << functionName << "' is not (yet) implemented for input type '" << typeNameArg1 << "'." << std::endl;
	}
	return nullptr;

}


antok::Function*
antok::generators::generateLog(const YAML::Node&               function,
                               const std::vector<std::string>& quantityNames,
                               const int                       index)
{
	const YAML::Node& functionName = function["Name"];
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << functionName << "'." << std::endl;
		return nullptr;
	}

	// Get type of argument
	antok::Data& data = antok::ObjectManager::instance()->getData();
	std::string  typeNameArg1;
	if (hasNodeKey(function, "Arg")) {
		typeNameArg1 = data.getType(antok::generators::mergeNameIndex(antok::YAMLUtils::getString(function["Arg"]), index));
	} else {
		std::cerr << "Argument 'Arg' not found (required for function '" << functionName << "')." << std::endl;
		return nullptr;
	}
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Arg", typeNameArg1));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	if        (typeNameArg1 == "double") {
		return new antok::functions::Log<double>(*data.getAddr<double>(args[0].first), *data.getAddr<double>(quantityName));
	} else if (typeNameArg1 == "int") {
		return new antok::functions::Log<int>   (*data.getAddr<int>(args[0].first),    *data.getAddr<double>(quantityName));
	} else {
		std::cerr << "'" << functionName << "' is not (yet) implemented for input type '" << typeNameArg1 << "'." << std::endl;
	}
	return nullptr;

}


antok::Function*
antok::generators::generateConvertIntToDouble(const YAML::Node&               function,
                                              const std::vector<std::string>& quantityNames,
                                              const int                       index)
{
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << function["Name"] << "'." << std::endl;
		return nullptr;
	}

	// Define type of argument
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Int", "int"));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::functions::ConvertIntToDouble(*data.getAddr<int>(args[0].first), *data.getAddr<double>(quantityName));
}


antok::Function*
antok::generators::generateDiff(const YAML::Node&               function,
                                const std::vector<std::string>& quantityNames,
                                const int                       index)
{
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << function["Name"] << "'." << std::endl;
		return nullptr;
	}

	// Define types of arguments
	//TODO support other types like done in antok::generators::generateQuotient()
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Minuend",    "double"));
	args.push_back(std::pair<std::string, std::string>("Subtrahend", "double"));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::functions::Diff(*data.getAddr<double>(args[0].first),
	                                  *data.getAddr<double>(args[1].first),
	                                  *data.getAddr<double>(quantityName));
}


namespace {

	template <typename T>
	antok::Function*
	__generateQuotientHelper(antok::Data&                                             data,
	                         const std::vector<std::pair<std::string, std::string> >& args,
	                         const std::vector<std::string>&                          quantityNames,
	                         const std::string&                                       quantityName)
	{
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
	const YAML::Node& functionName = function["Name"];
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << functionName << "'." << std::endl;
		return nullptr;
	}

	// Get types of arguments
	antok::Data& data = antok::ObjectManager::instance()->getData();
	std::string  typeNameArg1, typeNameArg2;
	if (hasNodeKey(function, "Dividend")) {
		typeNameArg1 = data.getType(antok::generators::mergeNameIndex(antok::YAMLUtils::getString(function["Dividend"]), index));
	} else {
		std::cerr << "Argument 'Dividend' not found (required for function '" << functionName << "')." << std::endl;
		return nullptr;
	}
	if (hasNodeKey(function, "Divisor")) {
		typeNameArg2 = data.getType(antok::generators::mergeNameIndex(antok::YAMLUtils::getString(function["Divisor"]), index));
	} else {
		std::cerr << "Argument 'Divisor' not found (required for function '" << functionName << "')." << std::endl;
		return nullptr;
	}
	if (typeNameArg1 != typeNameArg2) {
		std::cerr << "Argument 'Dividend' (" << typeNameArg1 << ") and 'Divisor' (" << typeNameArg2 << ") "
		          << "have different types (required for function '" << functionName << "')." << std::endl;
	}
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Dividend", typeNameArg1));
	args.push_back(std::pair<std::string, std::string>("Divisor",  typeNameArg2));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	if        (typeNameArg1 == "double") {
		return __generateQuotientHelper<double>(data, args, quantityNames, quantityName);
	} else if (typeNameArg1 == "int") {
		return __generateQuotientHelper<int>   (data, args, quantityNames, quantityName);
	} else {
		std::cerr << "'" << functionName << "' not implemented for type '" << typeNameArg1 << "' for "
		          << "variables '" << args[0].first << "' and '" << args[1].first << "'." << std::endl;
		return nullptr;
	}
}


namespace {

	template <typename T>
	antok::Function*
	__generateMulHelper(antok::Data&                                             data,
	                    const std::vector<std::pair<std::string, std::string> >& args,
	                    const std::vector<std::string>&                          quantityNames,
	                    const std::string& quantityName )
	{
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
	const YAML::Node& functionName = function["Name"];
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << functionName << "'." << std::endl;
		return nullptr;
	}

	// Get type of arguments
	antok::Data& data = antok::ObjectManager::instance()->getData();
	std::string  typeNameArg1, typeNameArg2;
	if (hasNodeKey(function, "Factor1")) {
		typeNameArg1 = data.getType(antok::generators::mergeNameIndex(antok::YAMLUtils::getString(function["Factor1"]), index));
	} else {
		std::cerr << "Argument 'Factor1' not found (required for function '" << functionName << "')." << std::endl;
		return nullptr;
	}
	if (hasNodeKey(function, "Factor2")) {
		typeNameArg2 = data.getType(antok::generators::mergeNameIndex(antok::YAMLUtils::getString(function["Factor2"]), index));
	} else {
		std::cerr << "Argument 'Factor2' not found (required for function '" << functionName << "')." << std::endl;
		return nullptr;
	}
	if (typeNameArg1 != typeNameArg2) {
		std::cerr << "Argument 'Factor1' (" << typeNameArg1 << ") and 'Factor2' (" << typeNameArg2 << ") have different types "
		          << "(required for function '" << functionName << "')." << std::endl;
	}
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Factor1", typeNameArg1));
	args.push_back(std::pair<std::string, std::string>("Factor2", typeNameArg2));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	if        (typeNameArg1 == "double") {
		return __generateMulHelper<double>(data, args, quantityNames, quantityName);
	} else if (typeNameArg1 == "int") {
		return __generateMulHelper<int>   (data, args, quantityNames, quantityName);
	} else {
		std::cerr << "'" << functionName << "' not (yet) implemented for type '" << typeNameArg1 << "' for "
		          << "variables '" << args[0].first << "' and '" << args[1].first << "'." << std::endl;
		return nullptr;
	}
}


antok::Function*
antok::generators::generateEnergy(const YAML::Node&               function,
                                  const std::vector<std::string>& quantityNames,
                                  const int                       index)
{
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << function["Name"] << "'." << std::endl;
		return nullptr;
	}

	// Define type of argument
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

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
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << function["Name"] << "'." << std::endl;
		return nullptr;
	}

	// Define types of arguments
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("dX",          "double"));
	args.push_back(std::pair<std::string, std::string>("dY",          "double"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const double* massBeamAddr   = nullptr;
	const double* massTargetAddr = nullptr;
	std::vector<std::pair<std::string, double*>> possible_const_args;
	const bool hasNodeKeyBeamMass   = antok::YAMLUtils::hasNodeKey(function, "BeamMass");
	const bool hasNodeKeyTargetMass = antok::YAMLUtils::hasNodeKey(function, "TargetMass");
	if (hasNodeKeyBeamMass)
		possible_const_args.push_back(std::pair<std::string, double*>("BeamMass", 0));
	if (hasNodeKeyTargetMass)
		possible_const_args.push_back(std::pair<std::string, double*>("TargetMass", 0));
	if (not antok::generators::functionArgumentHandlerPossibleConst(possible_const_args, function, 0)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	// keep this order!!!
	if (hasNodeKeyTargetMass) {
		massTargetAddr = possible_const_args.back().second;
		possible_const_args.pop_back();
	}
	if (hasNodeKeyBeamMass) {
		massBeamAddr = possible_const_args.back().second;
		possible_const_args.pop_back();
	}

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
	if (quantityNames.size() != 2) {
		std::cerr << "Need 2 names for function '" << function["Name"] << "'." << std::endl;
		return nullptr;
	}

	// Define type of argument
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

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
	if (quantityNames.size() != 5) {
		std::cerr << "Need 5 names for function '" << function["Name"] << "'." << std::endl;
		return nullptr;
	}

	// Define type of argument
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

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
	const YAML::Node& functionName = function["Name"];
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << functionName << "'." << std::endl;
		return nullptr;
	}

	// Define types of arguments
	const std::string& quantityName = quantityNames[0];
	int pType;
	if        (function["X"] and function["M"]) {
		pType = 0;
	} else if (function["Px"] and function["E"]) {
		pType = 1;
	} else if (function["Vec3"] and function["M"]) {
		pType = 2;
	} else if (function["Vec3"] and function["E"]) {
		pType = 3;
	} else {
		std::cerr << "Function '" << functionName << "' needs either variables '[X, Y, Z, M]' or '[Px, Py, Pz, E]' "
		          << "(variable '" << quantityName << "')." << std::endl;
		return nullptr;
	}
	std::vector<std::pair<std::string, std::string> > args;
	double* mAddr = nullptr;
	switch (pType) {
		case 0:
			args.push_back(std::pair<std::string, std::string>("X", "double"));
			args.push_back(std::pair<std::string, std::string>("Y", "double"));
			args.push_back(std::pair<std::string, std::string>("Z", "double"));
			try {
				function["M"].as<double>();
			} catch (const YAML::TypedBadConversion<double>& e) {
				std::cerr << "Argument 'M' in function '" << functionName << "' should be of type double (variable '" << quantityName << "')." << std::endl;
				return nullptr;
			}
			mAddr    = new double();
			(*mAddr) = function["M"].as<double>();
			break;
		case 1:
			args.push_back(std::pair<std::string, std::string>("Px", "double"));
			args.push_back(std::pair<std::string, std::string>("Py", "double"));
			args.push_back(std::pair<std::string, std::string>("Pz", "double"));
			args.push_back(std::pair<std::string, std::string>("E",  "double"));
			break;
		case 2:
			args.push_back(std::pair<std::string, std::string>("Vec3", "TVector3"));
			try {
				function["M"].as<double>();
			} catch (const YAML::TypedBadConversion<double>& e) {
				std::cerr << "Argument 'M' in function '" << functionName << "' should be of type double (variable '" << quantityName << "')." << std::endl;
				return nullptr;
			}
			mAddr    = new double();
			(*mAddr) = function["M"].as<double>();
			break;
		case 3:
			args.push_back(std::pair<std::string, std::string>("Vec", "TVector3"));
			args.push_back(std::pair<std::string, std::string>("E",   "double"));
			break;
	}
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();
	if (not data.insert<TLorentzVector>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	switch (pType) {
		case 0:
			return new antok::functions::GetLorentzVec(*data.getAddr<double>(args[0].first),         // x
			                                           *data.getAddr<double>(args[1].first),         // y
			                                           *data.getAddr<double>(args[2].first),         // z
			                                           *mAddr,                                       // m
			                                           *data.getAddr<TLorentzVector>(quantityName),  // out
			                                           pType);
		case 1:
			return new antok::functions::GetLorentzVec(*data.getAddr<double>(args[0].first),         // x
			                                           *data.getAddr<double>(args[1].first),         // y
			                                           *data.getAddr<double>(args[2].first),         // z
			                                           *data.getAddr<double>(args[3].first),         // m
			                                           *data.getAddr<TLorentzVector>(quantityName),  // out
			                                           pType);
		case 2:
			return new antok::functions::GetLorentzVec(*data.getAddr<TVector3>(args[0].first),       // (x, y, z)^T
			                                           *mAddr,                                       // m
			                                           *data.getAddr<TLorentzVector>(quantityName),  // out
			                                           pType);
		case 3:
			return new antok::functions::GetLorentzVec(*data.getAddr<TVector3>(args[0].first),       // (x, y, z)^T
			                                           *data.getAddr<double>(args[1].first),         // m
			                                           *data.getAddr<TLorentzVector>(quantityName),  // out
			                                           pType);
	}
	return nullptr;
}


antok::Function*
antok::generators::generateGetTs(const YAML::Node&               function,
                                 const std::vector<std::string>& quantityNames,
                                 const int                       index)
{
	if (quantityNames.size() != 3) {
		std::cerr << "Need 3 names for function '" << function["Name"] << "'" << std::endl;
		return nullptr;
	}

	// Define types of arguments
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec",    "TLorentzVector"));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	std::vector<std::pair<std::string, double*> > possible_const;
	possible_const.push_back(std::pair<std::string, double*>("TargetMass", nullptr));
	const double* targetMassAddr = nullptr;
	if (hasNodeKey(function, "TargetMass")) {
		if (not antok::generators::functionArgumentHandlerPossibleConst(possible_const, function, 0)) {
			std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
			return nullptr;
		}
		targetMassAddr = possible_const[0].second;
	}

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
	const YAML::Node& functionName = function["Name"];
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << functionName << "'." << std::endl;
		return nullptr;
	}

	// Define types of arguments
	std::vector<std::pair<std::string, std::string> > args;
	//TODO maybe it would be better to use an enum for the cases instead of 2 booleans
	bool fromTLorentzVector = false;
	bool fromVectors        = false;
	if (hasNodeKey(function, "X")) {
		args.push_back(std::pair<std::string, std::string>("X", "double"));
		args.push_back(std::pair<std::string, std::string>("Y", "double"));
		args.push_back(std::pair<std::string, std::string>("Z", "double"));
	} else if (hasNodeKey(function,"VectorX")) {
		fromTLorentzVector = false;
		fromVectors        = true;
		args.push_back(std::pair<std::string, std::string>("VectorX", "std::vector<double>"));
		args.push_back(std::pair<std::string, std::string>("VectorY", "std::vector<double>"));
		args.push_back(std::pair<std::string, std::string>("VectorZ", "std::vector<double>"));
	} else if (hasNodeKey(function,"LVector")) {
		fromTLorentzVector = true;
		fromVectors        = false;
		args.push_back(std::pair<std::string, std::string>("LVector", "TLorentzVector"));
	} else {
		std::cerr << "Unknown properties for function '" << functionName << "'." << std::endl;
		return nullptr;
	}
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (fromVectors) {
		if (fromTLorentzVector) {
			std::cerr << "Unclear what to do in function '" << functionName << "'." << std::endl;
			return nullptr;
		} else {
			if (not data.insert<std::vector<TVector3>>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
				return nullptr;
			}
			return new antok::functions::GetVectorTVector3(*data.getAddr<std::vector<double>>(args[0].first),    // x
			                                               *data.getAddr<std::vector<double>>(args[1].first),    // y
			                                               *data.getAddr<std::vector<double>>(args[2].first),    // z
			                                               *data.getAddr<std::vector<TVector3>>(quantityName));  // out
		}
	} else {
		if (not data.insert<TVector3>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
			return nullptr;
		}
		if (fromTLorentzVector) {
			return new antok::functions::GetTVector3FromTLorenzVector(*data.getAddr<TLorentzVector>(args[0].first), *data.getAddr<TVector3>(quantityName));
		} else {
			return new antok::functions::GetTVector3(*data.getAddr<double>(args[0].first),    // x
			                                         *data.getAddr<double>(args[1].first),    // y
			                                         *data.getAddr<double>(args[2].first),    // z
			                                         *data.getAddr<TVector3>(quantityName));  // out
		}
	}
}


antok::Function*
antok::generators::generateGetVectorEntry(const YAML::Node&               function,
                                          const std::vector<std::string>& quantityNames,
                                          const int                       index)
{
	const YAML::Node& functionName = function["Name"];
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << functionName << "'." << std::endl;
		return nullptr;
	}

	// Define types of arguments
	std::vector<std::pair<std::string, std::string> > args;
	int type;
	if        (hasNodeKey(function, "VectorInt")) {
		type = 0;
		args.push_back(std::pair<std::string, std::string>("VectorInt", "std::vector<int>"));
	} else if (hasNodeKey(function, "VectorDouble")) {
		type = 1;
		args.push_back(std::pair<std::string, std::string>("VectorDouble", "std::vector<double>"));
	} else if (hasNodeKey(function, "VectorTVector3")) {
		type = 2;
		args.push_back(std::pair<std::string, std::string>("VectorTVector3", "std::vector<TVector3>"));
	} else if (hasNodeKey(function, "VectorTLorentzVector")) {
		type = 3;
		args.push_back(std::pair<std::string, std::string>("VectorTLorentzVector", "std::vector<TLorentzVector>"));
	} else {
		std::cerr << "Unknown properties for function '" << functionName << "'." << std::endl;
		return nullptr;
	}
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	try {
		function["Entry"].as<int>();
	} catch (const YAML::InvalidNode &exception) {
		std::cerr << "Argument 'Entry' in function '" << functionName << "' not found for calculation of variables '" << quantityNames[0] << "'." << std::endl;
		return nullptr;
	} catch (const YAML::TypedBadConversion<int> &exception) {
		std::cerr << "Argument 'Entry' in function '" << functionName << "' should be of type int for calculation of variables '" << quantityNames[0] << "'." << std::endl;
		return nullptr;
	}
	const int entry = function["Entry"].as<int>();

	//TODO introduce template helper function to avoid repetition
	antok::Data& data = antok::ObjectManager::instance()->getData();
	switch (type) {
		case 0: {
			if (not data.insert<int>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
				return nullptr;
			}
			return new antok::functions::GetVectorEntry<int>(*data.getAddr<std::vector<int>>(args[0].first),  // vector
			                                                 entry,
			                                                 *data.getAddr<int>(quantityName));               // result
		}
		case 1: {
			if (not data.insert<double>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
				return nullptr;
			}
			return new antok::functions::GetVectorEntry<double>(*data.getAddr<std::vector<double>>(args[0].first),  // vector
			                                                    entry,
			                                                    *data.getAddr<double>(quantityName));               // result
		}
		case 2: {
			if (not data.insert<TVector3>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
				return nullptr;
			}
			return new antok::functions::GetVectorEntry<TVector3>(*data.getAddr<std::vector<TVector3>>(args[0].first),  // vector
			                                                      entry,
			                                                      *data.getAddr<TVector3>(quantityName));               // result
		}
		case 3: {
			if (not data.insert<TLorentzVector>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
				return nullptr;
			}
			return new antok::functions::GetVectorEntry<TLorentzVector>(*data.getAddr<std::vector<TLorentzVector>>(args[0].first),  // vector
			                                                            entry,
			                                                            *data.getAddr<TLorentzVector>(quantityName));               // result
		}
	}
	return nullptr;
}


antok::Function*
antok::generators::generateMass(const YAML::Node&               function,
                                const std::vector<std::string>& quantityNames,
                                const int                       index)
{
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << function["Name"] << "'." << std::endl;
		return nullptr;
	}

	// Define type of argument
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::functions::Mass(*data.getAddr<TLorentzVector>(args[0].first), *data.getAddr<double>(quantityName));
}


antok::Function*
antok::generators::generateRadToDegree(const YAML::Node&               function,
                                       const std::vector<std::string>& quantityNames,
                                       const int                       index)
{
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << function["Name"] << "'." << std::endl;
		return nullptr;
	}

	// Define types of arguments
	antok::Data& data = antok::ObjectManager::instance() ->getData();
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Angle", "double"));
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

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
	__getSumFunction(const std::vector<std::pair<std::string, std::string> >& summandNames,
	                 const std::vector<std::pair<std::string, std::string> >& subtrahendsNames,
	                 const std::string&                                       quantityName)
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

		// Produce the function
		if (not data.insert<T>(quantityName)) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
			return nullptr;
		}
		return new antok::functions::Sum<T>(inputAddrsSummands, inputAddrsSubtrahends, *data.getAddr<T>(quantityName));
	}


	std::vector<std::pair<std::string, std::string> >*
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
		std::vector<std::pair<std::string, std::string> >* summandNames
		  = new std::vector<std::pair<std::string, std::string> >();
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
				inner_indices = function[type.c_str()]["Indices"].as<std::vector<int> >();
			} catch (const YAML::TypedBadConversion<std::vector<int> >& e) {
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
	const YAML::Node& functionName = function["Name"];
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << functionName << "'." << std::endl;
		return nullptr;
	}

	// Define types of arguments
	const std::string& quantityName = quantityNames[0];
	std::vector<std::pair<std::string, std::string> >* summandNamesPtr    = __getSummandNames(function, quantityName, index, "Summands");
	std::vector<std::pair<std::string, std::string> >* subtrahendNamesPtr = __getSummandNames(function, quantityName, index, "Subtrahends");
	if ((summandNamesPtr == nullptr) and (subtrahendNamesPtr == nullptr)) {
		std::cerr << "Could not generate summands for function '" << functionName << "' when trying to register calculation of '" << quantityName << "'." << std::endl;
		return nullptr;
	}
	std::vector<std::pair<std::string, std::string> > summandNames;
	if (summandNamesPtr != nullptr) {
		summandNames = (*summandNamesPtr);
	}
	if (not antok::generators::functionArgumentHandler(summandNames, function, index, true)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	std::vector<std::pair<std::string, std::string> > subtrahendNames;
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
		std::cerr << "Type '" << typeName << "' is not (yet) supported by '" << functionName << "' "
		          << "(registering calculation of '" << quantityName << "')." << std::endl;
		return nullptr;
	}
	delete summandNamesPtr;
	//TODO why isn't subtrahendNamesPtr deleted too?
	return antokFunction;
}


antok::Function*
antok::generators::generateSum2(const YAML::Node&               function,
                                const std::vector<std::string>& quantityNames,
                                const int                       index)
{
	const YAML::Node& functionName = function["Name"];
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function '" << functionName << "'." << std::endl;
		return nullptr;
	}

	// Define arguments
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	std::vector<std::pair<std::string, std::string> >* summandNamesPtr = __getSummandNames(function, quantityName, index, "Summands");
	if (summandNamesPtr == nullptr) {
		std::cerr << "Could not generate summands for function '" << functionName << "' "
		          << "when trying to register calculation of '" << quantityName << "'." << std::endl;
		return nullptr;
	}
	std::vector<std::pair<std::string, std::string> >& summandNames = (*summandNamesPtr);
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

	// And produce the function
	if (not data.insert<double>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}
	delete summandNamesPtr;
	return new antok::functions::Sum2(doubleInputAddrs, *data.getAddr<double>(quantityName));
}
