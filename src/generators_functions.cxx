#include<generators_functions.h>

#include<assert.h>

#include<TVector3.h>
#include<TLorentzVector.h>

#include<data.h>
#include<functions.hpp>
#include<initializer.h>
#include<object_manager.h>
#include<yaml_utils.hpp>
#include<constants.h>

std::string antok::generators::mergeNameIndex( std::string const& name, int const index ){
		if(index > 0) {
			std::stringstream strStr;
			strStr<<name<<index;
			return strStr.str();
		}
		return name;
}

bool antok::generators::functionArgumentHandler(std::vector<std::pair<std::string, std::string> >& args,
                                                const YAML::Node& function,
                                                int index,
                                                bool argStringsAlreadyValues)
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
		argName = antok::generators::mergeNameIndex( argName, index );
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


/**
 * Sets the data pointers in the args vector to the address of the variable or to an constant if no variable name, but a number is given
 * @param args Vector of pairs where first: node/variable name, second: data pointer (will be set in this function)
 * @param function: Node of the function
 * @param index: Index of the function call (0 if this arguments have no index)
 * @return true if everything was ok
 */
template< typename T>
bool antok::generators::functionrgumentHandlerPossibleConst( std::vector< std::pair< std::string, T* > >& args,
		                                  const YAML::Node& function,
		                                  int index) {

	using antok::YAMLUtils::hasNodeKey;
	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector< std::pair< std::string, std::string > > given_args;
	std::map< int, int > map_args_given_args;

	// find all arguments which are given in the function node
	for( size_t i = 0; i < args.size(); ++i ){
		auto& arg = args[i];
		if( hasNodeKey(function, arg.first) ){
			const YAML::Node& node = function[arg.first];
			try {
				const T val = node.as<T>();
				arg.second = new T(val);
			} catch (const YAML::TypedBadConversion<T>& e) { // test if variable is a variable name
				std::string variable_name = antok::YAMLUtils::getString( node );
				if(variable_name == "") {
					std::cerr<<"Entry has to be either a variable name or a convertible type."<<std::endl;
					return false;
				}
				variable_name = antok::generators::mergeNameIndex(variable_name, index);
				arg.second = data.getAddr<T>(variable_name);
				if( arg.second == nullptr ){
					std::cerr<<"Can not find variable << \"" << variable_name << "\" (required for function \""<<function["Name"]<<"\")."<<std::endl;
					return false;
				}
			}

		} else {
			std::cerr<<"Argument \""<<arg.first<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return false;
		}
	}


	return true;
}
// initialize for a bunch of data types, add more if needed
template bool antok::generators::functionrgumentHandlerPossibleConst<double>( std::vector< std::pair< std::string, double* > >& args, const YAML::Node& function, int index);
template bool antok::generators::functionrgumentHandlerPossibleConst<int>( std::vector< std::pair< std::string, int* > >& args, const YAML::Node& function, int index);

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

namespace {

	template<typename T>
	antok::Function* __getSumFunction(std::vector<std::pair<std::string, std::string> >& summandNames,
	                                  std::vector<std::pair<std::string, std::string> >& subtrahendsNames,
	                                  std::string quantityName) {

		std::vector<T*> inputAddrsSummands;
		std::vector<T*> inputAddrsSubtrahends;
		antok::Data& data = antok::ObjectManager::instance()->getData();

		// Now do type checking and get all the addresses
		// for summands
		for(unsigned int summandNames_i = 0; summandNames_i < summandNames.size(); ++summandNames_i) {

			std::string variableName = summandNames[summandNames_i].first;
			inputAddrsSummands.push_back(data.getAddr<T>(variableName));

		}
		// for subtrahends
		for(unsigned int subtrahendsNames_i = 0; subtrahendsNames_i < subtrahendsNames.size(); ++subtrahendsNames_i) {

			std::string variableName = subtrahendsNames[subtrahendsNames_i].first;
			inputAddrsSubtrahends.push_back(data.getAddr<T>(variableName));

		}

		// And produce the function
		if(not data.insert<T>(quantityName)) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityName);
			return 0;
		}
		return (new antok::functions::Sum<T>(inputAddrsSummands, inputAddrsSubtrahends, data.getAddr<T>(quantityName)));

	};

	std::vector<std::pair<std::string, std::string> >* __getSummandNames(const YAML::Node& function, std::string& quantityName, int index, std::string type) {

		using antok::YAMLUtils::hasNodeKey;

		std::string typeName = "notInitialized";
		std::vector<std::pair<std::string, std::string> >* summandNames = new std::vector<std::pair<std::string, std::string> >();

		// Summing over one variable with indices
		if(not hasNodeKey(function, type.c_str())) {
			return 0;
		}

		antok::Data& data = antok::ObjectManager::instance()->getData();

		if(hasNodeKey(function[type.c_str()], "Indices") or hasNodeKey(function[type.c_str()], "Name")) {
			if(not (hasNodeKey(function[type.c_str()], "Indices") and hasNodeKey(function[type.c_str()], "Name"))) {
				std::cerr<<"Either \"Indices\" or \"Name\" found in sum function, but not both (Variable: \""<<quantityName<<"\")."<<std::endl;
				return 0;
			}
			if(index > 0) {
				std::cerr<<"Cannot have sum over indices for every particle (Variable: \""<<quantityName<<"\")."<<std::endl;
				return 0;
			}
			std::vector<int> inner_indices;
			try {
				inner_indices = function[type.c_str()]["Indices"].as<std::vector<int> >();
			} catch (const YAML::TypedBadConversion<std::vector<int> >& e) {
				std::cerr<<"Could not convert YAML sequence to std::vector<int> when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
				return 0;
			} catch (const YAML::TypedBadConversion<int>& e) {
				std::cerr<<"Could not convert entries in YAML sequence to int when parsing \"sum\"' \"Indices\" (for variable \""<<quantityName<<"\")."<<std::endl;
				return 0;
			}
			typeName = antok::YAMLUtils::getString(function[type.c_str()]["Name"]);
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
			typeName = antok::YAMLUtils::getString(*function[type.c_str()].begin());
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
			for(YAML::const_iterator summand_it = function[type.c_str()].begin(); summand_it != function[type.c_str()].end(); ++summand_it) {
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
	antok::Data& data = antok::ObjectManager::instance()->getData();

	using antok::YAMLUtils::hasNodeKey;
	std::string typeNameArg1;
	if( hasNodeKey(function, "Arg") )	typeNameArg1 = data.getType( antok::generators::mergeNameIndex( antok::YAMLUtils::getString( function["Arg"] ), index ) );
	else {
		std::cerr<<"Argument \"Arg\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Arg", typeNameArg1));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}



	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	if      ( typeNameArg1 == "double" )
		return (new antok::functions::Abs<double>(data.getAddr<double>(args[0].first), data.getAddr<double>(quantityName)));
	else if ( typeNameArg1 == "int" )
		return (new antok::functions::Abs<int>(data.getAddr<int>(args[0].first), data.getAddr<double>(quantityName)));
	else if ( typeNameArg1 == "TVector3" )
		return (new antok::functions::Abs<TVector3>(data.getAddr<TVector3>(args[0].first), data.getAddr<double>(quantityName)));
	else
		std::cerr<<"abs is not (yet) implemented for input type '" << typeNameArg1 << "'."<<std::endl;
	return 0;

};

antok::Function* antok::generators::generateLog(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];
	antok::Data& data = antok::ObjectManager::instance()->getData();

	using antok::YAMLUtils::hasNodeKey;
	std::string typeNameArg1;
	if( hasNodeKey(function, "Arg") )	typeNameArg1 = data.getType( antok::generators::mergeNameIndex( antok::YAMLUtils::getString( function["Arg"] ), index ) );
	else {
		std::cerr<<"Argument \"Arg\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Arg", typeNameArg1));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}



	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	if      ( typeNameArg1 == "double" )
		return (new antok::functions::Log<double>(data.getAddr<double>(args[0].first), data.getAddr<double>(quantityName)));
	else if ( typeNameArg1 == "int" )
		return (new antok::functions::Log<int>(data.getAddr<int>(args[0].first), data.getAddr<double>(quantityName)));
	else
		std::cerr<<"abs is not (yet) implemented for input type '" << typeNameArg1 << "'."<<std::endl;
	return 0;

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

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
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

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
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

template< typename T >
antok::Function* __generateQuotientHelper( antok::Data& data, const std::vector<std::pair<std::string, std::string> >& args,
		const std::vector<std::string>& quantityNames, const std::string& quantityName  ){

	if(not data.insert<T>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}
	return (new antok::functions::Quotient<T>( data.getAddr<T>(args[0].first),
													data.getAddr<T>(args[1].first),
													data.getAddr<T>(quantityName)));

}

antok::Function* antok::generators::generateQuotient(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::string quantityName = quantityNames[0];


	using antok::YAMLUtils::hasNodeKey;

	// Get type for arguments
	std::string  typeNameArg1, typeNameArg2;
	if( hasNodeKey(function, "Dividend") )	typeNameArg1 = data.getType( antok::generators::mergeNameIndex( antok::YAMLUtils::getString( function["Dividend"] ), index ) );
	else {
		std::cerr<<"Argument \"Dividend\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
		return 0;
	}
	if( hasNodeKey(function, "Divisor") )	typeNameArg2 = data.getType( antok::generators::mergeNameIndex( antok::YAMLUtils::getString( function["Divisor"] ), index ) );
	else {
		std::cerr<<"Argument \"Divisor\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
		return 0;
	}

	if ( typeNameArg1 != typeNameArg2 ){
		std::cerr<<"Argument \"Dividend\" (" << typeNameArg1 << ") and \"Divisor\" (" << typeNameArg2 << ") have different types (required for function \""<<function["Name"]<<"\")."<<std::endl;

	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Dividend", typeNameArg1));
	args.push_back(std::pair<std::string, std::string>("Divisor", typeNameArg2));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}


	if      (typeNameArg1 == "double") return __generateQuotientHelper<double>(data, args, quantityNames, quantityName);
	else if (typeNameArg1 == "int")    return __generateQuotientHelper<int>(data, args, quantityNames, quantityName);
	else {
		std::cerr << "Divide not implemented for type \"" << typeNameArg1 << "\" for variables \"" << args[0].first << "\" and \"" << args[1].first << "\"." << std::endl;
		return 0;
	}

	return 0;

};

template< typename T >
antok::Function* __generateMulHelper( antok::Data& data, const std::vector<std::pair<std::string, std::string> >& args,
		const std::vector<std::string>& quantityNames, const std::string& quantityName  ){

	if(not data.insert<T>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}
	return (new antok::functions::Mul<T>( data.getAddr<T>(args[0].first),
													data.getAddr<T>(args[1].first),
													data.getAddr<T>(quantityName)));

}

antok::Function* antok::generators::generateMul(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::string quantityName = quantityNames[0];


	using antok::YAMLUtils::hasNodeKey;

	// Get type for arguments
	std::string  typeNameArg1, typeNameArg2;
	if( hasNodeKey(function, "Factor1") )	typeNameArg1 = data.getType( antok::generators::mergeNameIndex( antok::YAMLUtils::getString( function["Factor1"] ), index ) );
	else {
		std::cerr<<"Argument \"Factor1\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
		return 0;
	}
	if( hasNodeKey(function, "Factor2") )	typeNameArg2 = data.getType( antok::generators::mergeNameIndex( antok::YAMLUtils::getString( function["Factor2"] ), index ) );
	else {
		std::cerr<<"Argument \"Factor2\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
		return 0;
	}


	if ( typeNameArg1 != typeNameArg2 ){
		std::cerr<<"Argument \"Factor1\" (" << typeNameArg1 << ") and \"Factor2\" (" << typeNameArg2 << ") have different types (required for function \""<<function["Name"]<<"\")."<<std::endl;

	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Factor1", typeNameArg1));
	args.push_back(std::pair<std::string, std::string>("Factor2", typeNameArg2));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}


	if      (typeNameArg1 == "double") return __generateMulHelper<double>(data, args, quantityNames, quantityName);
	else if (typeNameArg1 == "int")    return __generateMulHelper<int>(data, args, quantityNames, quantityName);
	else {
		std::cerr << "Mul not implemented for type \"" << typeNameArg1 << "\" for variables \"" << args[0].first << "\" and \"" << args[1].first << "\"." << std::endl;
		return 0;
	}

	return 0;

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

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
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

antok::Function* antok::generators::generateComponents(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 4) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName0 = quantityNames[0];
	std::string quantityName1 = quantityNames[1];
	std::string quantityName2 = quantityNames[2];
	std::string quantityName3 = quantityNames[3];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector", "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	TLorentzVector* arg1Addr = data.getAddr<TLorentzVector>(args[0].first);

	if(not data.insert<double>(quantityName0)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}
	if(not data.insert<double>(quantityName1)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}
	if(not data.insert<double>(quantityName2)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}
	if(not data.insert<double>(quantityName3)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::Components(arg1Addr, data.getAddr<double>(quantityName0),
					data.getAddr<double>(quantityName1), data.getAddr<double>(quantityName2), data.getAddr<double>(quantityName3)));
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
	std::vector<std::pair<std::string, double*>> possible_const_args;

	args.push_back(std::pair<std::string, std::string>("dX", "double"));
	args.push_back(std::pair<std::string, std::string>("dY", "double"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	const double* beam_mass = nullptr;
	const double* target_mass = nullptr;
	if( antok::YAMLUtils::hasNodeKey( function, "BeamMass") )
		possible_const_args.push_back(std::pair<std::string, double*>("BeamMass", 0));
	if( antok::YAMLUtils::hasNodeKey( function, "TargetMass") )
		possible_const_args.push_back(std::pair<std::string, double*>("TargetMass", 0));

	if(not antok::generators::functionrgumentHandlerPossibleConst(possible_const_args, function, 0)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	// keep this order!!!
	if( antok::YAMLUtils::hasNodeKey( function, "TargetMass") ){
		target_mass = possible_const_args.back().second;
		possible_const_args.pop_back();
	}
	if( antok::YAMLUtils::hasNodeKey( function, "BeamMass") ){
		beam_mass = possible_const_args.back().second;
		possible_const_args.pop_back();
	}

	double* dXaddr = data.getAddr<double>(args[0].first);
	double* dYaddr = data.getAddr<double>(args[1].first);
	TLorentzVector* xLorentzVecAddr = data.getAddr<TLorentzVector>(args[2].first);

	if(not data.insert<TLorentzVector>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::functions::GetBeamLorentzVec(dXaddr, dYaddr, xLorentzVecAddr, data.getAddr<TLorentzVector>(quantityName), beam_mass, target_mass));

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

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
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

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
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

	int pType;
	if(function["X"] and function["M"]) {
		pType = 0;
	} else if (function["Px"] and function["E"]) {
		pType = 1;
	} else if (function["Vec3"] and function["M"]) {
		pType = 2;
	} else if (function["Vec3"] and function["E"]) {
		pType = 3;
	} else {
		std::cerr<<"Function \"getLorentzVec\" needs either variables \"[X, Y, Z, M]\" or \"[Px, Py, Pz, E]\" (variable \""<<quantityName<<"\")."<<std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;

	double* xAddr;
	double* yAddr;
	double* zAddr;
	double* mAddr;
	TVector3* vec3Addr;

	switch(pType)
	{
		case 0:
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
			break;
		case 1:
			args.push_back(std::pair<std::string, std::string>("Px", "double"));
			args.push_back(std::pair<std::string, std::string>("Py", "double"));
			args.push_back(std::pair<std::string, std::string>("Pz", "double"));
			args.push_back(std::pair<std::string, std::string>("E", "double"));
			break;
		case 2:
			args.push_back(std::pair<std::string, std::string>("Vec3", "TVector3"));
			try {
				function["M"].as<double>();
			} catch(const YAML::TypedBadConversion<double>& e) {
				std::cerr<<"Argument \"M\" in function \"mass\" should be of type double (variable \""<<quantityName<<"\")."<<std::endl;
				return 0;
			}
			mAddr = new double();
			(*mAddr) = function["M"].as<double>();
			break;
		case 3:
			args.push_back(std::pair<std::string, std::string>("Vec", "TVector3"));
			args.push_back(std::pair<std::string, std::string>("E", "double"));
			break;
	}

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	switch(pType)
	{
		case 0:
			xAddr = data.getAddr<double>(args[0].first);
			yAddr = data.getAddr<double>(args[1].first);
			zAddr = data.getAddr<double>(args[2].first);
			break;
		case 1:
			xAddr = data.getAddr<double>(args[0].first);
			yAddr = data.getAddr<double>(args[1].first);
			zAddr = data.getAddr<double>(args[2].first);
			mAddr = data.getAddr<double>(args[3].first);
			break;
		case 2:
			vec3Addr = data.getAddr<TVector3>(args[0].first);
			mAddr = data.getAddr<double>(args[1].first);
			break;
		case 3:
			vec3Addr = data.getAddr<TVector3>(args[0].first);
			mAddr = data.getAddr<double>(args[1].first);
			break;
	}

	if(not data.insert<TLorentzVector>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	switch(pType)
	{
		case 0:
			return (new antok::functions::GetLorentzVec(xAddr, yAddr, zAddr, mAddr, data.getAddr<TLorentzVector>(quantityName), pType));
		case 1:
			return (new antok::functions::GetLorentzVec(xAddr, yAddr, zAddr, mAddr, data.getAddr<TLorentzVector>(quantityName), pType));
		case 2:
			return (new antok::functions::GetLorentzVec(vec3Addr, mAddr, data.getAddr<TLorentzVector>(quantityName), pType));
		case 3:
			return (new antok::functions::GetLorentzVec(vec3Addr, mAddr, data.getAddr<TLorentzVector>(quantityName), pType));
	}
	return 0;
};

antok::Function* antok::generators::generateGetTs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	using antok::YAMLUtils::hasNodeKey;
	if(quantityNames.size() != 3) {
		std::cerr<<"Need 3 names for function \"getTs\""<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;

	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));


	std::vector<std::pair<std::string, double*> > possible_const;
	possible_const.push_back(std::pair<std::string, double*>("TargetMass", nullptr));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}


	TLorentzVector* beamLVAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* xLVAddr = data.getAddr<TLorentzVector>(args[1].first);


	const double* targetMassAddr = NULL;
	if( hasNodeKey(function, "TargetMass") ){
		if(not antok::generators::functionrgumentHandlerPossibleConst(possible_const, function, 0)) {
			std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
			return 0;
		}
		targetMassAddr = possible_const[0].second;
	}

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return (new antok::functions::GetTs(xLVAddr, beamLVAddr, quantityAddrs[0], quantityAddrs[1], quantityAddrs[2], targetMassAddr));

};

antok::Function* antok::generators::generateGetVector3(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
	using antok::YAMLUtils::hasNodeKey;

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	const std::string& quantityName = quantityNames[0];

	antok::Data& data = antok::ObjectManager::instance()->getData();

	bool fromTLorentzVector;
	std::vector<std::pair<std::string, std::string> > args;
	if(hasNodeKey(function, "X")){
		fromTLorentzVector = false;
		args.push_back(std::pair<std::string, std::string>("X", "double"));
		args.push_back(std::pair<std::string, std::string>("Y", "double"));
		args.push_back(std::pair<std::string, std::string>("Z", "double"));
	} else {
		fromTLorentzVector = true;
		args.push_back(std::pair<std::string, std::string>("LVector", "TLorentzVector"));
	}

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}


	if(not data.insert<TVector3>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	TVector3* outVec = data.getAddr<TVector3>(quantityName);
	if(not fromTLorentzVector){
		double* xAddr = data.getAddr<double>(args[0].first);
		double* yAddr = data.getAddr<double>(args[1].first);
		double* zAddr = data.getAddr<double>(args[2].first);
		return (new antok::functions::GetTVector3(xAddr, yAddr, zAddr, outVec));
	} else {
		return new antok::functions::GetTVector3FromTLorenzVector(data.getAddr<TLorentzVector>(args[0].first), outVec);
	}

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

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
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

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
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

	std::vector<std::pair<std::string, std::string> >* summandNamesPtr = __getSummandNames(function, quantityName, index, "Summands");
	std::vector<std::pair<std::string, std::string> >* subtrahendNamesPtr = __getSummandNames(function, quantityName, index, "Subtrahends");
	if((summandNamesPtr == 0) and (subtrahendNamesPtr == 0)) {
		std::cerr<<"Could not generate summands for function \"Sum\" when trying to register calculation of \""<<quantityName<<"\"."<<std::endl;
		return 0;
	}
	std::vector<std::pair<std::string, std::string> > summandNames;
	std::vector<std::pair<std::string, std::string> > subtrahendNames;
	if(summandNamesPtr != 0)
		summandNames = (*summandNamesPtr);
	if(subtrahendNamesPtr != 0)
		subtrahendNames = (*subtrahendNamesPtr);

	if(not antok::generators::functionArgumentHandler(summandNames, function, index, true)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}
	if(not antok::generators::functionArgumentHandler(subtrahendNames, function, index, true)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}
	std::string typeName;
	if(summandNamesPtr != 0)
		typeName = summandNames[0].second;
	else
		typeName = subtrahendNames[0].second;

	antok::Function* antokFunction = 0;
	if(typeName == "double") {
		antokFunction = __getSumFunction<double>(summandNames, subtrahendNames, quantityName);
	} else if (typeName == "int") {
		antokFunction = __getSumFunction<int>(summandNames, subtrahendNames, quantityName);
	} else if (typeName == "Long64_t") {
		antokFunction = __getSumFunction<Long64_t>(summandNames, subtrahendNames, quantityName);
	} else if (typeName == "TLorentzVector") {
		antokFunction = __getSumFunction<TLorentzVector>(summandNames, subtrahendNames, quantityName);
	} else if (typeName == "TVector3") {
		antokFunction = __getSumFunction<TVector3>(summandNames, subtrahendNames, quantityName);
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

	std::vector<std::pair<std::string, std::string> >* summandNamesPtr = __getSummandNames(function, quantityName, index, "Summands");
	if(summandNamesPtr == 0) {
		std::cerr<<"Could not generate summands for function \"Sum\" when trying to register calculation of \""<<quantityName<<"\"."<<std::endl;
		return 0;
	}
	std::vector<std::pair<std::string, std::string> >& summandNames = (*summandNamesPtr);

	if(not antok::generators::functionArgumentHandler(summandNames, function, index, true)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
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

