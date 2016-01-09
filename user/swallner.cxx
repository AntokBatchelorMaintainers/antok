
#include<swallner.h>

#include<iostream>
#include<sstream>
#include<string>
#include<functions.hpp>
#include<generators_functions.h>
#include<data.h>
#include<yaml_utils.hpp>

#include<swallner_functions.hpp>



antok::Function* antok::user::stefan::getCalcLTProjections(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
	if(quantityNames.size() != 2) {
		std::cerr<<"Need 2 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return nullptr;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args_per_index;
	std::vector<std::pair<std::string, std::string> > args;
	args_per_index.push_back(std::pair<std::string, std::string>("Vector", "TVector3"));
	args.push_back(std::pair<std::string, std::string>("Direction", "TVector3"));

	if(not antok::generators::functionArgumentHandler(args, function, 0)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}
	if(not antok::generators::functionArgumentHandler(args_per_index, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}


	TVector3* vector = data.getAddr<TVector3>(args_per_index[0].first);
	TVector3* direction = data.getAddr<TVector3>(args[0].first);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return new antok::user::stefan::functions::CalcLTProjection(vector, direction, quantityAddrs[0], quantityAddrs[1]);
}

antok::Function* antok::user::stefan::getCalcArmenterosAlpha(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
	if(quantityNames.size() != 1) {
		std::cerr<<"Need 1 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return nullptr;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("LongitudinalMom1", "double"));
	args.push_back(std::pair<std::string, std::string>("LongitudinalMom2", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, 0)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}


	double* longitudinal_mom_1 = data.getAddr<double>(args[0].first);
	double* longitudinal_mom_2 = data.getAddr<double>(args[1].first);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return new antok::user::stefan::functions::CalcArmenterosAlpha(longitudinal_mom_1, longitudinal_mom_2, quantityAddrs[0]);
}


antok::Function* antok::user::stefan::getUserFunction(const YAML::Node& function,
                                                          std::vector<std::string>& quantityNames,
                                                          int index)
{
	std::string functionName = antok::YAMLUtils::getString(function["Name"]);
	antok::Function* antokFunctionPtr = 0;
	if(functionName == "calcLTProjections") {
		antokFunctionPtr = stefan::getCalcLTProjections(function, quantityNames, index);
	}
	if(functionName == "calcArmenterosAlpha") {
		antokFunctionPtr = stefan::getCalcArmenterosAlpha(function, quantityNames, index);
	}
	return antokFunctionPtr;
}
