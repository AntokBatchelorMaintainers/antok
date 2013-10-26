
#include<hubers.h>

#include<constants.h>
#include<data.h>
#include<functions.hpp>
#include<generators_functions.h>
#include<hubers_functions.hpp>
#include<yaml_utils.hpp>
#include<iostream>

antok::Function* antok::user::hubers::getUserFunction(const YAML::Node& function,
		std::vector<std::string>& quantityNames,
		int index)
{
	std::string functionName = antok::YAMLUtils::getString(function["Name"]);
	antok::Function* antokFunctionPtr = 0;
	if(functionName == "BeamNN") {
		// 		antokFunctionPtr = getUserFunction2(function, quantityNames, index);
		antokFunctionPtr = antok::user::hubers::generateGetNeuronalBeamEnergy(function, quantityNames, index);
	}
	return antokFunctionPtr;
}


antok::Function* antok::user::hubers::generateGetNeuronalBeamEnergy(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{
	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];


	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("X" , "double"));
	args.push_back(std::pair<std::string, std::string>("Y" , "double"));
	args.push_back(std::pair<std::string, std::string>("dX", "double"));
	args.push_back(std::pair<std::string, std::string>("dY", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* xAddr  = data.getAddr<double>(args[0].first);
	double* yAddr  = data.getAddr<double>(args[1].first);
	double* dxAddr = data.getAddr<double>(args[2].first);
	double* dyAddr = data.getAddr<double>(args[3].first);
	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetNeuronalBeamEnergy(xAddr,
				yAddr, dxAddr, dyAddr, data.getAddr<double>(quantityName)));
	return NULL;


}

antok::Function* getUserFunction2(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
	return 0;
}



