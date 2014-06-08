
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
	if(functionName == "sqrt")
		antokFunctionPtr = antok::user::hubers::generateSqrt(function, quantityNames, index);
	else if(functionName == "frac")
		antokFunctionPtr = antok::user::hubers::generateFrac(function, quantityNames, index);
	else if(functionName == "getPt")
		antokFunctionPtr = antok::user::hubers::generateGetPt(function, quantityNames, index);
	else if(functionName == "EnforceEConservation")
		antokFunctionPtr = antok::user::hubers::generateEnforceEConservation(function, quantityNames, index);
	else if(functionName == "BeamNN")
		antokFunctionPtr = antok::user::hubers::generateGetNeuronalBeam(function, quantityNames, index);
	else if(functionName == "Theta")
		antokFunctionPtr = antok::user::hubers::generateGetTheta(function, quantityNames, index);
	else if(functionName == "ThetaZ")
		antokFunctionPtr = antok::user::hubers::generateGetThetaZCut(function, quantityNames, index);
	return antokFunctionPtr;
}

antok::Function* antok::user::hubers::generateSqrt(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Arg", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* argAddr = data.getAddr<double>(args[0].first);

	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::Sqrt(argAddr, data.getAddr<double>(quantityName)));
};

antok::Function* antok::user::hubers::generateFrac(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Numerator", "double"));
	args.push_back(std::pair<std::string, std::string>("Denominator", "double"));

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

	return (new antok::user::hubers::functions::Frac(arg1Addr, arg2Addr, data.getAddr<double>(quantityName)));
};

antok::Function* antok::user::hubers::generateGetPt(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 1) {
		std::cerr<<"Need 1 names for function \"getTs\""<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;

	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
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

	return (new antok::user::hubers::functions::GetPt(xLVAddr, beamLVAddr, quantityAddrs[0]));
};

antok::Function* antok::user::hubers::generateEnforceEConservation(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("LVBeam", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("LVPion", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("LVGamma", "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	TLorentzVector* beamAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* pionAddr = data.getAddr<TLorentzVector>(args[1].first);
	TLorentzVector* gammaAddr = data.getAddr<TLorentzVector>(args[2].first);
	if(not data.insert<TLorentzVector>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::EnforceEConservation<TLorentzVector>(beamAddr, pionAddr, gammaAddr, data.getAddr<TLorentzVector>(quantityName)));
}



antok::Function* antok::user::hubers::generateGetNeuronalBeam(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 2) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<" needed two\"."<<std::endl;
		return 0;
	}
	std::string quantityNameD = quantityNames[0];
	std::string quantityNameLV = quantityNames[1];


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
	if(not data.insert<double>(quantityNameD)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}
	if(not data.insert<TLorentzVector>(quantityNameLV)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetNeuronalBeam(xAddr, yAddr, dxAddr, dyAddr, data.getAddr<double>(quantityNameD), data.getAddr<TLorentzVector>(quantityNameLV)));
}


//***********************************
//Calculates the angle theta between
//two TLorentzVectors
//***********************************
antok::Function* antok::user::hubers::generateGetTheta(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];


	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("OutLorentzVec",  "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	TLorentzVector* beamLVAddr  = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* outLVAddr   = data.getAddr<TLorentzVector>(args[1].first);
	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetTheta(beamLVAddr, outLVAddr, data.getAddr<double>(quantityName)));
}

//***********************************
//Calculates the condition for a
//theta dependend Z cut
//***********************************
antok::Function* antok::user::hubers::generateGetThetaZCut(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];


	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Theta" , "double"));
	args.push_back(std::pair<std::string, std::string>("Z" , "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* thetaAddr  = data.getAddr<double>(args[0].first);
	double* zAddr  = data.getAddr<double>(args[1].first);

	double *zMeanAddr;
	try {
		function["ZMean"].as<double>();
	} catch(const YAML::TypedBadConversion<double>& e) {
		std::cerr<<"Argument \"ZMean\" in function \""<<__func__<<"\" should be of type double (variable \""<<"ZMean\")."<<std::endl;
		return 0;
	}
	zMeanAddr = new double();
	(*zMeanAddr) = function["ZMean"].as<double>();


	if(not data.insert<int>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetThetaZCut(zAddr, thetaAddr, zMeanAddr, data.getAddr<int>(quantityName)));
}
