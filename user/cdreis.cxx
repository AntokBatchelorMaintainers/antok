
#include<cdreis.h>

#include<constants.h>
#include<data.h>
#include<functions.hpp>
#include<generators_functions.h>
#include<cdreis_functions.hpp>
#include<yaml_utils.hpp>
#include<iostream>

antok::Function* antok::user::cdreis::getUserFunction(const YAML::Node& function,
                                                      std::vector<std::string>& quantityNames,
                                                      int index)
{
	std::string functionName = antok::YAMLUtils::getString(function["Name"]);
	antok::Function* antokFunctionPtr = 0;
	if(functionName == "GetRecoilLorentzVec")
		antokFunctionPtr =  antok::user::cdreis::generateGetRecoilLorentzVec(function, quantityNames, index);
	return antokFunctionPtr;
}

antok::Function* antok::user::cdreis::generateGetRecoilLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}
	try {
		function["RecoilMass"].as<double>();
	} catch(const YAML::TypedBadConversion<double>& e) {
		std::cerr<<"Argument \"RecoilMass\" in function \"GetRecoilLorentzVec\" should be of type double (variable \""<<quantityName<<"\")."<<std::endl;
		return 0;
	}
	double* RecoilMass = new double();
	(*RecoilMass) = function["RecoilMass"].as<double>();

	antok::Data& data = antok::ObjectManager::instance()->getData();

	TLorentzVector* BeamLorentzVec = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* XLorentzVec = data.getAddr<TLorentzVector>(args[1].first);

	if(not data.insert<TLorentzVector>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::cdreis::functions::GetRecoilLorentzVec(BeamLorentzVec, XLorentzVec, RecoilMass, data.getAddr<TLorentzVector>(quantityName)));
};
