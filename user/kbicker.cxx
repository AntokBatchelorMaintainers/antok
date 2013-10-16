
#include<kbicker.h>

#include<data.h>
#include<functions.hpp>
#include<generators_functions.h>
#include<kbicker_functions.hpp>
#include<yaml_utils.hpp>

antok::Function* antok::user::kbicker::getUserFunction(const YAML::Node& function,
                                                       std::vector<std::string>& quantityNames,
                                                       int index)
{
	std::string functionName = antok::YAMLUtils::getString(function["Name"]);
	Function* antokFunctionPtr = 0;
	if(functionName == "getRpdExpectedHitsParameters") {
		antokFunctionPtr = antok::user::kbicker::generateGetRpdExpectedHitsParameters(function, quantityNames, index);
	} else if(functionName == "getRpdPhi") {
		antokFunctionPtr = antok::user::kbicker::generateGetRpdPhi(function, quantityNames, index);
	}
	return antokFunctionPtr;
}

antok::Function* antok::user::kbicker::generateGetRpdExpectedHitsParameters(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {

	using antok::YAMLUtils::hasNodeKey;

	if(quantityNames.size() != 4) {
		std::cerr<<"Need 4 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("Vertex", "TVector3"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* beamLorentzVecAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* XLorentzVecAddr = data.getAddr<TLorentzVector>(args[1].first);
	TVector3* vertexAddr = data.getAddr<TVector3>(args[2].first);

	std::vector<std::string> possiblyConstArgs;
	possiblyConstArgs.push_back("XOffset");
	possiblyConstArgs.push_back("YOffset");
	possiblyConstArgs.push_back("XAngle");
	possiblyConstArgs.push_back("YAngle");

	for(unsigned int i = 0; i < possiblyConstArgs.size(); ++i) {
		if(not hasNodeKey(function, possiblyConstArgs[i])) {
			std::cerr<<"Argument \""<<possiblyConstArgs[i]<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
	}

	double* xOffsetAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[0]]);
	double* yOffsetAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[1]]);
	double* xAngleAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[2]]);
	double* yAngleAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[3]]);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return (new antok::user::kbicker::functions::GetRPDExpectedHitsParameters(beamLorentzVecAddr,
	                                                                          XLorentzVecAddr,
	                                                                          vertexAddr,
	                                                                          xOffsetAddr,
	                                                                          yOffsetAddr,
	                                                                          xAngleAddr,
	                                                                          yAngleAddr,
	                                                                          quantityAddrs[0],
	                                                                          quantityAddrs[1],
	                                                                          quantityAddrs[2],
	                                                                          quantityAddrs[3]));

}

antok::Function* antok::user::kbicker::generateGetRpdPhi(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(not (quantityNames.size() == 2 or quantityNames.size() == 4)) {
		std::cerr<<"Need 2 names for function \"getRpdPhi\""<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::string method = antok::YAMLUtils::getString(function["Method"]);
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("RPDProtonLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* beamLVAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* RPDprotLVAddr = data.getAddr<TLorentzVector>(args[1].first);
	TLorentzVector* xLVAddr = data.getAddr<TLorentzVector>(args[2].first);

	TVector3* vertexAddr = 0;

	if(method == "Prediction") {
		args.clear();
		args.push_back(std::pair<std::string, std::string>("Vertex", "TVector3"));
		if(not antok::generators::functionArgumentHandler(args, function, index)) {
			std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
			std::cerr<<"\"Prediction\" method requires \"Vertex\"."<<std::endl;
			return 0;
		}
		vertexAddr = data.getAddr<TVector3>(args[0].first);
	}


	int methodSwitch = -1;
	if(method == "Projection") {
		if(quantityNames.size() != 2) {
			std::cerr<<"\"getRpdPhi\" with method \"Projection\" only works with 2 calculated quantities, found "
			         <<quantityNames.size()<<" when calculating variables \"["<<std::endl;
			for(unsigned int i = 0; i < quantityNames.size()-1; ++i) {
				std::cerr<<quantityNames[i]<<", ";
			}
			std::cerr<<quantityNames[quantityNames.size()-1]<<"]\"."<<std::endl;
			return 0;
		}
		methodSwitch = 0;
	} else if (method == "Rotation") {
		methodSwitch = 1;
	} else if (method == "Prediction") {
		methodSwitch = 2;
	} else {
		if(method == "") {
			std::cerr<<"Could not convert \"Method\" to std::string in function \"getRpdPhi\" when calculating variables \"["<<std::endl;
		} else {
			std::cerr<<"Method \""<<method<<"\" is not supported for function \"getRpdPhi\" when calculating variables \"["<<std::endl;
		}
		for(unsigned int i = 0; i < quantityNames.size()-1; ++i) {
			std::cerr<<quantityNames[i]<<", ";
		}
		std::cerr<<quantityNames[quantityNames.size()-1]<<"]\"."<<std::endl;
		return 0;
	}

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	double* quantityAddrs2 = 0;
	double* quantityAddrs3 = 0;

	if(quantityNames.size() == 4) {
		quantityAddrs2 = quantityAddrs[2];
		quantityAddrs3 = quantityAddrs[3];
	}
	return (new antok::user::kbicker::functions::GetRpdPhi(beamLVAddr,
	                                                       RPDprotLVAddr,
	                                                       xLVAddr,
	                                                       quantityAddrs[0],
	                                                       quantityAddrs[1],
	                                                       methodSwitch,
	                                                       quantityAddrs2,
	                                                       quantityAddrs3,
	                                                       vertexAddr));

};
