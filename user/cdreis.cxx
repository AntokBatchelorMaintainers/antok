
#include<cdreis.h>

#include<constants.h>
#include<data.h>
#include<functions.hpp>
#include<generators_functions.h>
#include<cdreis_functions.hpp>
#include<yaml_utils.hpp>

antok::Function *antok::user::cdreis::getUserFunction(const YAML::Node &function,
                                                      std::vector<std::string> &quantityNames,
                                                      int index) {
	std::string functionName = antok::YAMLUtils::getString(function["Name"]);
	antok::Function *antokFunctionPtr = 0;
	if (functionName == "GetRecoilLorentzVec")
		antokFunctionPtr = antok::user::cdreis::generateGetRecoilLorentzVec(function, quantityNames, index);
	else if (functionName == "getPhotonLorentzVecs")
		antokFunctionPtr = antok::user::cdreis::generateGetPhotonLorentzVecs(function, quantityNames, index);
	else if (functionName == "getVectorLorentzVectorAttributes")
		antokFunctionPtr = antok::user::cdreis::generateGetVectorLorentzVectorAttributes(function, quantityNames, index);
	else if (functionName == "getPi0s")
		antokFunctionPtr = antok::user::cdreis::generateGetPi0s(function, quantityNames, index);
	else if (functionName == "getCleanedEcalClusters")
		antokFunctionPtr = antok::user::cdreis::generateGetCleanedEcalClusters(function, quantityNames, index);
	else if (functionName == "getPi0Pair")
		antokFunctionPtr = antok::user::cdreis::generateGetPi0Pair(function, quantityNames, index);
	else if(functionName == "getOmega")
		antokFunctionPtr = antok::user::cdreis::generateGetOmega(function, quantityNames, index);

	return antokFunctionPtr;
}

antok::Function *
antok::user::cdreis::generateGetRecoilLorentzVec(const YAML::Node &function, std::vector<std::string> &quantityNames,
                                                 int index) {

	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}
	try {
		function["RecoilMass"].as<double>();
	} catch (const YAML::TypedBadConversion<double> &e) {
		std::cerr << "Argument \"RecoilMass\" in function \"GetRecoilLorentzVec\" should be of type double (variable \""
		          << quantityName << "\")." << std::endl;
		return 0;
	}
	double *RecoilMass = new double();
	(*RecoilMass) = function["RecoilMass"].as<double>();

	antok::Data &data = antok::ObjectManager::instance()->getData();

	TLorentzVector *BeamLorentzVec = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector *XLorentzVec = data.getAddr<TLorentzVector>(args[1].first);

	if (not data.insert<TLorentzVector>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::cdreis::functions::GetRecoilLorentzVec(BeamLorentzVec, XLorentzVec, RecoilMass,
	                                                                data.getAddr<TLorentzVector>(quantityName)));
};

antok::Function *
antok::user::cdreis::generateGetPhotonLorentzVecs(const YAML::Node &function, std::vector<std::string> &quantityNames,
                                                  int index) {
	if (quantityNames.size() > 2) {
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("VectorX", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorY", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorZ", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorE", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorTime", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("xPV", "double"));
	args.push_back(std::pair<std::string, std::string>("yPV", "double"));
	args.push_back(std::pair<std::string, std::string>("zPV", "double"));

	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<double> *VectorXAddr = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double> *VectorYAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double> *VectorZAddr = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double> *VectorEAddr = data.getAddr<std::vector<double> >(args[3].first);
	std::vector<double> *VectorTimeAddr = data.getAddr<std::vector<double> >(args[4].first);
	double *xPVAddr = data.getAddr<double>(args[5].first);
	double *yPVAddr = data.getAddr<double>(args[6].first);
	double *zPVAddr = data.getAddr<double>(args[7].first);

	std::string resultVec = quantityNames[0];
	std::string resultECALIndex = quantityNames[1];

	if (not data.insert<std::vector<TLorentzVector> >(resultVec)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[0]);
		return 0;
	}

	if (not data.insert<std::vector<int> >(resultECALIndex)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[1]);
		return 0;
	}

	return (new antok::user::cdreis::functions::GetPhotonLorentzVecs(VectorXAddr,
	                                                                 VectorYAddr,
	                                                                 VectorZAddr,
	                                                                 VectorEAddr,
	                                                                 VectorTimeAddr,
	                                                                 xPVAddr,
	                                                                 yPVAddr,
	                                                                 zPVAddr,
	                                                                 data.getAddr<std::vector<TLorentzVector> >(resultVec),
	                                                                 data.getAddr<std::vector<int> >(resultECALIndex))
	);
};


antok::Function *
antok::user::cdreis::generateGetCleanedEcalClusters(const YAML::Node &function, std::vector<std::string> &quantityNames,
                                                    int index) {
	if (quantityNames.size() > 6) {
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("VectorX", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorY", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorZ", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorE", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorT", "std::vector<double>"));

	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<double> *VectorXAddr = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double> *VectorYAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double> *VectorZAddr = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double> *VectorEAddr = data.getAddr<std::vector<double> >(args[3].first);
	std::vector<double> *VectorTAddr = data.getAddr<std::vector<double> >(args[4].first);

	std::string resultVectorX = quantityNames[0];
	std::string resultVectorY = quantityNames[1];
	std::string resultVectorZ = quantityNames[2];
	std::string resultVectorE = quantityNames[3];
	std::string resultVectorT = quantityNames[4];
	std::string resultVectorIndex = quantityNames[5];

	for (unsigned int i = 0; i < quantityNames.size() - 1; i++) {
		if (not data.insert<std::vector<double> >(quantityNames[i])) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return 0;
		}
	}
	if (not data.insert<std::vector<int> >(quantityNames[5])) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[5]);
		return 0;
	}

	return (new antok::user::cdreis::functions::GetCleanedEcalClusters(VectorXAddr,
	                                                                   VectorYAddr,
	                                                                   VectorZAddr,
	                                                                   VectorEAddr,
	                                                                   VectorTAddr,
	                                                                   data.getAddr<std::vector<double> >(resultVectorX),
	                                                                   data.getAddr<std::vector<double> >(resultVectorY),
	                                                                   data.getAddr<std::vector<double> >(resultVectorZ),
	                                                                   data.getAddr<std::vector<double> >(resultVectorE),
	                                                                   data.getAddr<std::vector<double> >(resultVectorT),
	                                                                   data.getAddr<std::vector<int> >(resultVectorIndex))
	);
};


antok::Function *antok::user::cdreis::generateGetVectorLorentzVectorAttributes(const YAML::Node &function,
                                                                               std::vector<std::string> &quantityNames,
                                                                               int index) {
	if (quantityNames.size() > 7) {
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("VectorLV", "std::vector<TLorentzVector>"));

	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<TLorentzVector> *VectorLVAddr = data.getAddr<std::vector<TLorentzVector> >(args[0].first);

	std::string resultVecX = quantityNames[0];
	std::string resultVecY = quantityNames[1];
	std::string resultVecZ = quantityNames[2];
	std::string resultVecE = quantityNames[3];
	std::string resultVecTheta = quantityNames[4];
	std::string resultVecPhi = quantityNames[5];
	std::string resultVecMag = quantityNames[6];

	for (unsigned int i = 0; i < quantityNames.size(); i++) {
		if (not data.insert<std::vector<double> >(quantityNames[i])) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return 0;
		}
	}

	return (new antok::user::cdreis::functions::GetVectorLorentzVectorAttributes(VectorLVAddr,
	                                                                             data.getAddr<std::vector<double> >(resultVecX),
	                                                                             data.getAddr<std::vector<double> >(resultVecY),
	                                                                             data.getAddr<std::vector<double> >(resultVecZ),
	                                                                             data.getAddr<std::vector<double> >(resultVecE),
	                                                                             data.getAddr<std::vector<double> >(resultVecTheta),
	                                                                             data.getAddr<std::vector<double> >(resultVecPhi),
	                                                                             data.getAddr<std::vector<double> >(resultVecMag))
	);
};

antok::Function *
antok::user::cdreis::generateGetPi0s(const YAML::Node &function, std::vector<std::string> &quantityNames, int index) {
	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("VectorLV", "std::vector<TLorentzVector>"));
	args.push_back(std::pair<std::string, std::string>("ECALIndex", "std::vector<int>"));

	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<TLorentzVector> *VectorLVAddr = data.getAddr<std::vector<TLorentzVector> >(args[0].first);
	std::vector<int> *ECALIndex = data.getAddr<std::vector<int> >(args[1].first);

	if (not antok::YAMLUtils::hasNodeKey(function, "Mass")) {
		std::cerr << "Argument \"" << "mass" << "\" not found (required for function \"" << function["Name"] << "\")."
		          << std::endl;
		return 0;
	}
	double *massAddr = antok::YAMLUtils::getAddress<double>(function["Mass"]);

	std::string resultVecLV = quantityNames[0];

	if (not data.insert<std::vector<TLorentzVector> >(resultVecLV)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultVecLV);
		return 0;
	}

	return (new antok::user::cdreis::functions::GetPi0s(VectorLVAddr,
	                                                    ECALIndex,
	                                                    massAddr,
	                                                    data.getAddr<std::vector<TLorentzVector> >(resultVecLV))

	);
};

antok::Function *
antok::user::cdreis::generateGetPi0Pair(const YAML::Node &function, std::vector<std::string> &quantityNames,
                                        int index) {
	if (quantityNames.size() > 4) {
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("VectorLV", "std::vector<TLorentzVector>"));

	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<TLorentzVector> *VectorLVAddr = data.getAddr<std::vector<TLorentzVector> >(args[0].first);

	std::string resultVecLV = quantityNames[0];
	std::string resultVecLV0 = quantityNames[1];
	std::string resultVecLV1 = quantityNames[2];
	std::string resultGoodPair = quantityNames[3];

	if (not data.insert<std::vector<TLorentzVector> >(resultVecLV)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultVecLV);
		return 0;
	}

	if (not data.insert<TLorentzVector>(resultVecLV0)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultVecLV0);
		return 0;
	}

	if (not data.insert<TLorentzVector>(resultVecLV1)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultVecLV1);
		return 0;
	}

	if (not data.insert<int>(resultGoodPair)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultGoodPair);
		return 0;
	}

	return (new antok::user::cdreis::functions::GetPi0Pair(VectorLVAddr,
	                                                       data.getAddr<std::vector<TLorentzVector> >(resultVecLV),
	                                                       data.getAddr<TLorentzVector>(resultVecLV0),
	                                                       data.getAddr<TLorentzVector>(resultVecLV1),
	                                                       data.getAddr<int>(resultGoodPair))
	);
};

antok::Function* antok::user::cdreis::generateGetOmega(const YAML::Node& function, std::vector<std::string>& quantityNames, int index){
	if(quantityNames.size() > 2) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Pi0_0", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("Pi0_1", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("Scattered0", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("Scattered1", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("Scattered2", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("Charge0", "int"));
	args.push_back(std::pair<std::string, std::string>("Charge1", "int"));
	args.push_back(std::pair<std::string, std::string>("Charge2", "int"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	TLorentzVector* Pi0_0Addr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* Pi0_1Addr = data.getAddr<TLorentzVector>(args[1].first);
	TLorentzVector* Scattered0Addr = data.getAddr<TLorentzVector>(args[2].first);
	TLorentzVector* Scattered1Addr = data.getAddr<TLorentzVector>(args[3].first);
	TLorentzVector* Scattered2Addr = data.getAddr<TLorentzVector>(args[4].first);

	int* Charged0Addr = data.getAddr<int>(args[5].first);
	int* Charged1Addr = data.getAddr<int>(args[6].first);
	int* Charged2Addr = data.getAddr<int>(args[7].first);

	std::string resultOmegaLV = quantityNames[0];
	std::string resultAccepted = quantityNames[1];

	if(not data.insert<TLorentzVector>(quantityNames[0])) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[0]);
		return 0;
	}

	if(not data.insert<int>(quantityNames[1])) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[1]);
		return 0;
	}

	return (new antok::user::cdreis::functions::GetOmega( Pi0_0Addr,
	                                                      Pi0_1Addr,
	                                                      Scattered0Addr,
	                                                      Scattered1Addr,
	                                                      Scattered2Addr,
	                                                      Charged0Addr,
	                                                      Charged1Addr,
	                                                      Charged2Addr,
	                                                      data.getAddr<TLorentzVector>(quantityNames[0]),
	                                                      data.getAddr<int>(quantityNames[1]) )
	);
};

antok::Function * antok::user::cdreis::generateGetECALCorrectedTiming(const YAML::Node &function, std::vector<std::string> &quantityNames,
                                                                      int index) {

	if (quantityNames.size() > 1) {
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Timing", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("Energy", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("ECALIndex", "std::vector<int>"));

	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<double>* Timing = data.getAddr<std::vector<double>>(args[0].first);
	std::vector<double>* Energy = data.getAddr<std::vector<double>>(args[1].first);
	std::vector<int>* ECALIndex = data.getAddr<std::vector<int>>(args[2].first);

	if (not data.insert<std::vector<double>>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::cdreis::functions::GetECALCorrectedTiming(Timing,
	                                                                   Energy,
	                                                                   ECALIndex,
	                                                                   data.getAddr<std::vector<double>>(quantityName)));
};

