
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

antok::Function* antok::user::stefan::getCalcRICHPID(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
	if(quantityNames.size() != 6 && quantityNames.size() != 7) {
		std::cerr<<"Need 6/7 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return nullptr;
	}

	bool determine_pid = (quantityNames.size() == 7 )? true : false;

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args_per_index;
	std::vector<std::pair<std::string, std::string> > args;
	std::vector<std::pair<std::string, double*> > possible_const;
	args_per_index.push_back(std::pair<std::string, std::string>("PidLRichPion", "double"));
	args_per_index.push_back(std::pair<std::string, std::string>("PidLRichKaon", "double"));
	args_per_index.push_back(std::pair<std::string, std::string>("PidLRichProton", "double"));
	args_per_index.push_back(std::pair<std::string, std::string>("PidLRichElectron", "double"));
	args_per_index.push_back(std::pair<std::string, std::string>("PidLRichMuon", "double"));
	args_per_index.push_back(std::pair<std::string, std::string>("PidLRichBackground", "double"));

	if (determine_pid ){
		args_per_index.push_back(std::pair<std::string, std::string>("Mom", "TVector3"));
		possible_const.push_back(std::pair<std::string, double* >("PRatioCut", nullptr));
		possible_const.push_back(std::pair<std::string, double* >("MomPionMin", nullptr));
		possible_const.push_back(std::pair<std::string, double* >("MomPionMax", nullptr));
		possible_const.push_back(std::pair<std::string, double* >("MomKaonMin", nullptr));
		possible_const.push_back(std::pair<std::string, double* >("MomKaonMax", nullptr));
		possible_const.push_back(std::pair<std::string, double* >("MomProtonMin", nullptr));
		possible_const.push_back(std::pair<std::string, double* >("MomProtonMax", nullptr));
		possible_const.push_back(std::pair<std::string, double* >("MomElectronMin", nullptr));
		possible_const.push_back(std::pair<std::string, double* >("MomElectronMax", nullptr));
		possible_const.push_back(std::pair<std::string, double* >("MomMuonMin", nullptr));
		possible_const.push_back(std::pair<std::string, double* >("MomMuonMax", nullptr));

		for( auto& para: possible_const){
			if( antok::YAMLUtils::hasNodeKey( function, para.first )){
				para.second = antok::YAMLUtils::getAddress<double>(function[para.first]);
			} else{
				std::cerr << "Argument \"" << para.first << "\" not given for function \"" << function["Name"] << "\"." << std::endl;
				return nullptr;
			}
		}
	}

	if(not antok::generators::functionArgumentHandler(args, function, 0)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}
	if(not antok::generators::functionArgumentHandler(args_per_index, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}



	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size() - ((determine_pid)? 1: 0); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	int* pidAddrs = nullptr;
	if (determine_pid ){
		const size_t i = quantityNames.size() - 1;
		if(not data.insert<int>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		pidAddrs = data.getAddr<int>(quantityNames[i]);
	}

	if( ! determine_pid ){
	return new antok::user::stefan::functions::CalcRICHProbabilities(
			                                                data.getAddr<double>(args_per_index[0].first),
															data.getAddr<double>(args_per_index[1].first),
															data.getAddr<double>(args_per_index[2].first),
															data.getAddr<double>(args_per_index[3].first),
															data.getAddr<double>(args_per_index[4].first),
															data.getAddr<double>(args_per_index[5].first),
															quantityAddrs[0],
															quantityAddrs[1],
															quantityAddrs[2],
															quantityAddrs[3],
															quantityAddrs[4],
															quantityAddrs[5]
														);
	} else{
	return new antok::user::stefan::functions::CalcRICHPID(
			                                                data.getAddr<double>(args_per_index[0].first),
															data.getAddr<double>(args_per_index[1].first),
															data.getAddr<double>(args_per_index[2].first),
															data.getAddr<double>(args_per_index[3].first),
															data.getAddr<double>(args_per_index[4].first),
															data.getAddr<double>(args_per_index[5].first),
			                                                data.getAddr<TVector3>(args_per_index[6].first),
															possible_const[0].second,
															possible_const[1].second,
															possible_const[2].second,
															possible_const[3].second,
															possible_const[4].second,
															possible_const[5].second,
															possible_const[6].second,
															possible_const[7].second,
															possible_const[8].second,
															possible_const[9].second,
															possible_const[10].second,
															quantityAddrs[0],
															quantityAddrs[1],
															quantityAddrs[2],
															quantityAddrs[3],
															quantityAddrs[4],
															quantityAddrs[5],
															pidAddrs
														);

	}
}

antok::Function* antok::user::stefan::getDetermineKaonPionLV(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
	if(quantityNames.size() != 3) {
		std::cerr<<"Need 3 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return nullptr;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args_per_index;
	std::vector<std::pair<std::string, std::string> > args;
	std::vector<std::pair<std::string, double*> > possible_const;
	args.push_back(std::pair<std::string, std::string>("MomCandidate1", "TVector3"));
	args.push_back(std::pair<std::string, std::string>("PidCandidate1", "int"));
	args.push_back(std::pair<std::string, std::string>("MomCandidate2", "TVector3"));
	args.push_back(std::pair<std::string, std::string>("PidCandidate2", "int"));
	possible_const.push_back(std::pair<std::string, double*>("MassChargedKaon", nullptr));
	possible_const.push_back(std::pair<std::string, double*>("MassChargedPion", nullptr));

	if(not antok::generators::functionArgumentHandler(args, function, 0)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}
	if(not antok::generators::functionArgumentHandler(args_per_index, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	for( auto& para: possible_const){
		if( antok::YAMLUtils::hasNodeKey( function, para.first )){
			para.second = antok::YAMLUtils::getAddress<double>(function[para.first]);
		} else{
			std::cerr << "Argument \"" << para.first << "\" not given for function \"" << function["Name"] << "\"." << std::endl;
			return nullptr;
		}
	}

	data.insert<TLorentzVector>( quantityNames[0] );
	data.insert<TLorentzVector>( quantityNames[1] );
	data.insert<int>( quantityNames[2] );
	TLorentzVector* kaon_lv = data.getAddr<TLorentzVector>(quantityNames[0]);
	TLorentzVector* pion_lv = data.getAddr<TLorentzVector>(quantityNames[1]);
	int* is_kp_pk = data.getAddr<int>(quantityNames[2]);


	return new antok::user::stefan::functions::DetermineKaonPionLV(
																	data.getAddr<TVector3>( args[0].first ),
																	data.getAddr<int>(      args[1].first ),
																	data.getAddr<TVector3>( args[2].first ),
																	data.getAddr<int>(      args[3].first ),
																	possible_const[0].second,
																	possible_const[1].second,
																	kaon_lv,
																	pion_lv,
																	is_kp_pk
																	);
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
	if(functionName == "calcRICHPID") {
		antokFunctionPtr = stefan::getCalcRICHPID(function, quantityNames, index);
	}
	if(functionName == "determineKaonPionLV") {
		antokFunctionPtr = stefan::getDetermineKaonPionLV(function, quantityNames, index);
	}
	return antokFunctionPtr;
}
