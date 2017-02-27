

#include<iostream>
#include<sstream>
#include<string>
#include<functions.hpp>
#include<generators_functions.h>
#include<data.h>
#include<yaml_utils.hpp>

#include<swallner.h>
#include<swallner_functions.hpp>

namespace {
	std::vector<double> getCalcCEDARPIDGetThresholds(const char* name, const YAML::Node& function);
}


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
	using antok::YAMLUtils::hasNodeKey;

	bool determine_pid = (quantityNames.size() == 7 )? true : false;

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args_per_index;
	std::vector<std::pair<std::string, std::string> > args;
	std::vector<std::pair<std::string, double*> > possible_const;
	std::vector<std::pair<std::string, double*> > possible_const_per_index;
	bool is_momv3_given = not hasNodeKey(function, "MomMag");

	// complete list of arguments
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichPion", nullptr));
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichKaon", nullptr));
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichProton", nullptr));
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichElectron", nullptr));
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichMuon", nullptr));
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichBackground", nullptr));


	if (determine_pid ){
		if( is_momv3_given )
			args_per_index.push_back(std::pair<std::string, std::string>("Mom", "TVector3"));
		else
			args_per_index.push_back(std::pair<std::string, std::string>("MomMag", "double"));
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

		if( not antok::generators::functionrgumentHandlerPossibleConst<double>(possible_const, function, 0 ) ){
			std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
			return 0;
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

	if(not antok::generators::functionrgumentHandlerPossibleConst<double>(possible_const_per_index, function, index ) ){
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
															possible_const_per_index[0].second,
															possible_const_per_index[1].second,
															possible_const_per_index[2].second,
															possible_const_per_index[3].second,
															possible_const_per_index[4].second,
															possible_const_per_index[5].second,
															quantityAddrs[0],
															quantityAddrs[1],
															quantityAddrs[2],
															quantityAddrs[3],
															quantityAddrs[4],
															quantityAddrs[5]
														);
	} else{
	return new antok::user::stefan::functions::CalcRICHPID(
															possible_const_per_index[0].second,
															possible_const_per_index[1].second,
															possible_const_per_index[2].second,
															possible_const_per_index[3].second,
															possible_const_per_index[4].second,
															possible_const_per_index[5].second,
			                                                (is_momv3_given)? data.getAddr<TVector3>(args_per_index[0].first): NULL,
			                                                (is_momv3_given)? NULL : data.getAddr<double>(args_per_index[0].first),
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
	using antok::YAMLUtils::hasNodeKey;


	if(quantityNames.size() != 5 and quantityNames.size() != 7) {
		std::cerr<<"Need 5/7 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return nullptr;
	}


	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args_per_index;
	std::vector<std::pair<std::string, std::string> > args;
	std::vector<std::pair<std::string, double*> > possible_const;
	std::vector<std::pair<std::string, int*> > possible_const_int;
	args.push_back(std::pair<std::string, std::string>("MomCandidate1", "TVector3"));
	args.push_back(std::pair<std::string, std::string>("PidCandidate1", "int"));
	args.push_back(std::pair<std::string, std::string>("MomCandidate2", "TVector3"));
	args.push_back(std::pair<std::string, std::string>("PidCandidate2", "int"));
	possible_const.push_back(std::pair<std::string, double*>("MassChargedKaon", nullptr));
	possible_const.push_back(std::pair<std::string, double*>("MassChargedPion", nullptr));

	possible_const_int.push_back(std::pair<std::string, int*>("Method", nullptr));
	if ( hasNodeKey(function, "PidCandidate1Mct")  or hasNodeKey(function, "PidCandidate2Mct")){
		if ( not (hasNodeKey(function, "PidCandidate1Mct") and hasNodeKey(function, "PidCandidate2Mct") ) ){
			std::cerr << "Pid MCT value only given for one of both candidates in function " << function["Name"] << std::endl;
			return 0;
		}
		if(quantityNames.size() != 7) {
			std::cerr << "Pid MCT values given, but no 7 return values given in function " << function["Name"] << std::endl;
			return 0;
		}
		possible_const_int.push_back(std::pair<std::string, int*>("PidCandidate1Mct", nullptr));
		possible_const_int.push_back(std::pair<std::string, int*>("PidCandidate2Mct", nullptr));
	} else {
		if(quantityNames.size() != 5) {
			std::cerr << "No Pid MCT values given, but no 5 return values given in function " << function["Name"] << std::endl;
			return 0;
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

	if( not antok::generators::functionrgumentHandlerPossibleConst<double>(possible_const, function, 0) ){
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}
	if( not antok::generators::functionrgumentHandlerPossibleConst<int>(possible_const_int, function, 0) ){
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}


	if (possible_const_int.size() == 1){ // add dummy entries if no Mct values are given
		possible_const_int.push_back(std::pair<std::string, int*>("PidCandidate1Mct", nullptr));
		possible_const_int.push_back(std::pair<std::string, int*>("PidCandidate2Mct", nullptr));
	}


	data.insert<TLorentzVector>( quantityNames[0] );
	data.insert<TLorentzVector>( quantityNames[1] );
	data.insert<int>( quantityNames[2] );
	data.insert<int>( quantityNames[3] );
	data.insert<int>( quantityNames[4] );
	TLorentzVector* kaon_lv = data.getAddr<TLorentzVector>(quantityNames[0]);
	TLorentzVector* pion_lv = data.getAddr<TLorentzVector>(quantityNames[1]);
	int* is_kp_pk = data.getAddr<int>(quantityNames[2]);
	int* pid_kaon = data.getAddr<int>(quantityNames[3]);
	int* pid_pion = data.getAddr<int>(quantityNames[4]);

	int* pid_kaon_mct = nullptr;
	int* pid_pion_mct = nullptr;
	if(quantityNames.size() == 7){
		data.insert<int>( quantityNames[5] );
		data.insert<int>( quantityNames[6] );
		pid_kaon_mct = data.getAddr<int>(quantityNames[5]);
		pid_pion_mct = data.getAddr<int>(quantityNames[6]);
	}


	return new antok::user::stefan::functions::DetermineKaonPionLV(
																	data.getAddr<TVector3>( args[0].first ),
																	data.getAddr<int>(      args[1].first ),
																	data.getAddr<TVector3>( args[2].first ),
																	data.getAddr<int>(      args[3].first ),
																	possible_const[0].second,
																	possible_const[1].second,
																	kaon_lv,
																	pion_lv,
																	is_kp_pk,
																	pid_kaon,
																	pid_pion,
																	pid_kaon_mct,
																	pid_pion_mct,
																	possible_const_int[0].second,
																	possible_const_int[1].second,
																	possible_const_int[2].second
																	);
}

antok::Function* antok::user::stefan::getDetermineKaonPionLVLikelihood(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
	if(quantityNames.size() != 5) {
		std::cerr<<"Need 5 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return nullptr;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args_per_index;
	std::vector<std::pair<std::string, std::string> > args;
	std::vector<std::pair<std::string, double*> > possible_const;
	std::vector<std::pair<std::string, int*> > possible_const_int;
	args.push_back(std::pair<std::string, std::string>("MomCandidate1", "TVector3"));
	args.push_back(std::pair<std::string, std::string>("LPionCandidate1", "double"));
	args.push_back(std::pair<std::string, std::string>("LKaonCandidate1", "double"));
	args.push_back(std::pair<std::string, std::string>("PidCandidate1", "int"));
	args.push_back(std::pair<std::string, std::string>("MomCandidate2", "TVector3"));
	args.push_back(std::pair<std::string, std::string>("LPionCandidate2", "double"));
	args.push_back(std::pair<std::string, std::string>("LKaonCandidate2", "double"));
	args.push_back(std::pair<std::string, std::string>("PidCandidate2", "int"));
	possible_const.push_back(std::pair<std::string, double*>("MassChargedKaon", nullptr));
	possible_const.push_back(std::pair<std::string, double*>("MassChargedPion", nullptr));
	possible_const.push_back(std::pair<std::string, double*>("ThresholdLogKPi", nullptr));
	possible_const.push_back(std::pair<std::string, double*>("ThresholdLogPiK", nullptr));


	if(not antok::generators::functionArgumentHandler(args, function, 0)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}
	if(not antok::generators::functionArgumentHandler(args_per_index, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	if( not antok::generators::functionrgumentHandlerPossibleConst<double>(possible_const, function, 0) ){
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}
	if( not antok::generators::functionrgumentHandlerPossibleConst<int>(possible_const_int, function, 0) ){
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}


	data.insert<TLorentzVector>( quantityNames[0] );
	data.insert<TLorentzVector>( quantityNames[1] );
	data.insert<int>( quantityNames[2] );
	data.insert<int>( quantityNames[3] );
	data.insert<int>( quantityNames[4] );
	TLorentzVector* kaon_lv = data.getAddr<TLorentzVector>(quantityNames[0]);
	TLorentzVector* pion_lv = data.getAddr<TLorentzVector>(quantityNames[1]);
	int* is_kp_pk = data.getAddr<int>(quantityNames[2]);
	int* pid_kaon = data.getAddr<int>(quantityNames[3]);
	int* pid_pion = data.getAddr<int>(quantityNames[4]);


	return new antok::user::stefan::functions::DetermineKaonPionLVLikelihood(
																	data.getAddr<TVector3>( args[0].first ),
																	data.getAddr<double>(   args[1].first),
																	data.getAddr<double>(   args[2].first),
																	data.getAddr<int>(      args[3].first ),
																	data.getAddr<TVector3>( args[4].first ),
																	data.getAddr<double>(   args[5].first),
																	data.getAddr<double>(   args[6].first),
																	data.getAddr<int>(      args[7].first ),
																	possible_const[0].second,
																	possible_const[1].second,
																	possible_const[2].second,
																	possible_const[3].second,
																	kaon_lv,
																	pion_lv,
																	is_kp_pk,
																	pid_kaon,
																	pid_pion
																	);
}

antok::Function* antok::user::stefan::getCalcCEDARPID(const YAML::Node& function, std::vector<std::string>& quantityNames, int index){
	using antok::YAMLUtils::hasNodeKey;
	if(hasNodeKey(function, "ThresholdsKaonDeltaLogLikeCedar1")){
		return antok::user::stefan::getCalcCEDARPIDMulitL(function, quantityNames, index);
	}else{
		return antok::user::stefan::getCalcCEDARPIDOneL(function, quantityNames, index);
	}
}

antok::Function* antok::user::stefan::getCalcCEDARPIDMulitL(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
	if(quantityNames.size() != 5 and quantityNames.size() != 3 ) {
		std::cerr<<"Need 3/5 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return nullptr;
	}
	using antok::YAMLUtils::hasNodeKey;

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, double*> > possible_const_double;
	std::vector<std::pair<std::string, int*> > possible_const_int;

	// complete list of arguments
	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodPionCedar1", nullptr));
	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodKaonCedar1", nullptr));
//	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodProtonCedar1", nullptr));
	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodPionCedar2", nullptr));
	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodKaonCedar2", nullptr));
//	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodProtonCedar2", nullptr));
	possible_const_int.push_back(std::pair<std::string, int* >("NHitsCedar1", nullptr));
	possible_const_int.push_back(std::pair<std::string, int* >("NHitsCedar2", nullptr));

	if( not antok::generators::functionrgumentHandlerPossibleConst<double>(possible_const_double, function, 0 ) ){
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	if( not antok::generators::functionrgumentHandlerPossibleConst<int>(possible_const_int, function, 0 ) ){
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	const std::vector<double> thresholds_kaon_CEDAR1 = getCalcCEDARPIDGetThresholds("ThresholdsKaonDeltaLogLikeCedar1", function);
	const std::vector<double> thresholds_pion_CEDAR1 = getCalcCEDARPIDGetThresholds("ThresholdsPionDeltaLogLikeCedar1", function);
	const std::vector<double> thresholds_kaon_CEDAR2 = getCalcCEDARPIDGetThresholds("ThresholdsKaonDeltaLogLikeCedar2", function);
	const std::vector<double> thresholds_pion_CEDAR2 = getCalcCEDARPIDGetThresholds("ThresholdsPionDeltaLogLikeCedar2", function);
	if( thresholds_kaon_CEDAR1.size() == 0 or thresholds_pion_CEDAR1.size() == 0 or thresholds_kaon_CEDAR2.size() == 0 or thresholds_pion_CEDAR2.size() == 0){
		return 0;
	}


	std::vector<int*> quantityAddrs_int;
	std::vector<double*> quantityAddrs_double;
	for(unsigned int i = 0; i < 3; ++i) {
		if(not data.insert<int>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs_int.push_back(data.getAddr<int>(quantityNames[i]));

	}
	if (quantityNames.size() == 5) {
		for (unsigned int i = 3; i < quantityNames.size(); ++i) {
			if (not data.insert<double>(quantityNames[i])) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
				return 0;
			}
			quantityAddrs_double.push_back(data.getAddr<double>(quantityNames[i]));
		}
	} else {
		for (unsigned int i = 3; i < quantityNames.size(); ++i) {
			quantityAddrs_double.push_back(nullptr);
		}

	}

	return new antok::user::stefan::functions::CalcCEDARPID(
	        possible_const_double[0].second,
	        possible_const_double[1].second,
	        new double(),
	        possible_const_int[0].second,
	        thresholds_kaon_CEDAR1,
	        thresholds_pion_CEDAR1,
	        possible_const_double[2].second,
	        possible_const_double[3].second,
	        new double(),
	        possible_const_int[1].second,
	        thresholds_kaon_CEDAR2,
	        thresholds_pion_CEDAR2,
	        quantityAddrs_int[0],
	        quantityAddrs_int[1],
	        quantityAddrs_int[2],
	        quantityAddrs_double[0],
	        quantityAddrs_double[1]
	        );
}


antok::Function* antok::user::stefan::getCalcCEDARPIDOneL(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
	if(quantityNames.size() != 4 and quantityNames.size() != 1) {
		std::cerr<<"Need 1/4 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return nullptr;
	}
	using antok::YAMLUtils::hasNodeKey;

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, double*> > possible_const_double;
	std::vector<std::pair<std::string, int*> > possible_const_int;

	// complete list of arguments
	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodPionCedar1", nullptr));
	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodKaonCedar1", nullptr));
//	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodProtonCedar1", nullptr));
	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodPionCedar2", nullptr));
	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodKaonCedar2", nullptr));
//	possible_const_double.push_back(std::pair<std::string, double* >("LikelihoodProtonCedar2", nullptr));
//	possible_const_int.push_back(std::pair<std::string, int* >("NHitsCedar1", nullptr));
//	possible_const_int.push_back(std::pair<std::string, int* >("NHitsCedar2", nullptr));
	possible_const_double.push_back(std::pair<std::string, double* >("ThresholdKaonDeltaLogLike", nullptr));
	possible_const_double.push_back(std::pair<std::string, double* >("ThresholdPionDeltaLogLike", nullptr));

	if( not antok::generators::functionrgumentHandlerPossibleConst<double>(possible_const_double, function, 0 ) ){
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	if( not antok::generators::functionrgumentHandlerPossibleConst<int>(possible_const_int, function, 0 ) ){
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	std::vector<int*> quantityAddrs_int;
	std::vector<double*> quantityAddrs_double;
	for(unsigned int i = 0; i < 1; ++i) {
		if(not data.insert<int>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs_int.push_back(data.getAddr<int>(quantityNames[i]));

	}
	if (quantityNames.size() == 4) {
		for (unsigned int i = 1; i < quantityNames.size(); ++i) {
			if (not data.insert<double>(quantityNames[i])) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
				return 0;
			}
			quantityAddrs_double.push_back(data.getAddr<double>(quantityNames[i]));
		}
	} else {
		for (unsigned int i = 1; i < quantityNames.size(); ++i) {
			quantityAddrs_double.push_back(nullptr);
		}

	}

	return new antok::user::stefan::functions::CalcCEDARPIDOneL(
	        possible_const_double[0].second,
	        possible_const_double[1].second,
	        new double(),
	        possible_const_double[2].second,
	        possible_const_double[3].second,
	        new double(),
	        possible_const_double[4].second,
	        possible_const_double[5].second,
	        quantityAddrs_int[0],
	        quantityAddrs_double[0],
	        quantityAddrs_double[1],
	        quantityAddrs_double[2]
	        );
}

antok::Function* antok::user::stefan::getCalcAngles3P(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
	if(quantityNames.size() != 4) {
		std::cerr<<"Need 4 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return nullptr;
	}
	using antok::YAMLUtils::hasNodeKey;

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	std::vector<std::pair<std::string, double*> > possible_const;
	args.push_back(std::pair<std::string, std::string>("LVBachelor", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("LVIsoDaughter1", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("LVIsoDaughter2", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("LVBeam", "TLorentzVector"));

	possible_const.push_back(std::pair<std::string, double*>("TargetMass", nullptr));


	if(not antok::generators::functionArgumentHandler(args, function, 0)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	if( not antok::generators::functionrgumentHandlerPossibleConst<double>(possible_const, function, 0 ) ){
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	std::vector<double*> quantityAddrs_double;
	for(unsigned int i = 0; i < 4; ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs_double.push_back(data.getAddr<double>(quantityNames[i]));

	}

	return new antok::user::stefan::functions::CalcAngles3P(
	        data.getAddr<TLorentzVector>(args[0].first),
	        data.getAddr<TLorentzVector>(args[1].first),
	        data.getAddr<TLorentzVector>(args[2].first),
	        data.getAddr<TLorentzVector>(args[3].first),
	        possible_const[0].second,
	        quantityAddrs_double[0],
	        quantityAddrs_double[1],
	        quantityAddrs_double[2],
	        quantityAddrs_double[3]
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
	if(functionName == "determineKaonPionLVLikelihood") {
		antokFunctionPtr = stefan::getDetermineKaonPionLVLikelihood(function, quantityNames, index);
	}
	if(functionName == "calcCEDARPID") {
		antokFunctionPtr = stefan::getCalcCEDARPID(function, quantityNames, index);
	}
	if(functionName == "calcAngles3P") {
		antokFunctionPtr = stefan::getCalcAngles3P(function, quantityNames, index);
	}
	return antokFunctionPtr;
}


namespace {
std::vector<double> getCalcCEDARPIDGetThresholds(const char* name, const YAML::Node& function) {
	using antok::YAMLUtils::hasNodeKey;
	if (not hasNodeKey(function, name)) {
		std::cerr << "Need " << name << " for function \"" << function["Name"] << "\"." << std::endl;
	}

	const YAML::Node& node = function[name];

	try {
		std::vector<double> thresholds = node.as<std::vector<double> >();
		if (thresholds.size() != 9) {
			std::cerr << "Length of thresholds of argument " << name << " is << " <<
			        thresholds.size() << " but should be 9 for function \"" << function["Name"] << "\"." << std::endl;
			return std::vector<double>();
		}
		return thresholds;
	} catch (const YAML::TypedBadConversion<std::vector<double> >& e) {
		std::cerr << "Can not convert argument of node " << name << " for function \"" << function["Name"] << "\"." << std::endl;
		return std::vector<double>();
	}
	return std::vector<double>();
	}
}
