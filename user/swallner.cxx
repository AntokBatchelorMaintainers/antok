
#include<swallner.h>

#include<iostream>
#include<sstream>
#include<string>
#include<functions.hpp>
#include<generators_functions.h>
#include<data.h>
#include<yaml_utils.hpp>

#include<swallner_functions.hpp>

namespace antok{
namespace user{
namespace stefan{

/**
 * Sets the data pointers in the args vector to the address of the variable or to an constant if no variable name, but a number is given
 * @param args Vector of pairs where first: node/variable name, second: data pointet (will be set in this function)
 * @param function: Node of the function
 * @param index: Index of the function call (0 if this arguments have no index)
 * @return true if everything was ok
 */
template< typename T>
bool functionrgumentHandlerPossibleConst( std::vector< std::pair< std::string, T* > >& args,
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

}
}
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

	// complete list of arguments
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichPion", nullptr));
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichKaon", nullptr));
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichProton", nullptr));
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichElectron", nullptr));
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichMuon", nullptr));
    possible_const_per_index.push_back(std::pair<std::string, double* >("PidLRichBackground", nullptr));


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

		if( not antok::user::stefan::functionrgumentHandlerPossibleConst<double>(possible_const, function, 0 ) ){
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

	if(not antok::user::stefan::functionrgumentHandlerPossibleConst<double>(possible_const_per_index, function, index ) ){
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
			                                                data.getAddr<TVector3>(args_per_index[0].first),
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
	if(quantityNames.size() != 5) {
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

	if( not antok::user::stefan::functionrgumentHandlerPossibleConst<double>(possible_const, function, 0) ){
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
																	pid_pion
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
