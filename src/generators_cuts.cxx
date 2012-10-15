#include<generators_cuts.h>

#include<iostream>

#include<TLorentzVector.h>

#include<cut.hpp>
#include<data.h>
#include<initializer.h>
#include<object_manager.h>
#include<yaml_utils.hpp>

namespace {

	template<typename T>
	antok::Cut* _getEqualityCut(const YAML::Node& cut,
	                            const std::string& shortName,
	                            const std::string& longName,
	                            const std::string& abbreviation,
	                            int mode)
	{

		T* value = antok::YAMLUtils::getAddress<T>(cut["Value"]);
		if(value == 0) {
			std::cerr<<"Problem processing \"Value\" entry in \"Equality\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		T* variable = antok::YAMLUtils::getAddress<T>(cut["Variable"]);
		if(variable == 0) {
			std::cerr<<"Problem processing \"Variable\" entry in \"Equality\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		return (new antok::cuts::EqualityCut<T>(shortName, longName, abbreviation, variable, value, mode));

	}

}

antok::Cut* antok::generators::generateRangeCut(const YAML::Node& cut,
                                                const std::string& shortName,
                                                const std::string& longName,
                                                const std::string& abbreviation)
{

	if(not (cut["Type"] and cut["Variable"])) {
		std::cerr<<"A required entry is missing for \"Range\" cut \""<<shortName<<"\" (either \"Type\" or \"Variable\")"<<std::endl;
		return 0;
	}
	int method = -1;
	double* lowerBound = 0;
	double* upperBound = 0;
	if(cut["LowerBound"] and cut["UpperBound"]) {
		lowerBound = antok::YAMLUtils::getAddress<double>(cut["LowerBound"]);
		upperBound = antok::YAMLUtils::getAddress<double>(cut["UpperBound"]);
		if(lowerBound == 0 or upperBound == 0) {
			std::cerr<<"Entries \"LowerBound\"/\"UpperBound\" invalid in \"Range\" cut \""<<shortName<<"\", has to be either a variable name or of type double."<<std::endl;
		}
		method = 0;
	} else if(cut["LowerBound"]) {
		lowerBound = antok::YAMLUtils::getAddress<double>(cut["LowerBound"]);
		if(lowerBound == 0) {
			std::cerr<<"Entry \"LowerBound\" invalid in \"Range\" cut \""<<shortName<<"\", has to be either a variable name or of type double."<<std::endl;
		}
		method = 4;
	} else if(cut["UpperBound"]) {
		upperBound = antok::YAMLUtils::getAddress<double>(cut["UpperBound"]);
		if(upperBound == 0) {
			std::cerr<<"Entries \"LowerBound\"/\"UpperBound\" invalid in \"Range\" cut \""<<shortName<<"\", has to be either a variable name or of type double."<<std::endl;
		}
		method = 2;
	} else {
		std::cerr<<"Either \"LowerBound\" or \"UpperBound\" has to be present in \"Range\" cut \""<<shortName<<"\"."<<std::endl;
		return 0;
	}
	std::string type = antok::YAMLUtils::getString(cut["Type"]);
	if(type == "Inclusive") {
		++method;
	} else if(type == "Exclusive") {
		// nothing to do
	} else if(type == "") {
		std::cerr<<"Could not convert \"Type\" to std::string for \"Range\" cut \""<<shortName<<"\"."<<std::endl;
		return 0;
	} else {
		std::cerr<<"\"Type\" entry in \"Range\" cut \""<<shortName<<"\" has either to be \"Exclusive\" or \"Inclusive\""<<std::endl;
		return 0;
	}
	std::string varName = antok::YAMLUtils::getString(cut["Variable"]);
	if(varName == "") {
		std::cerr<<"Could not convert \"Range\" cut \""<<shortName<<"\"'s \"Variable\" entry to std::string."<<std::endl;
		return 0;
	}
	antok::Data& data = antok::ObjectManager::instance()->getData();
	double* variable = data.getAddr<double>(varName);
	if(variable == 0) {
		std::cerr<<"Could not find \"Range\" cut \""<<shortName<<"\"'s \"Variable\" entry \""<<varName<<"\" in Data."<<std::endl;
		return 0;
	}
	return (new antok::cuts::RangeCut(shortName, longName, abbreviation, lowerBound, upperBound, variable, method));

};

antok::Cut* antok::generators::generateEqualityCut(const YAML::Node& cut,
                                                   const std::string& shortName,
                                                   const std::string& longName,
                                                   const std::string& abbreviation)
{

	if(not (cut["Type"] and cut["Value"] and cut["Variable"])) {
		std::cerr<<"One of the required arguments (\"Type\", \"Value\" and \"Variable\") for \"Equality\" cut \""<<shortName<<"\" is missing."<<std::endl;
		return 0;
	}

	std::string type = antok::YAMLUtils::getString(cut["Type"]);
	if(type == "") {
		std::cerr<<"Could not convert \"Type\" to std::string in \"Equality\" cut \""<<shortName<<"\"."<<std::endl;
		return 0;
	}

	int mode = -1;
	if(type == "==") {
		mode = 0;
	} else if (type == "!=") {
		mode = 1;
	} else {
		std::cerr<<"\"Type\" \""<<type<<"\" not supported by \"Equality\" cut."<<std::endl;
		return 0;
	}

	std::string variableName = antok::YAMLUtils::getString(cut["Variable"]);
	if(variableName == "") {
		std::cerr<<"Could not convert \"Equality\" cut \""<<shortName<<"\"'s \"Variable\" to std::string."<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();
	std::string typeName = data.getType(variableName);

	antok::Cut* antokCut = 0;
	if(typeName == "double") {
		antokCut = _getEqualityCut<double>(cut, shortName, longName, abbreviation, mode);
	} else if (typeName == "int") {
		antokCut = _getEqualityCut<int>(cut, shortName, longName, abbreviation, mode);
	} else if (typeName == "Long64_t") {
		antokCut = _getEqualityCut<Long64_t>(cut, shortName, longName, abbreviation, mode);
	} else if (typeName == "TLorentzVector") {
		antokCut = _getEqualityCut<TLorentzVector>(cut, shortName, longName, abbreviation, mode);
	} else {
		std::cerr<<"Type \""<<typeName<<"\" not supported in \"Equality\" cut \""<<shortName<<"\"."<<std::endl;
		return 0;
	}

	return antokCut;

};

