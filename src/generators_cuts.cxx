#include<generators_cuts.h>

#include<cut.hpp>
#include<data.h>
#include<initializer.h>
#include<object_manager.h>
#include<yaml_utils.hpp>

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

	return 0;
	return (new antok::cuts::EqualityCut<double>(shortName, longName, abbreviation, 0, 0, mode));

};

