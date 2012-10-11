#include<generators_cuts.h>

#include<cut.hpp>
#include<data.h>
#include<initializer.h>
#include<object_manager.h>

antok::Cut* antok::generators::generateCut(const YAML::Node& cutEntry) {

	std::string shortName = antok::Initializer::getYAMLStringSafe(cutEntry["ShortName"]);
	std::string longName = antok::Initializer::getYAMLStringSafe(cutEntry["LongName"]);
	std::string abbreviation = antok::Initializer::getYAMLStringSafe(cutEntry["Abbreviation"]);

	if(shortName == "" or longName == "" or abbreviation == "") {
		std::cerr<<"Did not find one of the cut's names (needed are \"ShortName\", \"LongName\" and \"Abbreviation\")."<<std::endl;
		return 0;
	}

	if(not cutEntry["Cut"]) {
		std::cerr<<"Cut \""<<shortName<<"\" does not have required entry \"Cut\"."<<std::endl;
		return 0;
	}

	const YAML::Node& cut = cutEntry["Cut"];

	std::string cutName = antok::Initializer::getYAMLStringSafe(cut["Name"]);
	if(cutName == "") {
		std::cerr<<"Could not get the cut's \"Cut\"->\"Name\" for cut \""<<shortName<<"\"."<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	antok::Cut* antokCut = 0;
	if(cutName == "Range") {

		if(not (cut["Type"] and cut["Variable"])) {
			std::cerr<<"A required entry is missing for \"Range\" cut \""<<shortName<<"\" (either \"Type\" or \"Variable\")"<<std::endl;
			return 0;
		}
		int method = -1;
		double* lowerBound = 0;
		double* upperBound = 0;
		if(cut["LowerBound"] and cut["UpperBound"]) {
			lowerBound = antok::Initializer::getYAMLDoubleAddress(cut["LowerBound"]);
			upperBound = antok::Initializer::getYAMLDoubleAddress(cut["UpperBound"]);
			if(lowerBound == 0 or upperBound == 0) {
				std::cerr<<"Entries \"LowerBound\"/\"UpperBound\" invalid in \"Range\" cut \""<<shortName<<"\", has to be either a variable name or of type double."<<std::endl;
			}
			method = 0;
		} else if(cut["LowerBound"]) {
			lowerBound = antok::Initializer::getYAMLDoubleAddress(cut["LowerBound"]);
			if(lowerBound == 0) {
				std::cerr<<"Entry \"LowerBound\" invalid in \"Range\" cut \""<<shortName<<"\", has to be either a variable name or of type double."<<std::endl;
			}
			method = 4;
		} else if(cut["UpperBound"]) {
			upperBound = antok::Initializer::getYAMLDoubleAddress(cut["UpperBound"]);
			if(upperBound == 0) {
				std::cerr<<"Entries \"LowerBound\"/\"UpperBound\" invalid in \"Range\" cut \""<<shortName<<"\", has to be either a variable name or of type double."<<std::endl;
			}
			method = 2;
		} else {
			std::cerr<<"Either \"LowerBound\" or \"UpperBound\" has to be present in \"Range\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		std::string type = antok::Initializer::getYAMLStringSafe(cut["Type"]);
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
		std::string varName = antok::Initializer::getYAMLStringSafe(cut["Variable"]);
		if(varName == "") {
			std::cerr<<"Could not convert \"Range\" cut \""<<shortName<<"\"'s \"Variable\" entry to std::string."<<std::endl;
			return 0;
		}
		double* variable = data.getDoubleAddr(varName);
		if(variable == 0) {
			std::cerr<<"Could not find \"Range\" cut \""<<shortName<<"\"'s \"Variable\" entry \""<<varName<<"\" in Data."<<std::endl;
			return 0;
		}
		antokCut = new antok::cuts::RangeCut(shortName, longName, abbreviation, lowerBound, upperBound, variable, method);

	} else if (cutName == "Equality") {
		std::cout<<"Equalitycut "<<shortName<<std::endl;
	} else if (cutName == "TriggerMask") {
		std::cerr<<"Trigger mask cut not implemented yet"<<std::endl;
	} else if (cutName == "Group") {
		std::cout<<"Groupcut "<<shortName<<std::endl;
	} else {
		std::cerr<<"Cut \""<<cutName<<"\" not supported."<<std::endl;
		return 0;
	}

	return antokCut;

};

